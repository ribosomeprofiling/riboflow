// Top-level genome (ribo-seq) workflow. Builds the input channels, runs
// preprocess → genome alignment (+ optional transcriptome dedup) → stats.

include { PREPROCESS }                from '../subworkflows/local/preprocess.nf'
include { GENOME_ALIGN }              from '../subworkflows/local/genome_align.nf'
include { STAR_TRANSCRIPTOME_DEDUP }  from '../subworkflows/local/star_transcriptome_dedup.nf'
include { TRANSCRIPTOME_ALIGN }       from '../subworkflows/local/transcriptome_align.nf'
include { ALIGNMENT_STATS }           from '../subworkflows/local/alignment_stats.nf'
include { TRANSCRIPTOME_STATS }       from '../subworkflows/local/transcriptome_stats.nf'
include { WRITE_CORRESPONDENCE }      from '../modules/local/write_correspondence.nf'

include { RNASEQ_PREPROCESS }                              from '../subworkflows/local/rnaseq_preprocess.nf'
include { RNASEQ_GENOME_ALIGN }                            from '../subworkflows/local/rnaseq_genome_align.nf'
include { RNASEQ_TRANSCRIPTOME_ALIGN }                     from '../subworkflows/local/rnaseq_transcriptome_align.nf'
include { ALIGNMENT_STATS     as RNASEQ_GENOME_STATS }     from '../subworkflows/local/alignment_stats.nf'
include { TRANSCRIPTOME_STATS as RNASEQ_TX_STATS }         from '../subworkflows/local/transcriptome_stats.nf'
include { RIBOPY_RNASEQ_SET }                              from '../modules/local/ribopy_rnaseq_set.nf'
include { RIBOPY_MERGE }                                   from '../modules/local/ribopy_merge.nf'
include { STAR_INDEX }                                     from '../modules/local/star_index.nf'

workflow RIBOFLOW {

    main:
    // ── Required inputs ───────────────────────────────────────────────────
    if (!(params.input instanceof Map) || !(params.input.fastq instanceof Map) || params.input.fastq.isEmpty()) {
        error "input.fastq is required: a map of sample → list of FASTQ paths (see example_*.yaml)."
    }
    if (!(params.input.reference instanceof Map)) {
        error "input.reference is required (filter index, and genome index / FASTA+GTF and/or transcriptome index)."
    }

    // ── Build per-lane input channel ──────────────────────────────────────
    // Ribo-seq lanes are single-end: one FASTQ per lane. (RiboFlow.groovy:171-177)
    // NOTE: `lane` is the 1-based POSITION in the YAML list, and every storeDir
    // keys on <sample>.<lane>. Re-ordering or inserting FASTQs for a sample
    // between runs re-points an existing lane at a different file while stored
    // artifacts are reused — clear the intermediates when you change a sample's
    // FASTQ list. index_fastq_correspondence.txt records the mapping of each run.
    def fastq_base = (params.input?.fastq_base ?: '')
    if (fastq_base && !fastq_base.endsWith('/')) fastq_base = "${fastq_base}/"

    def input_list = []
    params.input.fastq.each { sample, lanes ->
        lanes.eachWithIndex { fq, i ->
            input_list << [[id: sample, lane: (i + 1), strand: 'F'], file("${fastq_base}${fq}")]
        }
    }
    ch_reads = Channel.fromList(input_list)

    // ── Validate path selection ────────────────────────────────────────────
    def do_genome = Utils.as_bool(params.genome?.run, true)
    def do_tx     = Utils.as_bool(params.transcriptome?.run, false)
    if (!do_genome && !do_tx) {
        error "At least one alignment path must be enabled: set genome.run=true and/or transcriptome.run=true."
    }

    def do_rnaseq     = Utils.as_bool(params.do_rnaseq, false) && (params.rnaseq?.fastq != null)
    def do_rna_genome = do_rnaseq && do_genome
    def do_rna_tx     = do_rnaseq && do_tx

    // ── Validate the per-route quality filter ──────────────────────────────
    // Fails here, before any process is submitted, instead of inside a task.
    def filter_errors = Utils.validate_routes(params)
    if (filter_errors) {
        error "Invalid post-alignment filter configuration:\n  - " + filter_errors.join("\n  - ")
    }
    def active_routes = (do_genome ? ['genome'] : []) + (do_tx ? ['transcriptome'] : []) +
                        (do_rna_genome ? ['rnaseq_genome'] : []) + (do_rna_tx ? ['rnaseq_transcriptome'] : [])
    Utils.route_notes(params, active_routes).each { log.warn(it) }

    // ── Paired-end RNA-seq validation ──────────────────────────────────────
    // Must run before any subworkflow is wired: an error() raised after channels
    // exist races the stats stage's own "nothing produced" error and gets masked.
    // The genome path supports paired-end for dedup methods `none` and `umicollapse`.
    // Still unsupported: PE + `position` (needs fragment/BEDPE dedup) and the
    // RNA-seq transcriptome path under PE.
    if (do_rnaseq) {
        def has_pe = params.rnaseq.fastq.any { sample, lanes -> lanes.any { it instanceof List } }
        if (has_pe) {
            if (do_rna_tx) {
                error "Paired-end RNA-seq is not yet supported for the transcriptome path. Use single-end, or disable the RNA-seq transcriptome path."
            }
            if (Utils.resolve_rnaseq_dedup_method(params) == 'position') {
                error "PE RNA-seq with dedup_method=position is not yet supported. Use umicollapse or none."
            }
        }
    }

    // ── Required tool arguments (would otherwise interpolate `null` into a command)
    if (params.clip_arguments == null || params.clip_arguments.toString().trim().isEmpty()) {
        error "clip_arguments is required (cutadapt adapter/trim options for the ribo-seq reads)."
    }
    if (params.alignment_arguments?.filter == null) {
        error "alignment_arguments.filter is required (bowtie2 options for the rRNA/tRNA filter)."
    }
    if (do_genome && params.star?.ribo_arguments == null) {
        error "star.ribo_arguments is missing. If you passed `--star.<key>` on the command line, note that on Nextflow 26.x a dotted CLI param replaces the whole `star` map — set it in the params YAML instead."
    }
    if (do_rnaseq && params.star?.rnaseq_arguments == null) {
        error "star.rnaseq_arguments is missing (see the note above about dotted CLI params)."
    }
    if (do_rnaseq && (params.rnaseq?.clip_arguments == null || params.rnaseq.clip_arguments.toString().trim().isEmpty())) {
        error "rnaseq.clip_arguments is required when do_rnaseq is true."
    }

    // Detect build-from-FASTA mode vs pre-built index mode.
    def genome_fasta = params.input?.reference?.genome_fasta ?: null
    def genome_gtf   = params.input?.reference?.gtf          ?: null

    // ── Optional input existence checks (RiboFlow.groovy:200-224) ──────────
    if (Utils.as_bool(params.do_check_file_existence, false)) {
        if (do_genome) {
            if (!genome_fasta) {
                ['SA', 'SAindex', 'Genome', 'chrNameLength.txt'].each { f ->
                    assert file("${params.input.reference.genome}/${f}").exists() :
                        "Missing STAR index file: ${params.input.reference.genome}/${f}"
                }
            } else {
                assert file(genome_fasta).exists() : "Missing genome FASTA: ${genome_fasta}"
                assert file(genome_gtf  ).exists() : "Missing genome GTF: ${genome_gtf}"
            }
        }
        def filter_pref = params.input.reference.filter.replaceAll('\\*', '')
        ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2'].each { s ->
            assert file("${filter_pref}.${s}").exists() :
                "Missing bowtie2 filter index file: ${filter_pref}.${s}"
        }
        input_list.each { meta, f ->
            assert f.exists() : "Missing input FASTQ: ${f}"
        }
    }

    // Sample/lane → FASTQ correspondence table.
    ch_corr = ch_reads
        .map { meta, fq -> "${meta.id}\t${meta.lane}\t${fq}" }
        .collectFile(name: 'correspondence.txt', newLine: true)
    WRITE_CORRESPONDENCE(ch_corr)

    // ── Reference channels ─────────────────────────────────────────────────
    def filter_glob = params.input.reference.filter
    def filter_base = filter_glob.split('/')[-1].replaceAll('\\*$', '').replaceAll('\\.$', '')
    ch_filter_index = Channel.value([filter_base, files(filter_glob)])

    if (do_genome) {
        if (genome_fasta) {
            if (!genome_gtf) {
                error "input.reference.genome_fasta is set but input.reference.gtf is missing. Both are required to build a STAR index."
            }
            STAR_INDEX(
                Channel.value(file(genome_fasta)),
                Channel.value(file(genome_gtf))
            )
            ch_genome_index = STAR_INDEX.out.index.first()
        } else {
            ch_genome_index = Channel.value(file(params.input.reference.genome))
        }
    }

    if (do_tx) {
        def tx_glob = params.input.reference.transcriptome
        def tx_base = tx_glob.split('/')[-1].replaceAll('\\*$', '').replaceAll('\\.$', '')
        ch_tx_index  = Channel.value([tx_base, files(tx_glob)])
        ch_regions   = Channel.value(file(params.input.reference.regions))
        ch_lengths   = Channel.value(file(params.input.reference.transcript_lengths))
    }

    // Per-sample expmeta channel for RIBOPY_CREATE (optional).
    // Set ribo.metadata.files.<sample>: path in your YAML to embed per-sample metadata.
    def meta_base_raw = (params.ribo?.metadata?.base ?: '')
    def meta_base_pfx = meta_base_raw && !meta_base_raw.endsWith('/') ? "${meta_base_raw}/" : meta_base_raw
    ch_meta_files = params.ribo?.metadata?.files
        ? Channel.fromList(params.ribo.metadata.files.collect { s, f -> [ s, file("${meta_base_pfx}${f}") ] })
        : Channel.empty()

    // ── Pipeline ───────────────────────────────────────────────────────────
    PREPROCESS(ch_reads, ch_filter_index)

    if (do_genome) {
        GENOME_ALIGN(PREPROCESS.out.reads_for_genome, ch_genome_index)

        if (Utils.do_tx_dedup(params)) {
            STAR_TRANSCRIPTOME_DEDUP(GENOME_ALIGN.out.transcriptome_bam)
        }

        ALIGNMENT_STATS(
            PREPROCESS.out.clip_log,
            PREPROCESS.out.filter_log,
            GENOME_ALIGN.out.genome_log,
            GENOME_ALIGN.out.secondary_count,
            GENOME_ALIGN.out.qpass_total_count,
            GENOME_ALIGN.out.qpass_primary_count,
            GENOME_ALIGN.out.qpass_secondary_count,
            GENOME_ALIGN.out.qpass_unique_count,
            GENOME_ALIGN.out.individual_dedup_counts,
            GENOME_ALIGN.out.merged_dedup_counts,
        )
    }

    if (do_tx) {
        TRANSCRIPTOME_ALIGN(
            PREPROCESS.out.reads_for_genome,
            ch_tx_index, ch_regions, ch_lengths,
            ch_meta_files
        )

        TRANSCRIPTOME_STATS(
            PREPROCESS.out.clip_log,
            PREPROCESS.out.filter_log,
            TRANSCRIPTOME_ALIGN.out.bowtie2_log,
            TRANSCRIPTOME_ALIGN.out.qpass_total_count,
            TRANSCRIPTOME_ALIGN.out.individual_dedup_counts,
        )
    }

    // ── RNA-seq path ───────────────────────────────────────────────────────────
    if (do_rnaseq) {
        // Build RNA-seq input channel — SE entries are strings, PE entries are 2-element lists.
        def rna_base = (params.rnaseq?.fastq_base ?: '')
        if (rna_base && !rna_base.endsWith('/')) rna_base = "${rna_base}/"
        def rna_list = []
        params.rnaseq.fastq.each { sample, lanes ->
            lanes.eachWithIndex { lane_entry, i ->
                def is_pe = (lane_entry instanceof List)
                def reads = is_pe
                    ? [file("${rna_base}${lane_entry[0]}"), file("${rna_base}${lane_entry[1]}")]
                    : file("${rna_base}${lane_entry}")
                rna_list << [[id: sample, lane: (i + 1), strand: 'F', single_end: !is_pe], reads]
            }
        }
        ch_rna_reads = Channel.fromList(rna_list)

        RNASEQ_PREPROCESS(ch_rna_reads, ch_filter_index)

        if (do_rna_genome) {
            RNASEQ_GENOME_ALIGN(RNASEQ_PREPROCESS.out.reads_for_genome, ch_genome_index)
            RNASEQ_GENOME_STATS(
                RNASEQ_PREPROCESS.out.clip_log,
                RNASEQ_PREPROCESS.out.filter_log,
                RNASEQ_GENOME_ALIGN.out.genome_log,
                RNASEQ_GENOME_ALIGN.out.secondary_count,
                RNASEQ_GENOME_ALIGN.out.qpass_total_count,
                RNASEQ_GENOME_ALIGN.out.qpass_primary_count,
                RNASEQ_GENOME_ALIGN.out.qpass_secondary_count,
                RNASEQ_GENOME_ALIGN.out.qpass_unique_count,
                RNASEQ_GENOME_ALIGN.out.individual_dedup_counts,
                RNASEQ_GENOME_ALIGN.out.merged_dedup_counts,
            )
        }

        if (do_rna_tx) {
            RNASEQ_TRANSCRIPTOME_ALIGN(RNASEQ_PREPROCESS.out.reads_for_genome, ch_tx_index)
            RNASEQ_TX_STATS(
                RNASEQ_PREPROCESS.out.clip_log,
                RNASEQ_PREPROCESS.out.filter_log,
                RNASEQ_TRANSCRIPTOME_ALIGN.out.bowtie2_log,
                RNASEQ_TRANSCRIPTOME_ALIGN.out.qpass_total_count,
                RNASEQ_TRANSCRIPTOME_ALIGN.out.individual_dedup_counts,
            )

            // Merge the RNA-seq transcriptome BED into ribo-seq .ribo files
            // (only when the ribo-seq transcriptome path also ran).
            if (do_tx) {
                ch_ribopy_set_in = TRANSCRIPTOME_ALIGN.out.ribo
                    .map { meta, ribo -> [meta.id, ribo] }
                    .join(RNASEQ_TRANSCRIPTOME_ALIGN.out.ribo_bed.map { meta, bed -> [meta.id, bed] })
                    .map { id, ribo, bed -> [[id: id], ribo, bed] }
                RIBOPY_RNASEQ_SET(ch_ribopy_set_in)
            }
        }
    }

    // ── Merge all per-sample .ribo files → all.ribo ────────────────────────
    // Uses the most up-to-date .ribo per sample: RIBOPY_RNASEQ_SET output
    // (RNA-seq embedded) takes precedence; samples without RNA-seq fall back
    // to RIBOPY_CREATE output.
    if (do_tx) {
        if (do_rna_tx) {
            ch_ribo_base    = TRANSCRIPTOME_ALIGN.out.ribo.map { meta, ribo -> [meta.id, ribo] }
            ch_ribo_updated = RIBOPY_RNASEQ_SET.out.ribo.map  { meta, ribo -> [meta.id, ribo] }
            ch_final_ribo   = ch_ribo_base
                .join(ch_ribo_updated, remainder: true)
                .map { id, base, updated -> updated ?: base }
        } else {
            ch_final_ribo = TRANSCRIPTOME_ALIGN.out.ribo.map { meta, ribo -> ribo }
        }
        // ribopy merge requires >=2 inputs; skip it for single-sample runs.
        RIBOPY_MERGE(ch_final_ribo.collect().filter { it.size() > 1 })
    }
}
