// Genome alignment stats: per-lane rows → combined essential + per-sample sums →
// published stats.csv / individual_stats.csv. (RiboFlow.groovy:1486-1775)
//
// One subworkflow serves both genome routes: included as ALIGNMENT_STATS (ribo-seq)
// and as RNASEQ_GENOME_STATS (RNA-seq); conf/modules.config keys the per-process
// ext.route / storeDir on those two names.
//
// Every join is strict (failOnMismatch / failOnDuplicate) and the collects refuse
// an empty set: a lane missing from any input used to be dropped silently, and a
// run could finish "successfully" with no stats.csv at all.

include { STATS_INDIVIDUAL }                        from '../../modules/local/stats_individual.nf'
include { STATS_COMBINE as COMBINE_INDIVIDUAL }     from '../../modules/local/stats_combine.nf'
include { STATS_COMBINE as COMBINE_MERGED }         from '../../modules/local/stats_combine.nf'
include { STATS_SUM }                               from '../../modules/local/stats_sum.nf'
include { STATS_PUBLISH }                           from '../../modules/local/stats_publish.nf'

workflow ALIGNMENT_STATS {
    take:
    ch_clip_log                // [ meta, log ]
    ch_filter_log              // [ meta, log ]
    ch_genome_log              // [ meta, log ]
    ch_secondary_count         // [ meta, count ]
    ch_qpass_total_count       // [ meta, total ]
    ch_qpass_primary_count     // [ meta, primary ]
    ch_qpass_secondary_count   // [ meta, secondary ]
    ch_qpass_unique_count      // [ meta, unique ]
    ch_individual_dedup_counts // [ meta, t, p, s, u ]
    ch_merged_dedup_counts     // [ smeta, t, p, s, u ] (EMPTY when dedup none)

    main:
    def placeholder = file("${projectDir}/assets/NO_FILE.gz")

    // Per-lane stats input: join everything on the per-lane meta.
    ch_stats_in = ch_clip_log
        .join(ch_filter_log,              failOnMismatch: true, failOnDuplicate: true)
        .join(ch_genome_log,              failOnMismatch: true, failOnDuplicate: true)
        .join(ch_secondary_count,         failOnMismatch: true, failOnDuplicate: true)
        .join(ch_qpass_total_count,       failOnMismatch: true, failOnDuplicate: true)
        .join(ch_qpass_primary_count,     failOnMismatch: true, failOnDuplicate: true)
        .join(ch_qpass_secondary_count,   failOnMismatch: true, failOnDuplicate: true)
        .join(ch_individual_dedup_counts, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_qpass_unique_count,      failOnMismatch: true, failOnDuplicate: true)
    STATS_INDIVIDUAL(ch_stats_in)

    // Combined per-lane essential CSV.
    COMBINE_INDIVIDUAL(
        STATS_INDIVIDUAL.out.csv.map { meta, csv -> csv }.collect()
            .ifEmpty { error "No per-lane genome stats were produced (an upstream process failed, or every lane was dropped)." }
    )

    // Per-sample sums. The merged dedup counts exist only under position/umicollapse
    // (the channel is empty for none), so a remainder-join maps "no counts" onto the
    // placeholder that STATS_SUM ignores when ext.use_merged_counts is false.
    ch_grouped = STATS_INDIVIDUAL.out.csv
        .map { meta, csv -> [meta.id, csv] }
        .groupTuple()
        .map { id, csvs -> [[id: id, strand: 'F'], csvs] }

    // An unmatched left row comes back as [id, smeta, csvs, null] (one null, not four).
    ch_sum_in = ch_grouped
        .map { smeta, csvs -> [smeta.id, smeta, csvs] }
        .join(ch_merged_dedup_counts.map { smeta, t, p, s, u -> [smeta.id, t, p, s, u] }, remainder: true)
        .filter { row -> row[1] != null }
        .map { row ->
            def counts = (row.size() >= 7) ? row[3..6] : [null, null, null, null]
            [row[1], row[2]] + counts.collect { it ?: placeholder }
        }
    STATS_SUM(ch_sum_in)

    COMBINE_MERGED(
        STATS_SUM.out.csv.map { meta, csv -> csv }.collect()
            .ifEmpty { error "No per-sample genome stats were produced (an upstream process failed, or every lane was dropped)." }
    )

    STATS_PUBLISH(COMBINE_INDIVIDUAL.out.csv, COMBINE_MERGED.out.csv)

    emit:
    individual = STATS_PUBLISH.out.individual
    merged     = STATS_PUBLISH.out.merged
}
