// OPTIONAL: quality-filter, merge and dedup the STAR transcriptome-projected BAM
// (star.output_transcriptome_bam), emitting BAMs plus a per-lane distinct-read
// count for the stats (`transcriptome_projection`). No .ribo.
//
// Every BAM here keeps STAR's read-name grouping and carries no index
// (ext.unsorted on the processes): Salmon needs a read's alignments consecutive,
// and coordinate sorting the read×isoform BAM is both wrong for it and very slow.
//
// Shared by the ribo-seq and RNA-seq genome paths: included as
// STAR_TRANSCRIPTOME_DEDUP and RNASEQ_STAR_TRANSCRIPTOME_DEDUP, each with its
// own tree in conf/modules.config (`star_transcriptome/`, `rnaseq/star_transcriptome/`).
//
// Dedup follows the genome BAM of the same path:
//   none         merged qpass BAM is final
//   position     rfc dedup on the projected BEDs (ribo-seq only)
//   umicollapse  read-name subset of the genome post-dedup BAM (no second UMICollapse)

include { SAMTOOLS_QPASS }            from '../../modules/local/samtools_qpass.nf'
include { BAM_TO_BED }                from '../../modules/local/bam_to_bed.nf'
include { SAMTOOLS_MERGE }            from '../../modules/local/samtools_merge.nf'
include { CONCAT_SORT_BED }           from '../../modules/local/concat_sort_bed.nf'
include { RFC_DEDUP }                 from '../../modules/local/rfc_dedup.nf'
include { SEPARATE_BED }              from '../../modules/local/separate_bed.nf'
include { RFC_EXTRACT_DEDUP_READS }   from '../../modules/local/rfc_extract_dedup_reads.nf'
include { SUBSET_BAM_BY_GENOME_DEDUP } from '../../modules/local/subset_bam_by_genome_dedup.nf'
include { SPLIT_DEDUP_BAM }           from '../../modules/local/split_dedup_bam.nf'
include { COUNT_DISTINCT_READS }      from '../../modules/local/count_distinct_reads.nf'

workflow STAR_TRANSCRIPTOME_DEDUP {
    take:
    ch_tx_bam           // [ meta(id,lane,strand[,single_end]), transcriptome_bam ]
    ch_genome_dedup_bam // [ smeta, bam, bai ] merged post-dedup genome BAM (umicollapse branch)
    dedup               // 'none' | 'position' | 'umicollapse' — the genome dedup of this path

    main:
    // ext.unsorted/prefix set in conf/modules.config via the fully-qualified selector.
    SAMTOOLS_QPASS(ch_tx_bam)
    ch_qpass_bam = SAMTOOLS_QPASS.out.bam   // [ meta, bam, [] ]  (bai slot empty)

    ch_lane_meta    = ch_tx_bam.map { meta, bam -> [meta.id, meta] }
    ch_sample_lanes = ch_lane_meta.groupTuple()

    ch_merge_in = ch_qpass_bam
        .map { meta, bam, bai -> [meta.id, bam] }
        .groupTuple()
        .map { id, bams -> [[id: id, strand: 'F'], bams] }
    SAMTOOLS_MERGE(ch_merge_in)
    ch_merged_qpass_bam = SAMTOOLS_MERGE.out.bam

    ch_final_bam = Channel.empty()   // [ smeta, bam, [] ]
    ch_lane_bam  = Channel.empty()   // [ meta, bam ] per lane, post-dedup (none: the qpass BAMs)
    ch_split_in  = Channel.empty()   // [ smeta, bam, bai ] merged post-dedup BAM to fan out per lane

    if (dedup == 'none') {
        ch_final_bam = ch_merged_qpass_bam
        ch_lane_bam  = ch_qpass_bam.map { meta, bam, bai -> [meta, bam] }
    }
    else if (dedup == 'position') {
        BAM_TO_BED(ch_qpass_bam.map { meta, bam, bai -> [meta, bam] })   // append_index via modules.config
        CONCAT_SORT_BED(
            BAM_TO_BED.out.bed.map { meta, bed -> [meta.id, bed] }.groupTuple()
                .map { id, beds -> [[id: id, strand: 'F'], beds] }
        )
        RFC_DEDUP(CONCAT_SORT_BED.out.bed)

        SEPARATE_BED(
            RFC_DEDUP.out.bed.map { smeta, bed -> [smeta.id, smeta, bed] }
                .join(ch_sample_lanes)
                .map { id, smeta, bed, metas -> [smeta, bed, Utils.lane_ids(metas)] }
        )

        ch_extract_in = RFC_DEDUP.out.bed
            .map { smeta, bed -> [smeta.id, bed] }
            .join(ch_merged_qpass_bam.map { smeta, bam, bai -> [smeta.id, bam] })
            .map { id, bed, bam -> [[id: id, strand: 'F'], bed, bam] }
        RFC_EXTRACT_DEDUP_READS(ch_extract_in)
        ch_final_bam = RFC_EXTRACT_DEDUP_READS.out.bam
        ch_split_in  = ch_final_bam
    }
    else if (dedup == 'umicollapse') {
        // UMIs were already collapsed on the genome BAM; keep the same reads here
        // instead of running UMICollapse again on the much larger read×isoform BAM.
        ch_subset_in = ch_merged_qpass_bam
            .map { smeta, bam, bai -> [smeta.id, smeta, bam] }
            .join(ch_genome_dedup_bam.map { smeta, bam, bai -> [smeta.id, bam] })
            .map { id, smeta, tx_bam, g_bam -> [smeta, tx_bam, g_bam] }
        SUBSET_BAM_BY_GENOME_DEDUP(ch_subset_in)
        ch_final_bam = SUBSET_BAM_BY_GENOME_DEDUP.out.bam.map { smeta, bam -> [smeta, bam, []] }
        ch_split_in  = ch_final_bam
    }

    // Merged post-dedup BAM → per-lane BAMs (one samtools split per sample; empty
    // channel under dedup none, where the per-lane qpass BAMs are already final).
    SPLIT_DEDUP_BAM(
        ch_split_in.map { smeta, bam, bai -> [smeta.id, smeta, bam, bai] }
            .join(ch_sample_lanes)
            .map { id, smeta, bam, bai, metas -> [smeta, bam, bai, Utils.lane_ids(metas)] }
    )
    ch_lane_bam = ch_lane_bam.mix(
        SPLIT_DEDUP_BAM.out.bams
            .map { smeta, bams, bais -> [smeta.id, bams] }
            .join(ch_sample_lanes)
            .flatMap { id, files, metas -> Utils.lane_pairs(files, metas) }
    )

    // Distinct reads per lane → stats row `transcriptome_projection`. Lanes
    // partition the merged BAM, so the per-sample sum equals the merged count.
    COUNT_DISTINCT_READS(ch_lane_bam)

    emit:
    bam              = ch_final_bam
    lane_read_counts = COUNT_DISTINCT_READS.out.count   // [ meta, count ]
}
