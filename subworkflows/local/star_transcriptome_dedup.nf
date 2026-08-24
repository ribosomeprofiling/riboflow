// OPTIONAL: dedup the STAR transcriptome-projected BAM through the same
// position/umicollapse chain as the genome path, emitting BAM/BED only (no
// .ribo, no stats). Runs only when do_tx_dedup. (RiboFlow.groovy:675-975)
// Writes to its own `star_transcriptome` tree (see conf/modules.config).

include { SAMTOOLS_QPASS }            from '../../modules/local/samtools_qpass.nf'
include { BAM_TO_BED }                from '../../modules/local/bam_to_bed.nf'
include { SAMTOOLS_MERGE }            from '../../modules/local/samtools_merge.nf'
include { CONCAT_SORT_BED }           from '../../modules/local/concat_sort_bed.nf'
include { RFC_DEDUP }                 from '../../modules/local/rfc_dedup.nf'
include { SEPARATE_BED }              from '../../modules/local/separate_bed.nf'
include { RFC_EXTRACT_DEDUP_READS }   from '../../modules/local/rfc_extract_dedup_reads.nf'
include { UMICOLLAPSE_DEDUP }         from '../../modules/local/umicollapse_dedup.nf'
include { SPLIT_DEDUP_BAM }           from '../../modules/local/split_dedup_bam.nf'

workflow STAR_TRANSCRIPTOME_DEDUP {
    take:
    ch_tx_bam   // [ meta(id,lane,strand), transcriptome_bam ]

    main:
    def dedup = Utils.resolve_dedup_method(params)

    // qpass with pre-sort (STAR tx BAM is QNAME-sorted). ext.presort/prefix set
    // in conf/modules.config via the fully-qualified selector.
    SAMTOOLS_QPASS(ch_tx_bam)
    ch_qpass_bam = SAMTOOLS_QPASS.out.bam   // [ meta, bam, bai ]

    ch_lane_meta    = ch_tx_bam.map { meta, bam -> [meta.id, meta] }
    ch_sample_lanes = ch_lane_meta.groupTuple()

    ch_merge_in = ch_qpass_bam
        .map { meta, bam, bai -> [meta.id, bam] }
        .groupTuple()
        .map { id, bams -> [[id: id, strand: 'F'], bams] }
    SAMTOOLS_MERGE(ch_merge_in)
    ch_merged_qpass_bam = SAMTOOLS_MERGE.out.bam

    ch_final_bam = Channel.empty()

    if (dedup == 'position') {
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
    }
    else if (dedup == 'umicollapse') {
        UMICOLLAPSE_DEDUP(ch_merged_qpass_bam)
        SPLIT_DEDUP_BAM(
            UMICOLLAPSE_DEDUP.out.bam.map { smeta, bam, bai -> [smeta.id, smeta, bam, bai] }
                .join(ch_sample_lanes)
                .map { id, smeta, bam, bai, metas -> [smeta, bam, bai, Utils.lane_ids(metas)] }
        )
        ch_final_bam = UMICOLLAPSE_DEDUP.out.bam
    }

    emit:
    bam = ch_final_bam
}
