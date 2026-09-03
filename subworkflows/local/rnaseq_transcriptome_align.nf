// RNA-seq bowtie2 transcriptome alignment → qpass → dedup. Parallel to
// transcriptome_align.nf (ribo-seq), but does NOT call RIBOPY_CREATE — it emits
// the final merged/dedup BED (ribo_bed) for the parent workflow to feed into
// RIBOPY_RNASEQ_SET. Shared modules are disambiguated via
// .*:RNASEQ_TRANSCRIPTOME_ALIGN:MODULE selectors in conf/modules.config.

include { BOWTIE2_TRANSCRIPTOME_RNASEQ }                    from '../../modules/local/bowtie2_transcriptome_rnaseq.nf'
include { SAMTOOLS_QPASS }                                  from '../../modules/local/samtools_qpass.nf'
include { SAMTOOLS_MERGE }                                  from '../../modules/local/samtools_merge.nf'
include { BAM_TO_BED as RNASEQ_TX_QPASS_BAM_TO_BED }        from '../../modules/local/bam_to_bed.nf'
include { BAM_TO_BED as RNASEQ_TX_MERGED_QPASS_BED }        from '../../modules/local/bam_to_bed.nf'
include { BAM_TO_BED as RNASEQ_TX_MERGED_DEDUP_BED }        from '../../modules/local/bam_to_bed.nf'
include { CONCAT_SORT_BED as RNASEQ_MERGE_TX_PRE_DEDUP_BED } from '../../modules/local/concat_sort_bed.nf'
include { RFC_DEDUP }                                       from '../../modules/local/rfc_dedup.nf'
include { SEPARATE_BED }                                    from '../../modules/local/separate_bed.nf'
include { RFC_EXTRACT_DEDUP_READS }                         from '../../modules/local/rfc_extract_dedup_reads.nf'
include { UMICOLLAPSE_DEDUP }                               from '../../modules/local/umicollapse_dedup.nf'
include { SPLIT_DEDUP_BAM }                                 from '../../modules/local/split_dedup_bam.nf'

workflow RNASEQ_TRANSCRIPTOME_ALIGN {
    take:
    ch_reads        // [ meta(id,lane,strand,single_end), reads ]
    ch_tx_index     // value: [index_base, index_files]

    main:
    def dedup = Utils.resolve_rnaseq_dedup_method(params)

    BOWTIE2_TRANSCRIPTOME_RNASEQ(ch_reads, ch_tx_index)
    ch_tx_bam = BOWTIE2_TRANSCRIPTOME_RNASEQ.out.bam

    SAMTOOLS_QPASS(ch_tx_bam)
    ch_qpass_bam   = SAMTOOLS_QPASS.out.bam
    ch_qpass_total = SAMTOOLS_QPASS.out.total_count

    ch_lane_meta    = ch_tx_bam.map { meta, bam -> [meta.id, meta] }
    ch_sample_lanes = ch_lane_meta.groupTuple()

    ch_merge_in = ch_qpass_bam
        .map { meta, bam, bai -> [meta.id, bam] }
        .groupTuple()
        .map { id, bams -> [[id: id, strand: 'F'], bams] }
    SAMTOOLS_MERGE(ch_merge_in)
    ch_merged_qpass_bam = SAMTOOLS_MERGE.out.bam

    ch_ribo_bed             = Channel.empty()
    ch_individual_dedup_cnt = Channel.empty()

    if (dedup == 'position') {
        RNASEQ_TX_QPASS_BAM_TO_BED(ch_qpass_bam.map { meta, bam, bai -> [meta, bam] })
        RNASEQ_MERGE_TX_PRE_DEDUP_BED(
            RNASEQ_TX_QPASS_BAM_TO_BED.out.bed.map { meta, bed -> [meta.id, bed] }.groupTuple()
                .map { id, beds -> [[id: id, strand: 'F'], beds] }
        )
        RFC_DEDUP(RNASEQ_MERGE_TX_PRE_DEDUP_BED.out.bed)
        ch_ribo_bed = RFC_DEDUP.out.bed

        SEPARATE_BED(
            RFC_DEDUP.out.bed.map { smeta, bed -> [smeta.id, smeta, bed] }
                .join(ch_sample_lanes)
                .map { id, smeta, bed, metas -> [smeta, bed, Utils.lane_ids(metas)] }
        )
        ch_individual_dedup_cnt = SEPARATE_BED.out.total_counts.map { smeta, f -> [smeta.id, f] }
            .join(ch_sample_lanes)
            .flatMap { id, files, metas -> Utils.lane_pairs(files, metas) }

        ch_extract_in = RFC_DEDUP.out.bed
            .map { smeta, bed -> [smeta.id, bed] }
            .join(ch_merged_qpass_bam.map { smeta, bam, bai -> [smeta.id, bam] })
            .map { id, bed, bam -> [[id: id, strand: 'F'], bed, bam] }
        RFC_EXTRACT_DEDUP_READS(ch_extract_in)
    }
    else if (dedup == 'umicollapse') {
        UMICOLLAPSE_DEDUP(ch_merged_qpass_bam)

        SPLIT_DEDUP_BAM(
            UMICOLLAPSE_DEDUP.out.bam.map { smeta, bam, bai -> [smeta.id, smeta, bam, bai] }
                .join(ch_sample_lanes)
                .map { id, smeta, bam, bai, metas -> [smeta, bam, bai, Utils.lane_ids(metas)] }
        )
        ch_individual_dedup_cnt = SPLIT_DEDUP_BAM.out.total_counts.map { smeta, f -> [smeta.id, f] }
            .join(ch_sample_lanes)
            .flatMap { id, files, metas -> Utils.lane_pairs(files, metas) }

        RNASEQ_TX_MERGED_DEDUP_BED(UMICOLLAPSE_DEDUP.out.bam.map { smeta, bam, bai -> [smeta, bam] })
        ch_ribo_bed = RNASEQ_TX_MERGED_DEDUP_BED.out.bed
    }
    else {
        RNASEQ_TX_MERGED_QPASS_BED(ch_merged_qpass_bam.map { smeta, bam, bai -> [smeta, bam] })
        ch_ribo_bed = RNASEQ_TX_MERGED_QPASS_BED.out.bed
        ch_individual_dedup_cnt = ch_qpass_total
    }

    emit:
    bowtie2_log             = BOWTIE2_TRANSCRIPTOME_RNASEQ.out.log
    qpass_total_count       = ch_qpass_total
    individual_dedup_counts = ch_individual_dedup_cnt
    ribo_bed                = ch_ribo_bed   // [ smeta, bed ] per sample
}
