// RNA-seq preprocess: clip → [UMI extract (umicollapse only; SE or PE)] → rRNA/tRNA
// filter. Emits the genome/transcriptome-ready unaligned reads (single channel,
// SE single file or PE [R1,R2] list) plus the clip/filter logs the stats stage
// consumes. Parallel to preprocess.nf (ribo-seq).

include { CUTADAPT_CLIP_RNASEQ }        from '../../modules/local/cutadapt_clip_rnaseq.nf'
include { UMITOOLS_EXTRACT }            from '../../modules/local/umitools_extract.nf'
include { BOWTIE2_FILTER_RNASEQ }       from '../../modules/local/bowtie2_filter_rnaseq.nf'
include { FASTQC as RNASEQ_RAW_FASTQC }     from '../../modules/local/fastqc.nf'
include { FASTQC as RNASEQ_CLIPPED_FASTQC } from '../../modules/local/fastqc.nf'

workflow RNASEQ_PREPROCESS {
    take:
    ch_reads          // [ meta(id,lane,strand,single_end), reads ]
    ch_filter_index   // value: [ index_base, [index_files] ]

    main:
    def dedup = Utils.resolve_rnaseq_dedup_method(params)

    // Raw-read QC (SE file or PE [R1,R2]; gated on params.do_fastqc in the module).
    RNASEQ_RAW_FASTQC(ch_reads)

    CUTADAPT_CLIP_RNASEQ(ch_reads)

    // Rejoin PE R1+R2 clipped reads into a single [meta, reads] channel (SE stays a
    // single file). remainder:true keeps SE lanes whose reads2 is absent.
    ch_clipped = CUTADAPT_CLIP_RNASEQ.out.reads
        .join(CUTADAPT_CLIP_RNASEQ.out.reads2, remainder: true)
        .map { meta, r1, r2 -> r2 ? [meta, [r1, r2]] : [meta, r1] }

    // Post-clip QC.
    RNASEQ_CLIPPED_FASTQC(ch_clipped)

    // dedup=umicollapse extracts UMIs before filtering (SE or PE). Rejoin R1+R2 again.
    if (dedup == 'umicollapse') {
        UMITOOLS_EXTRACT(ch_clipped)
        ch_filter_in = UMITOOLS_EXTRACT.out.reads
            .join(UMITOOLS_EXTRACT.out.reads2, remainder: true)
            .map { meta, r1, r2 -> r2 ? [meta, [r1, r2]] : [meta, r1] }
    } else {
        ch_filter_in = ch_clipped
    }

    BOWTIE2_FILTER_RNASEQ(ch_filter_in, ch_filter_index)

    // Rejoin PE R1+R2 unaligned reads into a single [meta, reads] channel.
    ch_reads_joined = BOWTIE2_FILTER_RNASEQ.out.unaligned
        .join(BOWTIE2_FILTER_RNASEQ.out.unaligned2, remainder: true)
        .map { meta, r1, r2 -> r2 ? [meta, [r1, r2]] : [meta, r1] }

    emit:
    reads_for_genome = ch_reads_joined                  // [ meta, reads ] (SE file | PE [R1,R2])
    clip_log         = CUTADAPT_CLIP_RNASEQ.out.log     // [ meta, log ]
    filter_log       = BOWTIE2_FILTER_RNASEQ.out.log    // [ meta, log ]
}
