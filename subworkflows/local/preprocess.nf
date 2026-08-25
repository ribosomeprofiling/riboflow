// Preprocess: clip → [UMI extract] → rRNA/tRNA filter. Emits the genome-ready
// unaligned reads plus the clip/filter logs the stats stage consumes.
// (RiboFlow.groovy:257-387, dedup routing :302-307)

include { CUTADAPT_CLIP }              from '../../modules/local/cutadapt_clip.nf'
include { UMITOOLS_EXTRACT }           from '../../modules/local/umitools_extract.nf'
include { BOWTIE2_FILTER }             from '../../modules/local/bowtie2_filter.nf'
include { FASTQC as RAW_FASTQC }       from '../../modules/local/fastqc.nf'
include { FASTQC as CLIPPED_FASTQC }   from '../../modules/local/fastqc.nf'

workflow PREPROCESS {
    take:
    ch_reads          // [ meta(id,lane,strand), fastq ]
    ch_filter_index   // value: [ index_base, [index_files] ]

    main:
    def dedup = Utils.resolve_dedup_method(params)

    RAW_FASTQC(ch_reads)

    CUTADAPT_CLIP(ch_reads)
    CLIPPED_FASTQC(CUTADAPT_CLIP.out.reads)

    // dedup=umicollapse extracts UMIs before filtering; otherwise filter the
    // clipped reads directly.
    ch_filter_in = (dedup == 'umicollapse') \
        ? UMITOOLS_EXTRACT(CUTADAPT_CLIP.out.reads).reads \
        : CUTADAPT_CLIP.out.reads

    BOWTIE2_FILTER(ch_filter_in, ch_filter_index)

    emit:
    reads_for_genome = BOWTIE2_FILTER.out.unaligned   // [ meta, fastq ]
    clip_log         = CUTADAPT_CLIP.out.log           // [ meta, log ]
    filter_log       = BOWTIE2_FILTER.out.log          // [ meta, log ]
}
