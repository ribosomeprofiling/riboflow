// rRNA/tRNA contaminant filter for the RNA-seq path, single- or paired-end.
// Unaligned reads feed STAR genome and bowtie2 transcriptome alignment.
// Paired-end uses --un-conc-gz, whose % expands to _R1 and _R2; single-end uses
// --un-gz and writes _R1 only.
//
// A pair is removed only when it aligns concordantly to the contaminant index.
// --no-mixed and --no-discordant stop bowtie2 reporting mixed or discordant hits,
// so the log's "aligned concordantly 0 times" equals the number of pairs written
// to --un-conc. Pairs where only one mate hits rRNA are kept.
process BOWTIE2_FILTER_RNASEQ {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(reads)
    tuple val(index_base), path(index_files)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.filter.bam"),                          emit: bam
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.unaligned_R1.fastq.gz"),               emit: unaligned
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.unaligned_R2.fastq.gz"), optional: true, emit: unaligned2
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.filter.log"),                          emit: log

    script:
    def prefix       = "${meta.id}.${meta.lane}"
    def args         = task.ext.args ?: (params.rnaseq?.filter_arguments ?: '-L 15 --no-unal')
    def is_pe        = meta.single_end == false
    def aln_threads  = Math.min(task.cpus as int, 16)
    def sort_threads = Math.min(task.cpus as int, 8)
    def sort_mem     = Utils.samtools_sort_mem_per_thread_mb(task)
    def reads_arg    = is_pe ? "-1 ${reads[0]} -2 ${reads[1]} --no-mixed --no-discordant" : "-U ${reads}"
    def unal_arg     = is_pe \
        ? "--un-conc-gz ${prefix}.rnaseq.unaligned_R%.fastq.gz" \
        : "--un-gz ${prefix}.rnaseq.unaligned_R1.fastq.gz"
    """
    set -o pipefail
    bowtie2 ${args} \\
            -x ${index_base} ${reads_arg} \\
            --threads ${aln_threads} \\
            ${unal_arg} \\
                     2> ${prefix}.rnaseq.filter.log \\
            | samtools sort -@ ${sort_threads} -m ${sort_mem}M -o ${prefix}.rnaseq.filter.bam - \\
            && samtools index -@ ${sort_threads} ${prefix}.rnaseq.filter.bam
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    def is_pe  = meta.single_end == false
    def r2_cmd = is_pe ? "echo | gzip -c > ${prefix}.rnaseq.unaligned_R2.fastq.gz" : ''
    """
    touch ${prefix}.rnaseq.filter.bam
    echo | gzip -c > ${prefix}.rnaseq.unaligned_R1.fastq.gz
    ${r2_cmd}
    touch ${prefix}.rnaseq.filter.log
    """
}
