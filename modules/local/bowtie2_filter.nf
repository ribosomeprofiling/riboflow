// rRNA/tRNA contaminant filter. Ports `filter` (RiboFlow.groovy:344-385).
// Unaligned reads feed genome alignment. The sorted/indexed filter BAM and the
// aligned FASTQ are kept deliberately so users can inspect what was removed.
// bowtie2 reads gzip input natively — no decompression process substitution.
process BOWTIE2_FILTER {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(fastq)
    tuple val(index_base), path(index_files)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.filter.bam"),                  emit: bam
    tuple val(meta), path("${meta.id}.${meta.lane}.aligned.filter.fastq.gz"),     emit: aligned
    tuple val(meta), path("${meta.id}.${meta.lane}.unaligned.filter.fastq.gz"),   emit: unaligned
    tuple val(meta), path("${meta.id}.${meta.lane}.filter.log"),                  emit: log

    script:
    def prefix       = "${meta.id}.${meta.lane}"
    def args         = task.ext.args ?: (params.alignment_arguments?.filter ?: '')
    def aln_threads  = Math.min(task.cpus as int, 16)
    def sort_threads = Math.min(task.cpus as int, 8)
    def sort_mem     = Utils.samtools_sort_mem_per_thread_mb(task)
    """
    set -o pipefail
    bowtie2 ${args} \\
            -x ${index_base} -U ${fastq} \\
            --threads ${aln_threads} \\
            --al-gz ${prefix}.aligned.filter.fastq.gz \\
            --un-gz ${prefix}.unaligned.filter.fastq.gz \\
                     2> ${prefix}.filter.log \\
            | samtools sort -@ ${sort_threads} -m ${sort_mem}M -o ${prefix}.filter.bam - \\
            && samtools index -@ ${sort_threads} ${prefix}.filter.bam
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    """
    touch ${prefix}.filter.bam
    echo | gzip -c > ${prefix}.aligned.filter.fastq.gz
    echo | gzip -c > ${prefix}.unaligned.filter.fastq.gz
    touch ${prefix}.filter.log
    """
}
