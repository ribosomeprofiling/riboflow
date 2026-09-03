// Bowtie2 alignment against the transcriptome index. Outputs a coordinate-sorted
// BAM and the bowtie2 log, plus the aligned and unaligned FASTQs for inspection
// and optional FastQC. The BAM feeds the transcriptome dedup and ribopy steps.
process BOWTIE2_TRANSCRIPTOME {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(fastq)
    tuple val(index_base), path(index_files)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.transcriptome_alignment.bam"), emit: bam
    tuple val(meta), path("${meta.id}.${meta.lane}.transcriptome_alignment.aligned.fastq.gz"),   emit: aligned
    tuple val(meta), path("${meta.id}.${meta.lane}.transcriptome_alignment.unaligned.fastq.gz"), emit: unaligned
    tuple val(meta), path("${meta.id}.${meta.lane}.transcriptome_alignment.log"), emit: log

    script:
    def prefix       = "${meta.id}.${meta.lane}"
    def args         = task.ext.args ?: (params.alignment_arguments?.transcriptome ?: '')
    def aln_threads  = Math.min(task.cpus as int, 16)
    def sort_threads = Math.min(task.cpus as int, 8)
    def sort_mem     = Utils.samtools_sort_mem_per_thread_mb(task)
    """
    set -o pipefail
    bowtie2 ${args} \\
            -x ${index_base} -U ${fastq} \\
            --threads ${aln_threads} \\
            --un-gz ${prefix}.transcriptome_alignment.unaligned.fastq.gz \\
            --rg-id "${prefix}" \\
            --rg "SM:${meta.id}" \\
            --rg "LB:${prefix}" \\
            --rg "PL:ILLUMINA" \\
                     2> ${prefix}.transcriptome_alignment.log \\
            | samtools sort -@ ${sort_threads} -m ${sort_mem}M \\
                     -o ${prefix}.transcriptome_alignment.bam -
    samtools index -@ ${sort_threads} ${prefix}.transcriptome_alignment.bam

    # Aligned reads → FASTQ, primary records only. --un-gz above emits the unaligned.
    samtools fastq -@ ${sort_threads} -F 2308 ${prefix}.transcriptome_alignment.bam \\
        | pigz -p ${task.cpus} > ${prefix}.transcriptome_alignment.aligned.fastq.gz
    # bowtie2 may skip --un-gz output when there are no unaligned reads; ensure it exists.
    [ -f ${prefix}.transcriptome_alignment.unaligned.fastq.gz ] || echo -n | gzip > ${prefix}.transcriptome_alignment.unaligned.fastq.gz
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    """
    touch ${prefix}.transcriptome_alignment.bam
    echo | gzip -c > ${prefix}.transcriptome_alignment.aligned.fastq.gz
    echo | gzip -c > ${prefix}.transcriptome_alignment.unaligned.fastq.gz
    touch ${prefix}.transcriptome_alignment.log
    """
}
