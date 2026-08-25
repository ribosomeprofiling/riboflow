// Bowtie2 transcriptome alignment for the RNA-seq path. SE/PE aware. Mirrors
// bowtie2_transcriptome.nf; outputs a coord-sorted BAM + log, named with the
// rnaseq.transcriptome_alignment.* prefix (storeDir separation).
process BOWTIE2_TRANSCRIPTOME_RNASEQ {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(reads)
    tuple val(index_base), path(index_files)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.transcriptome_alignment.bam"), emit: bam
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.transcriptome_alignment.log"), emit: log

    script:
    def prefix       = "${meta.id}.${meta.lane}"
    def args         = task.ext.args ?: (params.rnaseq?.transcriptome?.alignment_arguments ?: '-L 15 --no-unal')
    def is_pe        = meta.single_end == false
    def reads_arg    = is_pe ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads}"
    def aln_threads  = Math.min(task.cpus as int, 16)
    def sort_threads = Math.min(task.cpus as int, 8)
    def sort_mem     = Utils.samtools_sort_mem_per_thread_mb(task)
    """
    set -o pipefail
    bowtie2 ${args} \\
            -x ${index_base} ${reads_arg} \\
            --threads ${aln_threads} \\
            --rg-id "${prefix}" \\
            --rg "SM:${meta.id}" \\
            --rg "LB:${prefix}" \\
            --rg "PL:ILLUMINA" \\
                     2> ${prefix}.rnaseq.transcriptome_alignment.log \\
            | samtools sort -@ ${sort_threads} -m ${sort_mem}M \\
                     -o ${prefix}.rnaseq.transcriptome_alignment.bam -
    samtools index -@ ${sort_threads} ${prefix}.rnaseq.transcriptome_alignment.bam
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    """
    touch ${prefix}.rnaseq.transcriptome_alignment.bam
    touch ${prefix}.rnaseq.transcriptome_alignment.log
    """
}
