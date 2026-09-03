// Merge the per-lane BAMs into one BAM per sample. ext.prefix sets the output
// basename, so the genome, transcriptome and dedup variants keep distinct names.
//
// prefix is resolved from task.ext in output:, script: and stub: separately, and
// is always def-scoped.
process SAMTOOLS_MERGE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta),
          path("${task.ext.prefix ?: meta.id + '.genome.qpass.merged'}.bam"),
          path("${task.ext.prefix ?: meta.id + '.genome.qpass.merged'}.bam.bai"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.genome.qpass.merged"
    """
    samtools merge -@ ${task.cpus} ${prefix}.bam ${bams}
    samtools index -@ ${task.cpus} ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.genome.qpass.merged"
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    """
}
