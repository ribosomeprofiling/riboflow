// Merge per-lane BAMs into one per-sample BAM. Ports `merge_genome_qpass_bam`
// (RiboFlow.groovy:634-651) and the transcriptome merge variants (:830,896).
// ext.prefix sets the output basename so the genome / transcriptome / dedup
// flavours keep their original filenames.
process SAMTOOLS_MERGE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam

    script:
    prefix = task.ext.prefix ?: "${meta.id}.genome.qpass.merged"
    """
    samtools merge -@ ${task.cpus} ${prefix}.bam ${bams}
    samtools index -@ ${task.cpus} ${prefix}.bam
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.genome.qpass.merged"
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    """
}
