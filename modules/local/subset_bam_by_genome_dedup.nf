// Dedup the STAR transcriptome-projected BAM by read name: keep exactly the
// reads that survived UMI dedup on the genome BAM. One dedup decision per read,
// shared by every isoform it projects to, and no second UMICollapse run — the
// transcriptome BAM (read × isoform) is several times the genome BAM and
// UMICollapse --two-pass crashes (NPE) or OOMs on it for large libraries.
process SUBSET_BAM_BY_GENOME_DEDUP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(tx_bam), path(genome_dedup_bam)

    // No index: the projection BAM keeps STAR's read-name grouping (Salmon input).
    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id + '.dedup'}.bam"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.dedup"
    """
    samtools view -@ ${task.cpus} ${genome_dedup_bam} | cut -f1 | sort -u > keep.names
    samtools view -@ ${task.cpus} -b -N keep.names -o ${prefix}.bam ${tx_bam}
    rm keep.names
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.dedup"
    """
    touch ${prefix}.bam
    """
}
