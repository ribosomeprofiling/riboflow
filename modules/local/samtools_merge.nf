// Merge the per-lane BAMs into one BAM per sample. ext.prefix sets the output
// basename, so the genome, transcriptome and dedup variants keep distinct names.
//
// With ext.unsorted the lanes are concatenated (samtools cat) under a header that
// carries every lane's @RG, and no index is written: the STAR transcriptome-
// projection BAM stays grouped by read name, as Salmon requires.
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
          path("${task.ext.prefix ?: meta.id + '.genome.qpass.merged'}.bam.bai", arity: '0..*'), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.genome.qpass.merged"
    if (task.ext.unsorted)
        """
        set -o pipefail
        first=\$(echo ${bams} | cut -d' ' -f1)
        { samtools view -H "\$first" | grep -v '^@RG'
          for b in ${bams}; do samtools view -H "\$b" | { grep '^@RG' || true; }; done; } > header.sam
        samtools cat -h header.sam -o ${prefix}.bam ${bams}
        """
    else
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
