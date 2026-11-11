// Number of distinct reads (fragments for paired-end: R1 records only) in a BAM.
// Used on the STAR transcriptome-projection BAMs, where every read appears once
// per compatible isoform, so `samtools view -c` would count records, not reads.
// Streams through an awk hash instead of sort -u: the projection BAMs reach
// >1e9 records and an external sort spills to disk for hours.
process COUNT_DISTINCT_READS {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.${task.ext.suffix ?: 'distinct'}.reads.count"), emit: count

    script:
    def suffix  = task.ext.suffix ?: 'distinct'
    def pe_flag = (meta.single_end == false) ? '-f 64' : ''
    """
    set -o pipefail
    samtools view -@ ${task.cpus} ${pe_flag} ${bam} \\
        | awk -F'\\t' '!s[\$1]++{n++} END{print n+0}' \\
        > ${meta.id}.${meta.lane}.${suffix}.reads.count
    """

    stub:
    def suffix = task.ext.suffix ?: 'distinct'
    """
    echo 0 > ${meta.id}.${meta.lane}.${suffix}.reads.count
    """
}
