// BAM to BED.
//
// With ext.append_index, each record gets the sample.lane id in column 7.
// With ext.indexed_prefix, that tagged copy is written in the same pass. An empty
// BAM yields empty BEDs.
//
// prefix and indexed_prefix are resolved from task.ext in output:, script: and
// stub: separately, and are always def-scoped.
process BAM_TO_BED {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id + '.' + meta.lane + '.genome.qpass'}.bed"), emit: bed
    tuple val(meta), path("${task.ext.indexed_prefix ?: (task.ext.prefix ?: meta.id + '.' + meta.lane + '.genome.qpass') + '.with_sample_index'}.bed"), optional: true, emit: indexed

    script:
    def prefix         = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome.qpass"
    def indexed_prefix = task.ext.indexed_prefix ?: "${prefix}.with_sample_index"
    def s          = "${meta.id}.${meta.lane}"
    def tag        = "awk -v s=${s} '{ print \$0\"\\t\"s }'"
    def cmd
    if (task.ext.indexed_prefix) {
        cmd = "bamToBed -i ${bam} | tee ${prefix}.bed | ${tag} > ${indexed_prefix}.bed"
    } else if (task.ext.append_index ?: false) {
        cmd = "bamToBed -i ${bam} | ${tag} > ${prefix}.bed"
    } else {
        cmd = "bamToBed -i ${bam} > ${prefix}.bed"
    }
    """
    ${cmd}
    """

    stub:
    def prefix         = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome.qpass"
    def indexed_prefix = task.ext.indexed_prefix ?: "${prefix}.with_sample_index"
    def extra      = task.ext.indexed_prefix ? "touch ${indexed_prefix}.bed" : ''
    """
    touch ${prefix}.bed
    ${extra}
    """
}
