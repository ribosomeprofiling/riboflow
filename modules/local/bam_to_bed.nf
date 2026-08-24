// BAM → BED. Ports `individual_genome_bam_to_bed` (RiboFlow.groovy:983-1002) and,
// with ext.append_index=true, `transcriptome_qpass_bam_to_bed` (:715-734) which
// tags each record with the sample.lane id in column 7.
//
// ext.indexed_prefix: additionally write a 7-column (sample.lane-tagged) copy in
// the same pass (`tee`), replacing the former ADD_SAMPLE_INDEX_COL re-read. An
// empty BAM simply yields empty BEDs — no separate `samtools view -c` pass.
process BAM_TO_BED {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bed"),         emit: bed
    tuple val(meta), path("${indexed_prefix}.bed"), optional: true, emit: indexed

    script:
    prefix         = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome.qpass"
    indexed_prefix = task.ext.indexed_prefix ?: "${prefix}.with_sample_index"
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
    prefix         = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome.qpass"
    indexed_prefix = task.ext.indexed_prefix ?: "${prefix}.with_sample_index"
    def extra      = task.ext.indexed_prefix ? "touch ${indexed_prefix}.bed" : ''
    """
    touch ${prefix}.bed
    ${extra}
    """
}
