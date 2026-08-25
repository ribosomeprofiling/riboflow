// Extract reads from a qpass merged BAM matching an rfc-dedup BED →
// sample-level post-dedup BAM. With ext.emit_counts=true also emits merged
// total (+ primary/secondary/unique in multi-mapper mode) counts. Ports
// `genome_convert_dedup_bed_to_bam_position` (RiboFlow.groovy:1162-1198) and
// `transcriptome_convert_dedup_bed_to_bam_position` (:859-883, no counts).
process RFC_EXTRACT_DEDUP_READS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(dedup_bed), path(qpass_bam)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.merged_dedup.total.count"),
                     optional: true, emit: total_count
    tuple val(meta), path("${meta.id}.merged_dedup.primary.count"),
                     path("${meta.id}.merged_dedup.secondary.count"),
                     path("${meta.id}.merged_dedup.unique.count"),
                     optional: true, emit: detail_counts

    script:
    prefix               = task.ext.prefix ?: "${meta.id}.post_dedup"
    def route            = task.ext.route ?: 'genome'
    def emit_counts      = task.ext.emit_counts ?: false
    def emit_full_counts = (task.ext.emit_full_counts != null) ? task.ext.emit_full_counts : !Utils.route_unique_only(params, route)
    def is_pe            = meta.single_end == false
    def counts           = emit_counts
        ? Utils.samtools_count_block(Utils.route_count_args(params, route), is_pe,
                                     "${prefix}.bam", "${meta.id}.merged_dedup", task.cpus as int,
                                     emit_full_counts, emit_full_counts)
        : ''
    """
    if [ -s ${dedup_bed} ]; then
        rfc extract-dedup-reads \\
            --bam ${qpass_bam} \\
            --bed ${dedup_bed} \\
            --output ${prefix}.bam
    else
        samtools view -H ${qpass_bam} | samtools view -b -o ${prefix}.bam
    fi
    samtools index -@ ${task.cpus} ${prefix}.bam
    ${counts}
    """

    stub:
    prefix               = task.ext.prefix ?: "${meta.id}.post_dedup"
    def route            = task.ext.route ?: 'genome'
    def emit_counts      = task.ext.emit_counts ?: false
    def emit_full_counts = (task.ext.emit_full_counts != null) ? task.ext.emit_full_counts : !Utils.route_unique_only(params, route)
    def total_cmd  = emit_counts ? "echo 0 > ${meta.id}.merged_dedup.total.count" : ''
    def detail_cmd = (emit_counts && emit_full_counts) ? "echo 0 > ${meta.id}.merged_dedup.primary.count; echo 0 > ${meta.id}.merged_dedup.secondary.count; echo 0 > ${meta.id}.merged_dedup.unique.count" : ''
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    ${total_cmd}
    ${detail_cmd}
    """
}
