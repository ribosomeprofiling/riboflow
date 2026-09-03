// Split a merged dedup BAM back into per-lane BAMs by read group, using one
// `samtools split` pass per sample. With ext.emit_bed_counts it also writes each
// lane's BED and total count, plus primary, secondary and unique in multi-mapper
// mode.
//
// `lanes` lists the sample.lane read-group ids expected. A lane with no reads
// left after dedup still gets a header-only BAM, so downstream joins never drop
// it.
//
// The output globs start with the sample id because every sample of a route
// shares one storeDir, and a bare *.count glob would let a sibling sample's
// stored files satisfy this task's outputs, so Nextflow would skip it.
process SPLIT_DEDUP_BAM {
    tag "${meta.id}"

    input:
    tuple val(meta), path(merged_bam), path(merged_bai), val(lanes)

    output:
    tuple val(meta), path("${meta.id}.*.${task.ext.suffix ?: 'post_dedup'}.bam"),
                     path("${meta.id}.*.${task.ext.suffix ?: 'post_dedup'}.bam.bai"), emit: bams
    tuple val(meta), path("${meta.id}.*.${task.ext.bed_suffix ?: 'genome.post_dedup'}.bed"), optional: true, emit: beds
    tuple val(meta), path("${meta.id}.*.dedup.total.count"),                        optional: true, emit: total_counts
    tuple val(meta), path("${meta.id}.*.dedup.primary.count"), path("${meta.id}.*.dedup.secondary.count"),
                     path("${meta.id}.*.dedup.unique.count"),                       optional: true, emit: detail_counts

    script:
    def suffix           = task.ext.suffix     ?: 'post_dedup'
    def bed_suffix       = task.ext.bed_suffix ?: 'genome.post_dedup'
    def route            = task.ext.route ?: 'genome'
    def emit_extra       = task.ext.emit_bed_counts ?: false
    def emit_full_counts = (task.ext.emit_full_counts != null) ? task.ext.emit_full_counts : !Utils.route_unique_only(params, route)
    def is_pe            = meta.single_end == false
    def per_lane_threads = Math.max(1, (task.cpus as int).intdiv(Math.max(1, lanes.size())))
    def counts           = emit_extra
        ? Utils.samtools_count_block(Utils.route_count_args(params, route), is_pe,
                                     "\"\$l\".${suffix}.bam", "\"\$l\".dedup", per_lane_threads,
                                     emit_full_counts, emit_full_counts)
        : ''
    def bed_cmd          = emit_extra ? "bamToBed -i \"\$l\".${suffix}.bam > \"\$l\".${bed_suffix}.bed" : ''
    """
    if ! samtools view -H ${merged_bam} | grep -q "^@RG"; then
        echo "ERROR: No read groups found in ${merged_bam}" >&2; exit 1
    fi
    samtools split -@ ${task.cpus} -f '%!.${suffix}.bam' -u unassigned_reads.bam ${merged_bam}
    for l in ${lanes.join(' ')}; do
        if [ ! -f "\$l".${suffix}.bam ]; then
            samtools view -H ${merged_bam} | samtools view -b -o "\$l".${suffix}.bam
        fi
    done
    lane_pids=()
    for l in ${lanes.join(' ')}; do
        (
        samtools index -@ ${per_lane_threads} "\$l".${suffix}.bam
        ${counts}
        ${bed_cmd}
        ) & lane_pids+=(\$!)
    done
    for p in \${lane_pids[@]+"\${lane_pids[@]}"}; do wait "\$p"; done
    rm -f unassigned_reads.bam
    """

    stub:
    def suffix           = task.ext.suffix     ?: 'post_dedup'
    def bed_suffix       = task.ext.bed_suffix ?: 'genome.post_dedup'
    def route            = task.ext.route ?: 'genome'
    def emit_extra       = task.ext.emit_bed_counts ?: false
    def emit_full_counts = (task.ext.emit_full_counts != null) ? task.ext.emit_full_counts : !Utils.route_unique_only(params, route)
    def extra = emit_extra ? "touch \"\$l\".${bed_suffix}.bed; echo 0 > \"\$l\".dedup.total.count" : ''
    def full  = (emit_extra && emit_full_counts) ? "echo 0 > \"\$l\".dedup.primary.count; echo 0 > \"\$l\".dedup.secondary.count; echo 0 > \"\$l\".dedup.unique.count" : ''
    """
    for l in ${lanes.join(' ')}; do
        touch "\$l".${suffix}.bam "\$l".${suffix}.bam.bai
        ${extra}
        ${full}
    done
    """
}
