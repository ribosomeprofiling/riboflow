// Split a merged post-dedup BED back into per-lane BEDs using the sample.lane id
// in column 7, in one pass per sample. With ext.emit_counts it also writes
// per-lane total counts, plus primary, secondary and unique in multi-mapper mode.
//
// Counts on a BED: total is the record count; primary is the number of distinct
// read names; secondary is total minus primary; unique is the number of records
// whose MAPQ (column 5) reaches the threshold from
// <route>.samtools_count_arguments.unique, which defaults to -q 255.
//
// The output globs start with the sample id because every sample of a route
// shares one storeDir, and a bare *.count glob would let a sibling sample's
// stored files satisfy this task's outputs, so Nextflow would skip it.
process SEPARATE_BED {
    tag "${meta.id}"

    input:
    tuple val(meta), path(merged_bed), val(lanes)

    output:
    tuple val(meta), path("${meta.id}.*.${task.ext.suffix ?: 'genome'}.post_dedup.bed"), emit: beds
    tuple val(meta), path("${meta.id}.*.dedup.total.count"),                        optional: true, emit: total_counts
    tuple val(meta), path("${meta.id}.*.dedup.primary.count"), path("${meta.id}.*.dedup.secondary.count"),
                     path("${meta.id}.*.dedup.unique.count"),                       optional: true, emit: detail_counts

    script:
    def suffix           = task.ext.suffix ?: 'genome'
    def route            = task.ext.route ?: 'genome'
    def emit_counts      = task.ext.emit_counts ?: false
    def emit_full_counts = (task.ext.emit_full_counts != null) ? task.ext.emit_full_counts : !Utils.route_unique_only(params, route)
    def uniq_mapq        = Utils.count_mapq_threshold(Utils.route_count_args(params, route).unique?.toString())
    def detail_cmd = (emit_counts && emit_full_counts) ? """
    pids=()
    for l in ${lanes.join(' ')}; do
        (
        p=\$(cut -f4 "\$l".${suffix}.post_dedup.bed | LC_ALL=C sort -u | wc -l | tr -d ' ')
        t=\$(cat "\$l".dedup.total.count)
        echo "\$p" > "\$l".dedup.primary.count
        echo \$((t - p)) > "\$l".dedup.secondary.count
        ) & pids+=(\$!)
    done
    for p in \${pids[@]+"\${pids[@]}"}; do wait "\$p"; done
    """ : ''
    def cleanup = !emit_counts ? 'rm -f *.dedup.total.count *.dedup.unique.count' \
                : (!emit_full_counts ? 'rm -f *.dedup.unique.count' : '')
    """
    awk -F'\\t' -v OFS='\\t' -v sfx='${suffix}' -v U=${uniq_mapq} '
        { f = \$7 "." sfx ".post_dedup.bed"
          print \$1, \$2, \$3, \$4, \$5, \$6 > f
          tot[\$7]++
          if (\$5 + 0 >= U) uq[\$7]++ }
        END { for (l in tot) {
                  print tot[l]     > (l ".dedup.total.count")
                  print uq[l] + 0  > (l ".dedup.unique.count") } }' ${merged_bed}
    for l in ${lanes.join(' ')}; do
        [ -f "\$l".${suffix}.post_dedup.bed ] || : > "\$l".${suffix}.post_dedup.bed
        [ -f "\$l".dedup.total.count ]        || echo 0 > "\$l".dedup.total.count
        [ -f "\$l".dedup.unique.count ]       || echo 0 > "\$l".dedup.unique.count
    done
    ${detail_cmd}
    ${cleanup}
    """

    stub:
    def suffix           = task.ext.suffix ?: 'genome'
    def route            = task.ext.route ?: 'genome'
    def emit_counts      = task.ext.emit_counts ?: false
    def emit_full_counts = (task.ext.emit_full_counts != null) ? task.ext.emit_full_counts : !Utils.route_unique_only(params, route)
    def total_cmd  = emit_counts ? "echo 0 > \"\$l\".dedup.total.count" : ''
    def detail_cmd = (emit_counts && emit_full_counts) ? "echo 0 > \"\$l\".dedup.primary.count; echo 0 > \"\$l\".dedup.secondary.count; echo 0 > \"\$l\".dedup.unique.count" : ''
    """
    for l in ${lanes.join(' ')}; do
        touch "\$l".${suffix}.post_dedup.bed
        ${total_cmd}
        ${detail_cmd}
    done
    """
}
