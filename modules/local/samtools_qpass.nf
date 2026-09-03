// Quality-pass filter: a MAPQ cutoff plus a user-supplied record filter, with one
// count file per read type. With ext.presort the BAM is coordinate-sorted after
// filtering, which the QNAME-sorted STAR transcriptome BAM needs.
//
// ext.route picks the params block: genome, transcriptome, rnaseq_genome,
// rnaseq_transcriptome or star_transcriptome. Utils resolves from it:
//   -q            255 when <route>.unique_only, else <route>.mapping_quality_cutoff
//   record filter <route>.samtools_filter_arguments, re-quoted token by token
//   counter flags <route>.samtools_count_arguments
//
// ext.emit_primary_secondary and ext.count_unique default to multi-mapper mode
// only on the genome routes; the transcriptome routes set both false.
//
// prefix is resolved from task.ext in output:, script: and stub: separately, and
// is always def-scoped.
process SAMTOOLS_QPASS {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(bam)

    // `prefix` is resolved independently here, in `script:` and in `stub:`, and is
    // always `def`-scoped. An undeclared `prefix = ...` lands in the session-global
    // script binding shared by every concurrent task, leaking one process's prefix
    // into another's output filenames.
    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id + '.' + meta.lane + '.genome_alignment'}.qpass.bam"),
                     path("${task.ext.prefix ?: meta.id + '.' + meta.lane + '.genome_alignment'}.qpass.bam.bai"), emit: bam
    tuple val(meta), path("${task.ext.prefix ?: meta.id + '.' + meta.lane + '.genome_alignment'}.qpass.total.count"),     emit: total_count
    tuple val(meta), path("${task.ext.prefix ?: meta.id + '.' + meta.lane + '.genome_alignment'}.qpass.primary.count"),   optional: true, emit: primary_count
    tuple val(meta), path("${task.ext.prefix ?: meta.id + '.' + meta.lane + '.genome_alignment'}.qpass.secondary.count"), optional: true, emit: secondary_count
    tuple val(meta), path("${task.ext.prefix ?: meta.id + '.' + meta.lane + '.genome_alignment'}.qpass.unique.count"),    optional: true, emit: unique_count

    script:
    def prefix       = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome_alignment"
    def route        = task.ext.route ?: 'genome'
    def presort      = task.ext.presort ?: false
    def unique_only  = Utils.route_unique_only(params, route)
    def mapq         = Utils.route_effective_mapq(params, route)
    def filter_args  = Utils.shell_quote_args(Utils.route_filter_args(params, route))
    def emit_ps      = (task.ext.emit_primary_secondary != null) ? task.ext.emit_primary_secondary : !unique_only
    def count_unique = (task.ext.count_unique != null) ? task.ext.count_unique : !unique_only
    def is_pe        = meta.single_end == false
    def sort_threads = Math.min(task.cpus as int, 8)
    def sort_mem     = Utils.samtools_sort_mem_per_thread_mb(task)
    // `-b -q`, not `-bq`, so the user's arguments append cleanly after the cutoff.
    def make_bam     = presort \
        ? "samtools view -h -b -q ${mapq} ${filter_args} ${bam} | samtools sort -@ ${sort_threads} -m ${sort_mem}M -o ${prefix}.qpass.bam -" \
        : "samtools view -@ ${task.cpus} -b -q ${mapq} ${filter_args} ${bam} > ${prefix}.qpass.bam"
    def counts       = Utils.samtools_count_block(Utils.route_count_args(params, route), is_pe,
                                                  "${prefix}.qpass.bam", "${prefix}.qpass", task.cpus as int,
                                                  emit_ps, count_unique)
    """
    ${make_bam}
    samtools index -@ ${task.cpus} ${prefix}.qpass.bam
    ${counts}
    """

    stub:
    def prefix       = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome_alignment"
    def route        = task.ext.route ?: 'genome'
    def unique_only  = Utils.route_unique_only(params, route)
    def emit_ps      = (task.ext.emit_primary_secondary != null) ? task.ext.emit_primary_secondary : !unique_only
    def count_unique = (task.ext.count_unique != null) ? task.ext.count_unique : !unique_only
    def ps_cmd   = emit_ps      ? "echo 0 > ${prefix}.qpass.primary.count; echo 0 > ${prefix}.qpass.secondary.count" : ''
    def uniq_cmd = count_unique ? "echo 0 > ${prefix}.qpass.unique.count" : ''
    """
    touch ${prefix}.qpass.bam ${prefix}.qpass.bam.bai
    echo 0 > ${prefix}.qpass.total.count
    ${ps_cmd}
    ${uniq_cmd}
    """
}
