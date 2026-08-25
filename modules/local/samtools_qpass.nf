// Quality-pass filter (MAPQ + user-supplied record filter) with per-type count files.
// Ports `genome_quality_filter` (RiboFlow.groovy:576-609) and, with
// ext.presort=true, `transcriptome_sort_and_filter` (:675-701) — the STAR
// transcriptome BAM is QNAME-sorted so it must be coord-sorted after filtering.
//
// ext.route  selects the params block (genome | transcriptome | rnaseq_genome |
//            rnaseq_transcriptome | star_transcriptome). From it Utils resolves:
//              -q           <route>.unique_only ⇒ 255, else <route>.mapping_quality_cutoff
//              record filter <route>.samtools_filter_arguments (re-quoted token by token)
//              counter flags <route>.samtools_count_arguments
// ext.emit_primary_secondary / ext.count_unique default to "only in multi-mapper
//            mode" on the genome routes; the transcriptome routes set both false.
process SAMTOOLS_QPASS {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.qpass.bam"), path("${prefix}.qpass.bam.bai"), emit: bam
    tuple val(meta), path("${prefix}.qpass.total.count"),   emit: total_count
    tuple val(meta), path("${prefix}.qpass.primary.count"), optional: true, emit: primary_count
    tuple val(meta), path("${prefix}.qpass.secondary.count"), optional: true, emit: secondary_count
    tuple val(meta), path("${prefix}.qpass.unique.count"),  optional: true, emit: unique_count

    script:
    prefix           = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome_alignment"
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
    prefix           = task.ext.prefix ?: "${meta.id}.${meta.lane}.genome_alignment"
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
