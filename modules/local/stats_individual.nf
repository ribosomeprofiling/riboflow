// Per-lane alignment-stats row for the genome routes — `rfc stats-individual
// --route genome` parses the cutadapt / bowtie2 filter / STAR logs by label,
// reads the qpass/dedup count files and checks the read accounting
// (aligned_once + aligned_many + unaligned == filter_kept) before writing the CSV.
// Ports `individual_genome_alignment_stats` (RiboFlow.groovy:1528-1617).
//
// Shared by the ribo-seq and RNA-seq genome paths:
//   ext.route       genome | rnaseq_genome — selects the stats mode via Utils
//                   (unique-only vs multi-mapper, i.e. <route>.unique_only).
//   ext.stats_label (default 'genome') sets the output CSV label.
process STATS_INDIVIDUAL {
    tag "${meta.id}.${meta.lane}"

    input:
    // qpass_* are staged under fixed distinct names: the unique-only proxy feeds the
    // same total-count file into the primary slot, so without renaming the two would
    // collide on a shared basename. dedup_* likewise never collide with qpass_* when
    // dedup_method=='none' aliases qpass counts into the dedup slot.
    tuple val(meta),
          path(clip_log), path(filter_log), path(genome_log), path(genome_secondary_count),
          path(qpass_total,     stageAs: 'qpass_total.count'),
          path(qpass_primary,   stageAs: 'qpass_primary.count'),
          path(qpass_secondary, stageAs: 'qpass_secondary.count'),
          path('dedup.total.count'), path('dedup.primary.count'), path('dedup.secondary.count'),
          path('dedup.unique.count'),
          path(qpass_unique, stageAs: 'qpass_unique.count')

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.${label}_individual.csv"), emit: csv

    script:
    def prefix      = "${meta.id}.${meta.lane}"
    label           = task.ext.stats_label ?: 'genome'
    def route       = task.ext.route ?: 'genome'
    def unique_only = Utils.route_unique_only(params, route)
    def mode_args   = unique_only ? '--unique-only' : '--multi --qpass-unique qpass_unique.count --dedup-unique dedup.unique.count'
    """
    rfc stats-individual --route genome --prefix ${prefix} \\
        --clip-log ${clip_log} --filter-log ${filter_log} --align-log ${genome_log} \\
        --genome-secondary-count ${genome_secondary_count} \\
        --qpass-total qpass_total.count --qpass-primary qpass_primary.count --qpass-secondary qpass_secondary.count \\
        --dedup-total dedup.total.count --dedup-primary dedup.primary.count --dedup-secondary dedup.secondary.count \\
        ${mode_args} \\
        -o ${prefix}.${label}_individual.csv
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    label = task.ext.stats_label ?: 'genome'
    """
    printf ',${prefix}\\ntotal_reads,0\\n' > ${prefix}.${label}_individual.csv
    """
}
