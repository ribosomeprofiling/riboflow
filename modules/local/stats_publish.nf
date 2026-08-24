// Publish the final stats CSVs. Ports `publish_stats` (RiboFlow.groovy:1758-1775).
process STATS_PUBLISH {

    input:
    path(individual_stats)
    path(merged_stats)

    output:
    path("${task.ext.prefix ?: ''}individual_stats.csv"), emit: individual
    path("${task.ext.prefix ?: ''}stats.csv"),            emit: merged

    script:
    def pfx = task.ext.prefix ?: ''
    """
    cp ${individual_stats} ${pfx}individual_stats.csv
    cp ${merged_stats} ${pfx}stats.csv
    """

    stub:
    def pfx = task.ext.prefix ?: ''
    """
    touch ${pfx}individual_stats.csv ${pfx}stats.csv
    """
}
