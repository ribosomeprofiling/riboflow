// Combine the per-row stats CSVs and add percentage rows. ext.prefix sets the
// output basename, for example genome_individual_essential.
process STATS_COMBINE {

    input:
    path(stat_tables)

    output:
    path("${task.ext.prefix ?: 'genome_individual_essential'}.csv"), emit: csv

    script:
    def prefix      = task.ext.prefix ?: 'genome_individual_essential'
    def stats_label    = task.ext.stats_label ?: 'genome'
    def unique_only    = (task.ext.unique_only != null) ? task.ext.unique_only : Utils.route_unique_only(params, task.ext.route ?: 'genome')
    def unique_flag    = unique_only ? '--unique-only' : ''
    def pct_cmd        = (stats_label == 'transcriptome')
        ? "rfc stats-percentage --label-prefix transcriptome -i raw_${prefix}.csv -o ${prefix}.csv"
        : "rfc genome-stats-percentage ${unique_flag} -i raw_${prefix}.csv -o ${prefix}.csv"
    def n           = stat_tables instanceof List ? stat_tables.size() : 1
    if (n == 0) {
        """
        echo "No statistics data available" > ${prefix}.csv
        """
    } else {
        """
        rfc merge overall-stats -o raw_${prefix}.csv ${stat_tables}
        ${pct_cmd}
        """
    }

    stub:
    def prefix = task.ext.prefix ?: 'genome_individual_essential'
    """
    touch ${prefix}.csv
    """
}
