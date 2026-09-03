// Per-lane transcriptome alignment stats row — `rfc stats-individual --route
// transcriptome` parses the bowtie2 transcriptome log + clip/filter logs and the
// qpass/dedup totals into a raw-count CSV (with the same read-accounting checks as
// the genome row). Includes clip/filter rows so transcriptome-only runs produce a
// self-contained stats file.
//
// Shared by the ribo-seq and RNA-seq transcriptome paths: `ext.stats_label`
// (default 'transcriptome') sets the output CSV label.
process TX_STATS_INDIVIDUAL {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta),
          path(tx_log), path(clip_log), path(filter_log),
          path(qpass_total),
          path('dedup.total.count')

    output:
    tuple val(meta),
          path("${meta.id}.${meta.lane}.${task.ext.stats_label ?: 'transcriptome'}_individual.csv"), emit: csv

    script:
    def prefix = "${meta.id}.${meta.lane}"
    def label = task.ext.stats_label ?: 'transcriptome'
    """
    rfc stats-individual --route transcriptome --prefix ${prefix} \\
        --clip-log ${clip_log} --filter-log ${filter_log} --align-log ${tx_log} \\
        --qpass-total ${qpass_total} --dedup-total dedup.total.count \\
        -o ${prefix}.${label}_individual.csv
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    def label = task.ext.stats_label ?: 'transcriptome'
    """
    printf ',${prefix}\\ntotal_reads,0\\n' > ${prefix}.${label}_individual.csv
    """
}
