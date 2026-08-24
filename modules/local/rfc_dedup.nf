// Coordinate-based dedup via `rfc dedup`. Ports `genome_deduplicate_position`
// (RiboFlow.groovy:1069-1088) and `transcriptome_deduplicate_position` (:764-784).
// Every producer (CONCAT_SORT_BED and its aliases) already emits a coordinate-
// sorted BED, so the input is passed straight through — this used to re-sort the
// whole-sample BED a second time.
process RFC_DEDUP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${prefix}.post_dedup.bed"), emit: bed

    script:
    prefix = task.ext.prefix ?: "${meta.id}.genome"
    """
    rfc dedup -i ${bed} -o ${prefix}.post_dedup.bed
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.genome"
    """
    touch ${prefix}.post_dedup.bed
    """
}
