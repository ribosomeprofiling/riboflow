// Coordinate-based dedup via `rfc dedup`. CONCAT_SORT_BED already emits a
// coordinate-sorted BED, so the input is passed straight through.
process RFC_DEDUP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bed)

    // `prefix` is resolved independently in `output:`, `script:` and `stub:`, and is
    // always `def`-scoped. An undeclared `prefix = ...` lands in the session-global
    // script binding shared by every concurrent task, leaking one process's prefix
    // into another's output filenames.
    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id + '.genome'}.post_dedup.bed"), emit: bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.genome"
    """
    rfc dedup -i ${bed} -o ${prefix}.post_dedup.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.genome"
    """
    touch ${prefix}.post_dedup.bed
    """
}
