// Concatenate the per-lane BEDs and coordinate-sort them. ext.prefix sets the
// output basename.
//
// LC_ALL=C keeps the sort order locale-independent, which is what rfc dedup and
// the downstream tools expect.
process CONCAT_SORT_BED {
    tag "${meta.id}"

    input:
    tuple val(meta), path(beds)

    // `prefix` is resolved independently in `output:`, `script:` and `stub:`, and is
    // always `def`-scoped. An undeclared `prefix = ...` lands in the session-global
    // script binding shared by every concurrent task, leaking one process's prefix
    // into another's output filenames.
    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id + '.genome.merged.pre_dedup'}.bed"), emit: bed

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}.genome.merged.pre_dedup"
    def sort_mem = Utils.sort_mem_mb(task)
    """
    cat ${beds} | LC_ALL=C sort --parallel=${task.cpus} -S ${sort_mem}M -T . -k1,1 -k2,2n -k3,3n > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.genome.merged.pre_dedup"
    """
    touch ${prefix}.bed
    """
}
