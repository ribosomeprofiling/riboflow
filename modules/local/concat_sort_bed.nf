// Concatenate per-lane BEDs and coordinate-sort. Ports the several
// `cat <beds> | sort -k1,1 -k2,2n -k3,3n` steps: merge_genome_bed_for_position_dedup
// (RiboFlow.groovy:1048-1064), concat_genome_qpass_bed_for_publish_none (:1356-1374),
// concat_genome_post_dedup_bed_umi (:1319-1335) and the transcriptome BED merge
// (:745-762). ext.prefix sets the output basename.
//
// LC_ALL=C: byte-order collation is locale-independent (reproducible between the
// container and a conda run) and what `rfc dedup` and downstream tools expect.
process CONCAT_SORT_BED {
    tag "${meta.id}"

    input:
    tuple val(meta), path(beds)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed

    script:
    prefix       = task.ext.prefix ?: "${meta.id}.genome.merged.pre_dedup"
    def sort_mem = Utils.sort_mem_mb(task)
    """
    cat ${beds} | LC_ALL=C sort --parallel=${task.cpus} -S ${sort_mem}M -T . -k1,1 -k2,2n -k3,3n > ${prefix}.bed
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.genome.merged.pre_dedup"
    """
    touch ${prefix}.bed
    """
}
