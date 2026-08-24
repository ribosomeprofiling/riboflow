// UMI-aware dedup. Ports `genome_deduplicate_umicollapse` (RiboFlow.groovy:1208-1246)
// and `transcriptome_deduplicate_umicollapse` (:916-941, no counts).
//
// ENV CONSOLIDATION: the original called a hand-shipped jar via a `java11`
// wrapper. We now use the bioconda `umicollapse` jar from the single conda env.
// The JVM heap (ext.jvm_opts, set per profile to scale with task.attempt) must fit
// the container; the fallback derives it from task.memory for the same reason.
process UMICOLLAPSE_DEDUP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.merged_dedup.total.count"),
                     optional: true, emit: total_count
    tuple val(meta), path("${meta.id}.merged_dedup.primary.count"),
                     path("${meta.id}.merged_dedup.secondary.count"),
                     path("${meta.id}.merged_dedup.unique.count"),
                     optional: true, emit: detail_counts

    script:
    prefix               = task.ext.prefix ?: "${meta.id}.dedup"
    def route            = task.ext.route ?: 'genome'
    def args             = task.ext.args ?: (params.umicollapse_arguments ?: '')
    def algo             = task.ext.algo ?: 'cc'
    def heap_gb          = task.memory ? Math.max(1, ((task.memory.toGiga() * 0.75) as int)) : 6
    def jvm_opts         = task.ext.jvm_opts ?: "-Xms512m -Xmx${heap_gb}g -Xss256m"
    def emit_counts      = task.ext.emit_counts ?: false
    def emit_full_counts = (task.ext.emit_full_counts != null) ? task.ext.emit_full_counts : !Utils.route_unique_only(params, route)
    def is_pe            = meta.single_end == false
    def paired_arg       = is_pe ? '--paired' : ''
    def counts           = emit_counts
        ? Utils.samtools_count_block(Utils.route_count_args(params, route), is_pe,
                                     "${prefix}.bam", "${meta.id}.merged_dedup", task.cpus as int,
                                     emit_full_counts, emit_full_counts)
        : ''
    """
    ulimit -s unlimited
    JAR_DIR=\$(dirname \$(realpath \$(which umicollapse)))
    java ${jvm_opts} -jar \${JAR_DIR}/umicollapse.jar bam \\
        -i ${bam} \\
        -o ${prefix}.bam \\
        --umi-sep "_" \\
        --algo ${algo} \\
        --merge mapqual \\
        --two-pass \\
        ${paired_arg} \\
        ${args}
    samtools index -@ ${task.cpus} ${prefix}.bam
    ${counts}
    """

    stub:
    prefix               = task.ext.prefix ?: "${meta.id}.dedup"
    def route            = task.ext.route ?: 'genome'
    def emit_counts      = task.ext.emit_counts ?: false
    def emit_full_counts = (task.ext.emit_full_counts != null) ? task.ext.emit_full_counts : !Utils.route_unique_only(params, route)
    def total_cmd  = emit_counts ? "echo 0 > ${meta.id}.merged_dedup.total.count" : ''
    def detail_cmd = (emit_counts && emit_full_counts) ? "echo 0 > ${meta.id}.merged_dedup.primary.count; echo 0 > ${meta.id}.merged_dedup.secondary.count; echo 0 > ${meta.id}.merged_dedup.unique.count" : ''
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    ${total_cmd}
    ${detail_cmd}
    """
}
