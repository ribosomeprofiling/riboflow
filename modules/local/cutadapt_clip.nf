// Adapter trimming. The output is written with -o so cutadapt compresses it with
// its own worker threads. In that mode the report goes to stdout, hence the
// `> log 2>&1`.
process CUTADAPT_CLIP {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.clipped.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.${meta.lane}.clipped.log"),      emit: log

    script:
    def prefix = "${meta.id}.${meta.lane}"
    def args   = task.ext.args ?: (params.clip_arguments ?: '')
    """
    cutadapt --cores=${task.cpus} ${args} \\
        -o ${prefix}.clipped.fastq.gz ${fastq} > ${prefix}.clipped.log 2>&1
    """

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    """
    echo | gzip -c > ${prefix}.clipped.fastq.gz
    touch ${prefix}.clipped.log
    """
}
