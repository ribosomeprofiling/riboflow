// Adapter trimming for the RNA-seq path, single- or paired-end. Outputs are named
// with the rnaseq.clipped_R1/R2 prefix so the storeDir stays separate from
// ribo-seq. meta.single_end == false means paired-end, where `reads` is staged as
// [R1, R2]; otherwise `reads` is one file. Outputs use -o and -p, so the report
// goes to stdout.
process CUTADAPT_CLIP_RNASEQ {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.clipped_R1.fastq.gz"),               emit: reads
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.clipped_R2.fastq.gz"), optional: true, emit: reads2
    tuple val(meta), path("${meta.id}.${meta.lane}.rnaseq.clipped.log"),                       emit: log

    script:
    def prefix = "${meta.id}.${meta.lane}"
    def args   = task.ext.args ?: (params.rnaseq?.clip_arguments ?: '')
    def is_pe  = meta.single_end == false
    if (is_pe) {
        """
        cutadapt --cores=${task.cpus} ${args} \\
            -o ${prefix}.rnaseq.clipped_R1.fastq.gz \\
            -p ${prefix}.rnaseq.clipped_R2.fastq.gz \\
            ${reads[0]} ${reads[1]} > ${prefix}.rnaseq.clipped.log 2>&1
        """
    } else {
        """
        cutadapt --cores=${task.cpus} ${args} \\
            -o ${prefix}.rnaseq.clipped_R1.fastq.gz ${reads} > ${prefix}.rnaseq.clipped.log 2>&1
        """
    }

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    def is_pe  = meta.single_end == false
    def r2_cmd = is_pe ? "echo | gzip -c > ${prefix}.rnaseq.clipped_R2.fastq.gz" : ''
    """
    echo | gzip -c > ${prefix}.rnaseq.clipped_R1.fastq.gz
    ${r2_cmd}
    touch ${prefix}.rnaseq.clipped.log
    """
}
