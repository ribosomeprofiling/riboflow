// Adapter trimming for the RNA-seq path. SE/PE aware. Mirrors cutadapt_clip.nf
// but always names outputs with the rnaseq.clipped_R1/R2 prefix so the RNA-seq
// storeDir never collides with ribo-seq. meta.single_end == false => paired-end:
// `reads` is staged as a 2-element list [R1, R2]; otherwise `reads` is one file.
// Outputs are written with -o/-p (cutadapt's own parallel gzip); the report goes
// to stdout in that mode.
process CUTADAPT_CLIP_RNASEQ {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.rnaseq.clipped_R1.fastq.gz"),               emit: reads
    tuple val(meta), path("${prefix}.rnaseq.clipped_R2.fastq.gz"), optional: true, emit: reads2
    tuple val(meta), path("${prefix}.rnaseq.clipped.log"),                       emit: log

    script:
    prefix     = "${meta.id}.${meta.lane}"
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
    prefix     = "${meta.id}.${meta.lane}"
    def is_pe  = meta.single_end == false
    def r2_cmd = is_pe ? "echo | gzip -c > ${prefix}.rnaseq.clipped_R2.fastq.gz" : ''
    """
    echo | gzip -c > ${prefix}.rnaseq.clipped_R1.fastq.gz
    ${r2_cmd}
    touch ${prefix}.rnaseq.clipped.log
    """
}
