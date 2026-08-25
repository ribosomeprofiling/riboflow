// UMI extraction (umicollapse path only). Ports `extract_umi_via_umi_tools`
// (RiboFlow.groovy:280-298). SE/PE aware: for PE the UMI is extracted from both
// mates via --read2-in/--read2-out. The R1 output name is identical in SE and PE
// so the `reads` glob is stable (ribo-seq, always SE, keeps its filename).
process UMITOOLS_EXTRACT {
    tag "${meta.id}.${meta.lane}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.${meta.lane}.umi_extracted.fastq.gz"),    emit: reads
    tuple val(meta), path("${meta.id}.${meta.lane}.umi_extracted_R2.fastq.gz"), optional: true, emit: reads2
    tuple val(meta), path("${meta.id}.${meta.lane}.umi_extracted.log"),         emit: log

    script:
    def prefix = "${meta.id}.${meta.lane}"
    def args   = task.ext.args ?: (params.umi_tools_extract_arguments ?: '')
    def is_pe  = meta.single_end == false
    if (is_pe) {
        """
        umi_tools extract -I ${reads[0]} --read2-in ${reads[1]} \\
            -S ${prefix}.umi_extracted.fastq.gz \\
            --read2-out ${prefix}.umi_extracted_R2.fastq.gz \\
            -L ${prefix}.umi_extracted.log \\
            ${args}
        """
    } else {
        """
        umi_tools extract -I ${reads} -S ${prefix}.umi_extracted.fastq.gz \\
            -L ${prefix}.umi_extracted.log \\
            ${args}
        """
    }

    stub:
    def prefix = "${meta.id}.${meta.lane}"
    def is_pe  = meta.single_end == false
    def r2_cmd = is_pe ? "echo | gzip -c > ${prefix}.umi_extracted_R2.fastq.gz" : ''
    """
    echo | gzip -c > ${prefix}.umi_extracted.fastq.gz
    ${r2_cmd}
    touch ${prefix}.umi_extracted.log
    """
}
