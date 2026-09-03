// Split the final sample BAM into plus- and minus-strand BAM and BED files.
// Gated on do_strand_split.
process SPLIT_STRANDED_BAM {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.ribo.plus.bam"),  path("${meta.id}.ribo.plus.bam.bai"),
                     path("${meta.id}.ribo.minus.bam"), path("${meta.id}.ribo.minus.bam.bai"),
                     path("${meta.id}.ribo.plus.bed"),  path("${meta.id}.ribo.minus.bed"), emit: stranded

    when:
    // Plain string test: classes in lib/ are not visible inside a `when:` block.
    (params.do_strand_split?.toString()?.toLowerCase() in ['true', 'yes', 'on', '1'])

    script:
    def strand_arg = meta.strand ?: 'F'
    def fwd        = (strand_arg == 'F' || strand_arg == 'FR')
    def plus_flags  = fwd ? '-F 2064'       : '-f 16 -F 2048'
    def minus_flags = fwd ? '-f 16 -F 2048' : '-F 2064'
    def t           = Math.max(1, (task.cpus as int).intdiv(2))
    // The two strands are independent: run them side by side.
    def block = Utils.parallel_block([
        "( samtools view -@ ${t} -b ${plus_flags}  -o ${meta.id}.ribo.plus.bam  ${bam} && samtools index -@ ${t} ${meta.id}.ribo.plus.bam  && bamToBed -i ${meta.id}.ribo.plus.bam  > ${meta.id}.ribo.plus.bed )",
        "( samtools view -@ ${t} -b ${minus_flags} -o ${meta.id}.ribo.minus.bam ${bam} && samtools index -@ ${t} ${meta.id}.ribo.minus.bam && bamToBed -i ${meta.id}.ribo.minus.bam > ${meta.id}.ribo.minus.bed )",
    ])
    """
    ${block}
    """

    stub:
    """
    touch ${meta.id}.ribo.plus.bam ${meta.id}.ribo.plus.bam.bai
    touch ${meta.id}.ribo.minus.bam ${meta.id}.ribo.minus.bam.bai
    touch ${meta.id}.ribo.plus.bed ${meta.id}.ribo.minus.bed
    """
}
