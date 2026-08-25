// Strand-specific bigWigs. Ports `genome_create_strand_specific_bigwigs`
// (RiboFlow.groovy:1391-1430). deepTools --filterRNAstrand assumes reverse-
// stranded libraries, so for forward-stranded ribo-seq (strand F/FR) the
// plus/minus filters are swapped.
process DEEPTOOLS_BAMCOVERAGE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)
    val assay   // filename label: 'ribo' for ribo-seq, 'rna' for RNA-seq

    output:
    tuple val(meta), path("${meta.id}.${assay}.plus.bigWig"),
                     path("${meta.id}.${assay}.minus.bigWig"), emit: bigwig

    when:
    // Plain string test: lib/ classes are not visible in a `when:` block.
    (params.do_bigwig?.toString()?.toLowerCase() in ['true', 'yes', 'on', '1'])

    script:
    def strand_arg = meta.strand ?: 'F'
    def fwd        = (strand_arg == 'F' || strand_arg == 'FR')
    // deepTools --filterRNAstrand assumes reverse-stranded libraries, so for
    // forward-stranded ribo-seq (strand F/FR) the plus/minus filters are swapped.
    def plus_filter  = fwd ? 'reverse' : 'forward'
    def minus_filter = fwd ? 'forward' : 'reverse'
    // PE: extend each pair to its full fragment so coverage reflects fragments.
    def extend     = (meta.single_end == false) ? '--extendReads' : ''
    // The two strands are independent: run both, each with half the threads.
    def bw_threads = Math.max(1, Math.min(task.cpus as int, 8).intdiv(2))
    def common     = "${extend} --binSize 1 -p ${bw_threads} --minMappingQuality 0 --outFileFormat bigwig"
    def block = Utils.parallel_block([
        "bamCoverage -b ${bam} -o ${meta.id}.${assay}.plus.bigWig  --filterRNAstrand ${plus_filter}  ${common}",
        "bamCoverage -b ${bam} -o ${meta.id}.${assay}.minus.bigWig --filterRNAstrand ${minus_filter} ${common}",
    ])
    """
    ${block}
    """

    stub:
    """
    touch ${meta.id}.${assay}.plus.bigWig ${meta.id}.${assay}.minus.bigWig
    """
}
