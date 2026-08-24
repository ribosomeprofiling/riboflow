// Parametrized FastQC, gated on params.do_fastqc. One module replaces
// raw_fastqc / clipped_fastqc / genome_aligned_fastqc / genome_unaligned_fastqc
// (RiboFlow.groovy:233,312,494,523). `stage` (via ext.prefix) sets the output
// basename; the >20-byte guard skips empty FASTQs (matches the genome variants).
// PE mates are handed to ONE fastqc call (-t = files in parallel).
process FASTQC {
    tag "${prefix}"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*_fastqc.html"), path("*_fastqc.zip"), emit: report

    when:
    // Plain string test: lib/ classes are not visible in a `when:` block.
    (params.do_fastqc?.toString()?.toLowerCase() in ['true', 'yes', 'on', '1'])

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.lane}"
    def reads = (fastq instanceof List) ? fastq : [fastq]
    def names = reads.size() > 1 ? reads.indices.collect { "${prefix}_R${it + 1}" } : [prefix]
    def links = [reads, names].transpose().collect { f, nm ->
        // Rename so the report is named after this stage — but skip when the input
        // basename already equals the target (ln -sf X X = broken self-link).
        "[ \"${f}\" != \"${nm}.fastq.gz\" ] && ln -sf ${f} ${nm}.fastq.gz || true"
    }.join('\n    ')
    """
    ${links}
    todo=()
    for nm in ${names.join(' ')}; do
        if [ \$(stat -L -c%s \${nm}.fastq.gz) -gt 20 ]; then
            todo+=(\${nm}.fastq.gz)
        else
            echo "\${nm}.fastq.gz is empty, skipping FastQC"
            touch \${nm}_fastqc.html \${nm}_fastqc.zip
        fi
    done
    if [ \${#todo[@]} -gt 0 ]; then
        fastqc "\${todo[@]}" --outdir=\$PWD -t ${Math.min(task.cpus as int, names.size())}
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.lane}"
    def reads = (fastq instanceof List) ? fastq : [fastq]
    if (reads.size() > 1) {
        """
        touch ${prefix}_R1_fastqc.html ${prefix}_R1_fastqc.zip ${prefix}_R2_fastqc.html ${prefix}_R2_fastqc.zip
        """
    } else {
        """
        touch ${prefix}_fastqc.html ${prefix}_fastqc.zip
        """
    }
}
