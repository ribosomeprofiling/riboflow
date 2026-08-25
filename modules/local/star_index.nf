// STAR genome index generation (genomeGenerate mode).
// Triggered when the user provides genome_fasta + gtf instead of a pre-built
// index directory. Output directory is piped directly into STAR_ALIGN via
// ch_genome_index = STAR_INDEX.out.index.first() in workflows/riboflow.nf.
process STAR_INDEX {
    tag "star_index"

    input:
    path fasta
    path gtf

    output:
    path "star_index", emit: index

    script:
    def overhang = task.ext.sjdb_overhang ?: (params.star?.sjdb_overhang ?: 100)
    def extra    = task.ext.args          ?: (params.star?.index_args ?: '')
    """
    mkdir -p star_index

    # --genomeSAindexNbases must be scaled to the genome or genomeGenerate either dies
    # (too large for a small reference) or wastes memory. STAR's documented rule is
    # min(14, log2(genomeLength)/2 - 1), and STAR does NOT apply it itself. Derive it
    # from the FASTA so swapping the reference — chrM to chr1 to a whole assembly —
    # needs no config change. An explicit --genomeSAindexNbases in star.index_args
    # still wins.
    AUTO_SA_INDEX_NBASES=''
    case '${extra}' in
        *genomeSAindexNbases*)
            echo "STAR_INDEX: --genomeSAindexNbases supplied via index_args, not deriving one"
            ;;
        *)
            GENOME_LEN=\$(awk '!/^>/ { n += length(\$0) } END { print n+0 }' ${fasta})
            [ "\$GENOME_LEN" -gt 0 ] || { echo "STAR_INDEX: no sequence found in ${fasta}" >&2; exit 1; }
            NBASES=\$(awk -v L="\$GENOME_LEN" 'BEGIN {
                v = int(log(L) / log(2) / 2 - 1)
                if (v > 14) v = 14
                if (v < 2)  v = 2
                print v
            }')
            echo "STAR_INDEX: genome length \$GENOME_LEN bp -> --genomeSAindexNbases \$NBASES"
            AUTO_SA_INDEX_NBASES="--genomeSAindexNbases \$NBASES"
            ;;
    esac

    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang ${overhang} \\
        \$AUTO_SA_INDEX_NBASES \\
        ${extra}
    """

    stub:
    """
    mkdir -p star_index
    touch star_index/SA star_index/SAindex star_index/Genome star_index/chrNameLength.txt
    """
}
