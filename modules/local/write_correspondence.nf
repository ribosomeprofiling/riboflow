// Publish the sample and lane to FASTQ correspondence table.
process WRITE_CORRESPONDENCE {

    input:
    path(correspondence)

    output:
    path('index_fastq_correspondence.txt'), emit: txt

    script:
    """
    cat ${correspondence} > index_fastq_correspondence.txt
    """

    stub:
    """
    touch index_fastq_correspondence.txt
    """
}
