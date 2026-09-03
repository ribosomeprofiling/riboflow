// Merge all per-sample .ribo files into a single all.ribo. Runs after
// RIBOPY_CREATE, and after RIBOPY_RNASEQ_SET when do_rnaseq is set, so every
// sample's .ribo is complete first.
process RIBOPY_MERGE {
    tag "all"

    input:
    path(ribo_files)   // collected list of all per-sample .ribo files

    output:
    path "all.ribo", emit: ribo

    script:
    """
    ribopy merge all.ribo ${ribo_files}
    """

    stub:
    """
    touch all.ribo
    """
}
