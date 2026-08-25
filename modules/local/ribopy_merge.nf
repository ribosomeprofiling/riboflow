// Merge all per-sample .ribo files into a single all.ribo.
// Ports `ribopy_merge` (upstream RiboFlow). Runs after RIBOPY_CREATE (and
// RIBOPY_RNASEQ_SET when do_rnaseq) so each sample's .ribo is fully populated.
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
