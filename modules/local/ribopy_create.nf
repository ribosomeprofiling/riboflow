// Creates a per-sample .ribo file from a merged post-dedup BED.
// BED is stripped to 6 columns before passing to ribopy: the position-dedup
// path appends a 7th sample.lane column for dedup; ribopy only wants cols 1-6.
// All tuning flags come from params.ribo.*.
//
// ribometa and the global expmeta are staged as inputs rather than passed as host
// paths, so they are visible inside the container.
process RIBOPY_CREATE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(dedup_bed), path(expmeta_file)
    path regions_bed
    path lengths_tsv
    path(ribometa_file,   stageAs: 'ribometa/*')       // [] when params.ribo.ribometa is unset
    path(global_expmeta,  stageAs: 'global_expmeta/*') // [] when params.ribo.expmeta is unset

    output:
    tuple val(meta), path("${meta.id}.ribo"), emit: ribo

    script:
    // `!= null`, not Elvis: 0 is a legitimate span / radius / minimum length.
    def r            = (params.ribo ?: [:])
    def ref_name     = r.ref_name ?: 'appris_human'
    def radius       = (r.metagene_radius  != null) ? r.metagene_radius  : 50
    def left_span    = (r.left_span        != null) ? r.left_span        : 35
    def right_span   = (r.right_span       != null) ? r.right_span       : 15
    def len_min      = (r.read_length?.min != null) ? r.read_length.min  : 15
    def len_max      = (r.read_length?.max != null) ? r.read_length.max  : 35
    def nocov_flag   = Utils.as_bool(r.coverage, true) ? '' : '--nocoverage'
    // The per-sample expmeta wins; otherwise fall back to the global one.
    def expmeta_arg  = expmeta_file   ? "--expmeta ${expmeta_file}"
                     : global_expmeta ? "--expmeta ${global_expmeta}"
                     : ''
    def ribometa_arg = ribometa_file  ? "--ribometa ${ribometa_file}" : ''
    """
    ribopy create \\
        --name ${meta.id} \\
        --alignmentfile <(cut -f1-6 ${dedup_bed}) \\
        --reference ${ref_name} \\
        --lengths ${lengths_tsv} \\
        --annotation ${regions_bed} \\
        --metageneradius ${radius} \\
        -l ${left_span} -r ${right_span} \\
        --lengthmin ${len_min} --lengthmax ${len_max} \\
        ${nocov_flag} \\
        ${expmeta_arg} \\
        ${ribometa_arg} \\
        -n ${task.cpus} \\
        ${meta.id}.ribo
    """

    stub:
    """
    touch ${meta.id}.ribo
    """
}
