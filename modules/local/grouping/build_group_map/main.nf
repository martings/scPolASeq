process BUILD_GROUP_MAP {
    tag "build_group_map"
    label 'process_single'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/labels", mode: params.publish_dir_mode

    input:
    path annotation_files
    val grouping_levels

    output:
    path "cell_annotations.tsv", emit: cell_annotations
    path "group_map.tsv", emit: group_map
    path "group_summary.tsv", emit: group_summary

    script:
    def annotationArgs = annotation_files.collect { it.toString() }.join(' ')
    """
    python3 ${projectDir}/bin/build_group_map.py \\
        --annotation-files ${annotationArgs} \\
        --grouping-levels "${grouping_levels}" \\
        --out-cell-annotations cell_annotations.tsv \\
        --out-group-map group_map.tsv \\
        --out-group-summary group_summary.tsv
    """

    stub:
    """
    printf "sample_id\\tbarcode_raw\\tbarcode_corrected\\tcell_id\\tcluster_id\\tcell_type\\tcondition\\tbatch\\tlabel_source\\n" \\
        > cell_annotations.tsv
    printf "pbmc_1k_v3\\tpbmc_1k_v3_L001-CELL001\\tpbmc_1k_v3_L001-CELL001\\tpbmc_1k_v3:pbmc_1k_v3_L001-CELL001\\tcluster_1\\tunlabeled\\t\\t\\tstub\\n" \\
        >> cell_annotations.tsv
    printf "pbmc_1k_v3\\tpbmc_1k_v3_L002-CELL001\\tpbmc_1k_v3_L002-CELL001\\tpbmc_1k_v3:pbmc_1k_v3_L002-CELL001\\tcluster_2\\tunlabeled\\t\\t\\tstub\\n" \\
        >> cell_annotations.tsv
    printf "sample_id\\tbarcode_corrected\\tgroup_level\\tgroup_id\\n" \\
        > group_map.tsv
    printf "pbmc_1k_v3\\tpbmc_1k_v3_L001-CELL001\\tcluster\\tcluster_1\\n" \\
        >> group_map.tsv
    printf "pbmc_1k_v3\\tpbmc_1k_v3_L002-CELL001\\tcluster\\tcluster_2\\n" \\
        >> group_map.tsv
    printf "group_level\\tgroup_id\\tn_barcodes\\n" \\
        > group_summary.tsv
    printf "cluster\\tcluster_1\\t1\\n" >> group_summary.tsv
    printf "cluster\\tcluster_2\\t1\\n" >> group_summary.tsv
    """
}
