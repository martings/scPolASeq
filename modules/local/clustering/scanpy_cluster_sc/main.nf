process SCANPY_CLUSTER_SC {
    tag "$meta.library_id"
    label 'process_medium'
    label 'process_python'

    conda "${projectDir}/envs/clustering.yml"
    publishDir "${params.outdir}/labels", mode: params.publish_dir_mode, pattern: "*.cell_annotations.tsv"
    publishDir "${params.outdir}/labels", mode: params.publish_dir_mode, pattern: "*.cluster_report.tsv"
    publishDir "${params.outdir}/labels", mode: params.publish_dir_mode, pattern: "*.clustered.h5ad"

    input:
    tuple val(meta), path(matrix_dir), path(cell_metadata), val(enable_internal_clustering)

    output:
    tuple val(meta), path("${meta.library_id}.cell_annotations.tsv"), emit: cell_annotations
    tuple val(meta), path("${meta.library_id}.cluster_report.tsv"), emit: report
    tuple val(meta), path("${meta.library_id}.clustered.h5ad"), emit: h5ad

    script:
    """
    python ${projectDir}/bin/scanpy_cluster_sc.py \\
        --matrix-dir ${matrix_dir} \\
        --cell-metadata ${cell_metadata} \\
        --sample-id ${meta.sample_id} \\
        --library-id ${meta.library_id} \\
        --enable-internal-clustering ${enable_internal_clustering} \\
        --out-annotations ${meta.library_id}.cell_annotations.tsv \\
        --out-report ${meta.library_id}.cluster_report.tsv \\
        --out-h5ad ${meta.library_id}.clustered.h5ad
    """

    stub:
    """
    printf "sample_id\\tbarcode_raw\\tbarcode_corrected\\tcell_id\\tcluster_id\\tcell_type\\tcondition\\tbatch\\tlabel_source\\n${meta.sample_id}\\t${meta.sample_id}-CELL001\\t${meta.sample_id}-CELL001\\t${meta.sample_id}:${meta.sample_id}-CELL001\\tcluster_1\\tunlabeled\\t\\t\\tstub\\n" > ${meta.library_id}.cell_annotations.tsv
    printf "sample_id\\tlibrary_id\\tn_cells\\tclustering_mode\\n${meta.sample_id}\\t${meta.library_id}\\t1\\tstub\\n" > ${meta.library_id}.cluster_report.tsv
    echo "{}" > ${meta.library_id}.clustered.h5ad
    """
}
