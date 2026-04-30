process SCANPY_CLUSTER_SC {
    tag "$meta.library_id"
    label 'process_medium'
    label 'process_python'

    conda "${projectDir}/envs/clustering.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-clustering.sif" : null
    publishDir "${params.outdir}/labels", mode: params.publish_dir_mode, pattern: "*.cell_annotations.tsv"
    publishDir "${params.outdir}/labels", mode: params.publish_dir_mode, pattern: "*.cluster_report.tsv"
    publishDir "${params.outdir}/labels", mode: params.publish_dir_mode, pattern: "*.clustered.h5ad"
    publishDir "${params.outdir}/labels", mode: params.publish_dir_mode, pattern: "*.clustered.*.png"

    input:
    tuple val(meta), path(matrix_dir), path(cell_metadata), val(enable_internal_clustering), path(celltypist_model, stageAs: 'celltypist_model_input'), path(reference_h5ad, stageAs: 'reference_h5ad_input'), val(reference_label_col)

    output:
    tuple val(meta), path("${meta.library_id}.cell_annotations.tsv"), emit: cell_annotations
    tuple val(meta), path("${meta.library_id}.cluster_report.tsv"), emit: report
    tuple val(meta), path("${meta.library_id}.clustered.h5ad"), emit: h5ad
    tuple val(meta), path("${meta.library_id}.clustered.*.png"), optional: true, emit: embedding_plots

    script:
    def ct_model_arg = (celltypist_model  && celltypist_model.name  != 'NO_FILE') ? "--celltypist-model celltypist_model_input"   : ''
    def ref_h5ad_arg = (reference_h5ad    && reference_h5ad.name    != 'NO_FILE') ? "--reference-h5ad reference_h5ad_input"       : ''
    def ref_col_arg  = reference_label_col ? "--reference-label-col ${reference_label_col}" : ''
    """
    mkdir -p \$PWD/.numba_cache
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache
    python ${projectDir}/bin/scanpy_cluster_sc.py \\
        --matrix-dir ${matrix_dir} \\
        --cell-metadata ${cell_metadata} \\
        --sample-id ${meta.sample_id} \\
        --library-id ${meta.library_id} \\
        --enable-internal-clustering ${enable_internal_clustering} \\
        ${ct_model_arg} \\
        ${ref_h5ad_arg} \\
        ${ref_col_arg} \\
        --out-annotations ${meta.library_id}.cell_annotations.tsv \\
        --out-report ${meta.library_id}.cluster_report.tsv \\
        --out-h5ad ${meta.library_id}.clustered.h5ad \\
        --out-plot-prefix ${meta.library_id}.clustered
    """

    stub:
    """
    printf "sample_id\\tlibrary_id\\tbarcode_raw\\tbarcode_corrected\\tcell_id\\tcluster_id\\tcell_type\\tcondition\\tbatch\\tlabel_source\\n${meta.sample_id}\\t${meta.library_id}\\t${meta.library_id}-CELL001\\t${meta.library_id}-CELL001\\t${meta.sample_id}:${meta.library_id}:${meta.library_id}-CELL001\\t${meta.library_id.contains('L001') ? 'cluster_1' : 'cluster_2'}\\tunlabeled\\t\\t\\tstub\\n" > ${meta.library_id}.cell_annotations.tsv
    printf "sample_id\\tlibrary_id\\tn_cells\\tclustering_mode\\n${meta.sample_id}\\t${meta.library_id}\\t1\\tstub\\n" > ${meta.library_id}.cluster_report.tsv
    echo "{}" > ${meta.library_id}.clustered.h5ad
    touch ${meta.library_id}.clustered.pca.png
    touch ${meta.library_id}.clustered.umap.png
    """
}
