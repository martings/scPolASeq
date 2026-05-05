process HARMONIZE_CELL_METADATA {
    tag "harmonize_cell_metadata"
    label 'process_single'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/labels", mode: params.publish_dir_mode, pattern: "cell_metadata.harmonized.tsv"
    publishDir "${params.outdir}/provenance", mode: params.publish_dir_mode, pattern: "cell_metadata_manifest.json"

    input:
    path cell_metadata
    path cluster_assignments
    path cell_type_labels

    output:
    path "cell_metadata.harmonized.tsv", emit: cell_metadata
    path "cell_metadata_manifest.json", emit: manifest

    script:
    """
    python ${projectDir}/bin/harmonize_cell_metadata.py \\
        --cell-metadata ${cell_metadata} \\
        --cluster-assignments ${cluster_assignments} \\
        --cell-type-labels ${cell_type_labels} \\
        --out cell_metadata.harmonized.tsv \\
        --manifest cell_metadata_manifest.json
    """

    stub:
    """
    printf "sample_id\\tlibrary_id\\tbarcode_raw\\tbarcode_corrected\\tcell_id\\tcluster_id\\tcell_type\\tcondition\\tbatch\\tlabel_source\\n" > cell_metadata.harmonized.tsv
    echo "{}" > cell_metadata_manifest.json
    """
}
