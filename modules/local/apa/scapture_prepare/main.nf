process SCAPTURE_PREPARE {
    tag "${meta.library_id}"
    label 'process_low'

    conda "${projectDir}/envs/scapture.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-scapture.sif" : null

    input:
    tuple val(meta), path(gtf), path(prepare_script)

    output:
    tuple val(meta), path("${meta.library_id}.scapture.annotation.gtf"), emit: gtf
    tuple val(meta), path("${meta.library_id}.scapture.gene_name_map.tsv"), emit: gene_name_map
    tuple val(meta), path("${meta.library_id}.annotation.duplicate_gene_names.tsv"), emit: duplicate_gene_names
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.library_id
    """
    set -euo pipefail
    python ${prepare_script} \
        --gtf ${gtf} \
        --out-gtf ${prefix}.scapture.annotation.gtf \
        --out-map ${prefix}.scapture.gene_name_map.tsv \
        --out-duplicate-report ${prefix}.annotation.duplicate_gene_names.tsv
    test -s ${prefix}.scapture.annotation.gtf
    test -s ${prefix}.scapture.gene_name_map.tsv
    printf '"SCAPTURE_PREPARE":\n  python: "%s"\n' "\$(python --version 2>&1 | awk '{print \$2}')" > versions.yml
    """
}
