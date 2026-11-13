process SCAPTURE_PREPARE {
    tag "${meta.library_id}"
    label 'process_low'

    conda "${projectDir}/envs/scapture.yml"
    container "ghcr.io/martings/scapture:f056f02"

    input:
    tuple val(meta), path(gtf), path(prepare_script)

    output:
    tuple val(meta), path("${prefix}.scapture.annotation.gtf"), emit: gtf
    tuple val(meta), path("${prefix}.scapture.gene_name_map.tsv"), emit: gene_name_map
    tuple val(meta), path("${prefix}.annotation.duplicate_gene_names.tsv"), emit: duplicate_gene_names
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.library_id}"
    """
    set -euo pipefail
    python ${prepare_script} \
        --gtf ${gtf} \
        --out-gtf ${prefix}.scapture.annotation.gtf \
        --out-map ${prefix}.scapture.gene_name_map.tsv \
        --out-duplicate-report ${prefix}.annotation.duplicate_gene_names.tsv
    test -s ${prefix}.scapture.annotation.gtf
    test -s ${prefix}.scapture.gene_name_map.tsv
    printf '"%s":\n  python: "%s"\n' "${task.process}" "\$(python --version 2>&1 | awk '{print \$2}')" > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.library_id}"
    """
    cp ${gtf} ${prefix}.scapture.annotation.gtf 2>/dev/null || : > ${prefix}.scapture.annotation.gtf
    printf 'gene_id\tscapture_gene_name\toriginal_gene_name\n' > ${prefix}.scapture.gene_name_map.tsv
    printf 'gene_count\toriginal_gene_name\n' > ${prefix}.annotation.duplicate_gene_names.tsv
    printf '"%s":\n  python: "3.8.19"\n' "${task.process}" > versions.yml
    """
}
