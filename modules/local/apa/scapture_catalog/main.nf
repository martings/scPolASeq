process SCAPTURE_CATALOG {
    tag "${meta.library_id}"
    label 'process_low'

    conda "${projectDir}/envs/scapture.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-scapture.sif" : null

    input:
    tuple val(meta), path(exonic), path(intronic), path(extended), path(catalog_helper)

    output:
    tuple val(meta), path("${meta.library_id}.scapture.site_catalog.tsv"), emit: catalog
    tuple val(meta), path("${meta.library_id}.scapture.catalog.metrics.tsv"), emit: metrics
    path 'versions.yml', emit: versions

    script:
    def prefix = meta.library_id?.toString()
    if (!(prefix ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error "SCAPTURE_CATALOG requires a safe library_id: ${prefix}"
    }
    def threshold = (params.scapture_deeppass_threshold != null ? params.scapture_deeppass_threshold : 0.5) as Double
    if (!Double.isFinite(threshold) || threshold < 0 || threshold > 1) {
        error "SCAPTURE_CATALOG requires a DeepPASS threshold between 0 and 1: ${threshold}"
    }
    """
    set -euo pipefail
    python ${catalog_helper} \
        --exonic ${exonic} --intronic ${intronic} --extended ${extended} \
        --threshold ${threshold} --out-tsv ${prefix}.scapture.site_catalog.tsv
    rows=\$(awk 'END { print (NR > 0 ? NR - 1 : 0) }' ${prefix}.scapture.site_catalog.tsv)
    test "\$rows" -gt 0 || { echo 'ERROR: empty SCAPTURE site catalog' >&2; exit 1; }
    printf 'metric\tvalue\nsite_catalog_rows\t%s\ndeeppass_threshold\t%s\n' \
        "\$rows" '${threshold}' > ${prefix}.scapture.catalog.metrics.tsv
    printf '"SCAPTURE_CATALOG":\n  python: "%s"\n  catalog_helper_sha256: "%s"\n' \
        "\$(python --version 2>&1 | awk '{print \$NF}')" \
        "\$(sha256sum ${catalog_helper} | cut -d' ' -f1)" > versions.yml
    """
}
