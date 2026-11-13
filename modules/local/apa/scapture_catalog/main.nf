process SCAPTURE_CATALOG {
    tag "${meta.library_id}"
    label 'process_low'

    conda "${projectDir}/envs/scapture.yml"
    container "ghcr.io/martings/scapture:f056f02"

    input:
    tuple val(meta), path(exonic), path(intronic), path(extended), path(catalog_helper)

    output:
    tuple val(meta), path("${prefix}.scapture.site_catalog.tsv"), emit: catalog
    tuple val(meta), path("${prefix}.scapture.catalog.metrics.tsv"), emit: metrics
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.library_id?.toString()}"
    if (!(prefix ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error "SCAPTURE_CATALOG requires a safe library_id/prefix: ${prefix}"
    }
    def args        = task.ext.args != null ? task.ext.args.toString().trim() : '--threshold 0.5'
    def thr_match   = (args =~ /(?:--threshold|-t)\s+([0-9]+(?:\.[0-9]+)?)/)
    def thr_str     = thr_match.find() ? thr_match.group(1) : '0.5'
    """
    set -euo pipefail
    python ${catalog_helper} \
        --exonic ${exonic} --intronic ${intronic} --extended ${extended} \
        ${args} --out-tsv ${prefix}.scapture.site_catalog.tsv
    rows=\$(awk 'END { print (NR > 0 ? NR - 1 : 0) }' ${prefix}.scapture.site_catalog.tsv)
    test "\$rows" -gt 0 || { echo 'ERROR: empty SCAPTURE site catalog' >&2; exit 1; }
    printf 'metric\tvalue\nsite_catalog_rows\t%s\ndeeppass_threshold\t%s\n' \
        "\$rows" '${thr_str}' > ${prefix}.scapture.catalog.metrics.tsv
    printf '"%s":\n  python: "%s"\n  catalog_helper_sha256: "%s"\n' \
        "${task.process}" "\$(python --version 2>&1 | awk '{print \$NF}')" \
        "\$(sha256sum ${catalog_helper} | cut -d' ' -f1)" > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.library_id?.toString()}"
    """
    printf 'site_id\tgene_id\tchrom\tstart\tend\tstrand\tsite_class\tsite_source\tpriming_flag\n' > ${prefix}.scapture.site_catalog.tsv
    printf 'chr1:120:+\tG1\tchr1\t100\t120\t+\texonic\tSCAPTURE\tfalse\n' >> ${prefix}.scapture.site_catalog.tsv
    printf 'metric\tvalue\nsite_catalog_rows\t1\ndeeppass_threshold\t0.5\n' > ${prefix}.scapture.catalog.metrics.tsv
    printf '"%s":\n  python: "3.8.19"\n' "${task.process}" > versions.yml
    """
}
