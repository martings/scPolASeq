process SCAPTURE_MANIFEST {
    tag "${meta.library_id}"
    label 'process_low'

    conda "${projectDir}/envs/scapture.yml"
    container "ghcr.io/martings/scapture:f056f02"

    input:
    tuple val(meta), path(bam), path(gtf), path(fasta), path(fai), path(prepared_gtf), path(gene_map), path(duplicate_report), path(pascall_metrics), path(reference_metrics), path(catalog_metrics), path(quant_manifest), path(manifest_helper)

    output:
    path("${prefix}.scapture.manifest.tsv"), emit: manifest
    path('versions.yml'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix   = task.ext.prefix ?: "${meta.library_id?.toString()}"
    def args = task.ext.args != null ? task.ext.args.toString().trim() : "--species human --read-length 90 --peak-width 400 --normality-seed 1 --threshold 0.5"
    if (!(prefix ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error "SCAPTURE_MANIFEST requires a safe library_id/prefix: ${prefix}"
    }
    """
    set -euo pipefail
    python ${manifest_helper} \\
      --library-id '${prefix}' --bam '${bam.name}' --gtf '${gtf.name}' \\
      --fasta '${fasta.name}' --fai '${fai.name}' \\
      --prepared-gtf '${prepared_gtf.name}' --gene-map '${gene_map.name}' \\
      --duplicate-report ${duplicate_report} \\
      --pascall-metrics ${pascall_metrics} \\
      --reference-metrics ${reference_metrics} \\
      --catalog-metrics ${catalog_metrics} \\
      --quant-manifest ${quant_manifest} \\
      ${args} \\
      --output '${prefix}.scapture.manifest.tsv'
    printf '"%s":\\n  python: "%s"\\n  manifest_helper_sha256: "%s"\\n' \\
      "${task.process}" \\
      "\$(python --version 2>&1 | awk '{print \$NF}')" \\
      "\$(sha256sum ${manifest_helper} | cut -d' ' -f1)" > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.library_id?.toString()}"
    """
    printf 'field\\tvalue\\nlibrary_id\\t%s\\nsource\\tSCAPTURE\\nmode\\tstub\\n' '${prefix}' > '${prefix}.scapture.manifest.tsv'
    printf '"%s":\\n  python: "3.8.19"\\n' "${task.process}" > versions.yml
    """
}
