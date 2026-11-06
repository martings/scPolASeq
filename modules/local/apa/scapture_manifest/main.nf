process SCAPTURE_MANIFEST {
    tag "${meta.library_id}"
    label 'process_low'

    conda "${projectDir}/envs/scapture.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-scapture.sif" : null

    input:
    tuple val(meta), path(bam), path(gtf), path(fasta), path(fai), path(prepared_gtf), path(gene_map), path(duplicate_report), path(pascall_metrics), path(reference_metrics), path(catalog_metrics), path(quant_manifest), path(manifest_helper)

    output:
    path("${meta.library_id}.scapture.manifest.tsv"), emit: manifest
    path('versions.yml'), emit: versions

    script:
    def prefix = meta.library_id?.toString()
    def threshold = params.scapture_deeppass_threshold != null ? params.scapture_deeppass_threshold : 0.5
    def normality_seed = params.scapture_normality_seed == null ? 1 : params.scapture_normality_seed
    if (!(prefix ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error "SCAPTURE_MANIFEST requires a safe library_id: ${prefix}"
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
      --species '${params.scapture_species ?: 'human'}' \\
      --read-length '${params.scapture_read_length ?: 90}' \\
      --peak-width '${params.scapture_peak_width ?: 400}' \\
      --normality-seed '${normality_seed}' \\
      --threshold '${threshold}' \\
      --output '${prefix}.scapture.manifest.tsv'
    printf '"SCAPTURE_MANIFEST":\\n  python: "%s"\\n  manifest_helper_sha256: "%s"\\n' \\
      "\$(python --version 2>&1 | awk '{print \$NF}')" \\
      "\$(sha256sum ${manifest_helper} | cut -d' ' -f1)" > versions.yml
    """
}
