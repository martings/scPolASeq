process INTERNAL_PRIMING_FILTER {
    tag "internal_priming_filter"
    label 'process_single'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/apa", mode: params.publish_dir_mode, pattern: "site_catalog.filtered.tsv"
    publishDir "${params.outdir}/apa", mode: params.publish_dir_mode, pattern: "priming_filter_metrics.tsv"

    input:
    tuple path(site_catalog), path(genome_fasta), path(priming_blacklist)

    output:
    path "site_catalog.filtered.tsv", emit: site_catalog
    path "priming_filter_metrics.tsv", emit: metrics

    script:
    """
    python ${projectDir}/bin/internal_priming_filter.py \\
        --site-catalog ${site_catalog} \\
        --genome-fasta ${genome_fasta} \\
        --priming-blacklist ${priming_blacklist} \\
        --out-tsv site_catalog.filtered.tsv \\
        --out-metrics priming_filter_metrics.tsv
    """

    stub:
    """
    cp ${site_catalog} site_catalog.filtered.tsv
    printf "total_sites\\tflagged_sites\\n1\\t0\\n" > priming_filter_metrics.tsv
    """
}
