process APA_QUANTIFY_SITE_USAGE {
    tag "quantify_site_usage"
    label 'process_medium'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/apa", mode: params.publish_dir_mode

    input:
    path coverage_files
    path site_catalog
    path group_map
    path cell_annotations

    output:
    path "apa_usage.tsv", emit: apa_usage
    path "site_usage_matrix.tsv", emit: site_usage_matrix
    path "gene_summary.tsv", emit: gene_summary

    script:
    def coverageArgs = coverage_files.collect { it.toString() }.join(' ')
    """
    python ${projectDir}/bin/quantify_site_usage.py \\
        --coverage-files ${coverageArgs} \\
        --site-catalog ${site_catalog} \\
        --group-map ${group_map} \\
        --cell-annotations ${cell_annotations} \\
        --out-usage apa_usage.tsv \\
        --out-matrix site_usage_matrix.tsv \\
        --out-gene-summary gene_summary.tsv
    """

    stub:
    """
    printf "group_level\\tgroup_id\\tgene_id\\tsite_id\\tread_count\\tumi_count\\tusage_fraction\\tpdui_like\\ncluster\\tcluster_1\\tgeneA\\tgeneA:chr1:100:+\\t1\\t1\\t1.000000\\t1.000000\\n" > apa_usage.tsv
    printf "site_id\\tcluster_1\\ngeneA:chr1:100:+\\t1.000000\\n" > site_usage_matrix.tsv
    printf "group_level\\tgroup_id\\tgene_id\\tn_sites\\ttotal_umi_count\\ncluster\\tcluster_1\\tgeneA\\t1\\t1\\n" > gene_summary.tsv
    """
}
