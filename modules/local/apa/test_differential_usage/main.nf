process APA_TEST_DIFFERENTIAL_USAGE {
    tag "test_differential_usage"
    label 'process_single'
    label 'process_r'

    conda "${projectDir}/envs/r.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-r.sif" : null
    publishDir "${params.outdir}/apa", mode: params.publish_dir_mode, pattern: "apa_stats.tsv"
    publishDir "${params.outdir}/apa", mode: params.publish_dir_mode, pattern: "model_summary.tsv"

    input:
    path apa_usage
    path cell_annotations

    output:
    path "apa_stats.tsv", emit: apa_stats
    path "model_summary.tsv", emit: model_summary

    script:
    """
    Rscript ${projectDir}/bin/test_differential_usage.R \\
        --apa-usage ${apa_usage} \\
        --cell-annotations ${cell_annotations} \\
        --out-stats apa_stats.tsv \\
        --out-summary model_summary.tsv
    """

    stub:
    """
    printf "contrast_id\\tgene_id\\tsite_id\\tgroup_level\\tlogFC\\tdelta_usage\\tpvalue\\tfdr\\ttest_class\\ndescriptive__cluster\\tgeneA\\tgeneA:chr1:100:+\\tcluster\\t0\\t0\\tNA\\tNA\\tdescriptive_only\\n" > apa_stats.tsv
    printf "metric\\tvalue\\nn_usage_rows\\t1\\nn_stats_rows\\t1\\n" > model_summary.tsv
    """
}
