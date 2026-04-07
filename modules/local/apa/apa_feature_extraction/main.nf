// Stage 6 — APA feature extraction around candidate polyA sites.
// Extracts coverage-based and sequence-based features per site × group pair.
process APA_FEATURE_EXTRACTION {
    tag "apa_feature_extraction"
    label 'process_medium'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/apa_features", mode: params.publish_dir_mode

    input:
    path site_catalog
    path bedgraphs, stageAs: '?/*'  // collected: all *.bedGraph files, staged in unique subdirs
    path cell_annotations
    path known_polya

    output:
    path "apa_features.tsv",       emit: feature_table
    path "feature_extraction.log", emit: log

    script:
    """
    python ${projectDir}/bin/apa_feature_extraction.py \\
        --site-catalog     ${site_catalog}     \\
        --bedgraph-dir     .                   \\
        --bedgraph-glob    '**/*.bedGraph'      \\
        --cell-annotations ${cell_annotations} \\
        --known-polya      ${known_polya}       \\
        --out-features     apa_features.tsv     \\
        --log              feature_extraction.log
    """

    stub:
    """
    printf "site_id\tgroup_level\tgroup_id\tcoverage_at_site\tupstream_cov\tdownstream_cov\t" \
           "read_end_density\tproximal_distal_ratio\tdist_to_known_polya\tumi_support\tcluster_support\n" \
        > apa_features.tsv
    printf "geneA:chr1:1000:+\tcluster\tcluster_1\t50.0\t45.0\t10.0\t8.5\t0.82\t150\t12\t3\n" \
        >> apa_features.tsv
    printf "geneA:chr1:2500:+\tcluster\tcluster_1\t30.0\t28.0\t25.0\t5.1\t0.55\t0\t10\t3\n"  \
        >> apa_features.tsv
    printf "geneA:chr1:1000:+\tcluster\tcluster_2\t40.0\t35.0\t15.0\t7.2\t0.65\t150\t10\t3\n" \
        >> apa_features.tsv
    printf "geneA:chr1:2500:+\tcluster\tcluster_2\t55.0\t50.0\t8.0\t9.3\t0.88\t0\t14\t3\n"  \
        >> apa_features.tsv
    echo "Feature extraction complete (stub)" > feature_extraction.log
    """
}
