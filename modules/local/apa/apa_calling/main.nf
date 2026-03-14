// Stage 7 — Statistical APA detection (DaPars-like PDUI metric).
// Compares proximal vs distal site usage between groups; outputs significant events.
process APA_CALLING {
    tag "apa_calling"
    label 'process_medium'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/apa_calls", mode: params.publish_dir_mode

    input:
    path feature_table
    path site_catalog
    path cell_annotations
    val  min_coverage
    val  min_pdui_delta

    output:
    path "apa_events.tsv",        emit: apa_events
    path "pdui_usage_matrix.tsv", emit: pdui_matrix
    path "apa_calling.log",       emit: log

    script:
    """
    python ${projectDir}/bin/apa_calling.py \\
        --features         ${feature_table}    \\
        --site-catalog     ${site_catalog}     \\
        --cell-annotations ${cell_annotations} \\
        --min-coverage     ${min_coverage}     \\
        --min-pdui-delta   ${min_pdui_delta}   \\
        --out-events       apa_events.tsv       \\
        --out-pdui         pdui_usage_matrix.tsv \\
        --log              apa_calling.log
    """

    stub:
    """
    printf "gene_id\tsite_id\tgroup_a\tgroup_b\tpdui_a\tpdui_b\tdelta_pdui\tp_value\tadj_p_value\tis_significant\n" \
        > apa_events.tsv
    printf "geneA\tgeneA:chr1:2500:+\tcluster_1\tcluster_2\t0.75\t0.40\t0.35\t0.001\t0.01\tTRUE\n" \
        >> apa_events.tsv
    printf "gene_id\tsite_id\tcluster_1\tcluster_2\n" > pdui_usage_matrix.tsv
    printf "geneA\tgeneA:chr1:2500:+\t0.75\t0.40\n"   >> pdui_usage_matrix.tsv
    echo "APA calling complete (stub)" > apa_calling.log
    """
}
