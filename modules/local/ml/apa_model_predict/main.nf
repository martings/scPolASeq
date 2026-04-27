// Stage 9 — Score candidate APA events using the trained model.
// Outputs per-site-per-group APA probability scores.
process APA_MODEL_PREDICT {
    tag "apa_model_predict"
    label 'process_medium'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/apa_calls", mode: params.publish_dir_mode

    input:
    path feature_table
    path trained_model
    val  group_level
    val  enable_single_cell_projection

    output:
    path "scored_apa_events.tsv", emit: scored_events
    path "apa_predictions.log",   emit: log

    script:
    def scFlag = enable_single_cell_projection ? '--single-cell-projection' : ''
    """
    python3 ${projectDir}/bin/apa_model_predict.py \\
        --features     ${feature_table} \\
        --model        ${trained_model} \\
        --group-level  ${group_level}   \\
        ${scFlag}                        \\
        --out-scored   scored_apa_events.tsv \\
        --log          apa_predictions.log
    """

    stub:
    """
    printf "site_id\tgroup_id\tgroup_level\tapa_score\tis_apa\n"       > scored_apa_events.tsv
    printf "geneA:chr1:2500:+\tcluster_1\tcluster\t0.91\tTRUE\n"     >> scored_apa_events.tsv
    printf "geneA:chr1:2500:+\tcluster_2\tcluster\t0.43\tFALSE\n"    >> scored_apa_events.tsv
    echo "APA prediction complete (stub)" > apa_predictions.log
    """
}
