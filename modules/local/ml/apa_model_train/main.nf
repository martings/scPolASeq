// Stage 8 — Train ML model to classify APA events from coverage features.
// Supports: random_forest (default), xgboost, logistic_regression.
process APA_MODEL_TRAIN {
    tag "apa_model_train:${model_type}"
    label 'process_medium'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/models", mode: params.publish_dir_mode

    input:
    path feature_table
    path apa_events
    val  model_type     // 'random_forest' | 'xgboost' | 'logistic_regression'
    val  group_level    // 'cluster' | 'cell_type'

    output:
    path "apa_model.pkl",          emit: trained_model
    path "feature_importance.tsv", emit: feature_importance
    path "model_metrics.tsv",      emit: metrics
    path "model_train.log",        emit: log

    script:
    """
    python3 ${projectDir}/bin/apa_model_train.py \\
        --features       ${feature_table} \\
        --apa-events     ${apa_events}    \\
        --model-type     ${model_type}    \\
        --group-level    ${group_level}   \\
        --out-model      apa_model.pkl    \\
        --out-importance feature_importance.tsv \\
        --out-metrics    model_metrics.tsv      \\
        --log            model_train.log
    """

    stub:
    """
    touch apa_model.pkl
    printf "feature\timportance\n"                   > feature_importance.tsv
    printf "coverage_at_site\t0.35\n"               >> feature_importance.tsv
    printf "proximal_distal_ratio\t0.28\n"          >> feature_importance.tsv
    printf "read_end_density\t0.20\n"               >> feature_importance.tsv
    printf "upstream_cov\t0.17\n"                   >> feature_importance.tsv
    printf "metric\tvalue\n"                         > model_metrics.tsv
    printf "auc_roc_cv_mean\t0.87\n"               >> model_metrics.tsv
    printf "accuracy_train\t0.82\n"                >> model_metrics.tsv
    printf "model_type\t${model_type}\n"            >> model_metrics.tsv
    echo "Model training complete (stub): ${model_type}" > model_train.log
    """
}
