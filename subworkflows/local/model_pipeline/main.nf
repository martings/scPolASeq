// Stages 8 + 9 subworkflow — ML model training and APA event scoring.
include { APA_MODEL_TRAIN   } from '../../../modules/local/ml/apa_model_train/main'
include { APA_MODEL_PREDICT } from '../../../modules/local/ml/apa_model_predict/main'

workflow MODEL_PIPELINE {
    take:
    feature_table                  // path: apa_features.tsv
    apa_events                     // path: apa_events.tsv
    model_type                     // val: 'random_forest'|'xgboost'|'logistic_regression'
    group_level                    // val: 'cluster'|'cell_type'
    enable_single_cell_projection  // val: boolean

    main:
    APA_MODEL_TRAIN(
        feature_table,
        apa_events,
        model_type,
        group_level
    )

    APA_MODEL_PREDICT(
        feature_table,
        APA_MODEL_TRAIN.out.trained_model,
        group_level,
        enable_single_cell_projection
    )

    emit:
    trained_model      = APA_MODEL_TRAIN.out.trained_model      // path: apa_model.pkl
    feature_importance = APA_MODEL_TRAIN.out.feature_importance  // path: feature_importance.tsv
    model_metrics      = APA_MODEL_TRAIN.out.metrics             // path: model_metrics.tsv
    scored_events      = APA_MODEL_PREDICT.out.scored_events     // path: scored_apa_events.tsv
}
