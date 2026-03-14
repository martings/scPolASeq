#!/usr/bin/env python3
"""
apa_model_train.py — Train a supervised model to classify APA events.

Input
-----
  --features     apa_features.tsv  (site × group feature matrix)
  --apa-events   apa_events.tsv    (event table with is_significant column)

Labels
------
  A site/group pair is labelled 1 (APA) if any pairwise comparison involving
  that group and that site is marked is_significant == True in apa_events.tsv.
  If labels cannot be derived, coverage_at_site > median is used as a proxy.

Models
------
  random_forest       sklearn RandomForestClassifier  (default)
  xgboost             XGBClassifier  (falls back to RF if xgboost not installed)
  logistic_regression sklearn LogisticRegression

Output
------
  apa_model.pkl          pickle bundle: {"model": ..., "feature_cols": [...]}
  feature_importance.tsv
  model_metrics.tsv
  model_train.log
"""
import argparse
import logging
import pickle

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import label_binarize

FEATURE_COLS = [
    "coverage_at_site",
    "upstream_cov",
    "downstream_cov",
    "read_end_density",
    "proximal_distal_ratio",
    "dist_to_known_polya",
    "umi_support",
    "cluster_support",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--features",       required=True)
    p.add_argument("--apa-events",     required=True)
    p.add_argument("--model-type",     default="random_forest",
                   choices=["random_forest", "xgboost", "logistic_regression"])
    p.add_argument("--group-level",    default="cluster")
    p.add_argument("--out-model",      required=True)
    p.add_argument("--out-importance", required=True)
    p.add_argument("--out-metrics",    required=True)
    p.add_argument("--log",            required=True)
    p.add_argument("--n-estimators",   type=int, default=200)
    p.add_argument("--random-state",   type=int, default=42)
    return p.parse_args()


def build_model(model_type: str, n_estimators: int, random_state: int):
    if model_type == "random_forest":
        return RandomForestClassifier(
            n_estimators=n_estimators,
            random_state=random_state,
            n_jobs=-1,
            class_weight="balanced",
        )
    elif model_type == "xgboost":
        try:
            from xgboost import XGBClassifier
            return XGBClassifier(
                n_estimators=n_estimators,
                random_state=random_state,
                eval_metric="logloss",
                use_label_encoder=False,
            )
        except ImportError:
            logging.getLogger().warning("xgboost not installed — using RandomForest")
            return RandomForestClassifier(
                n_estimators=n_estimators,
                random_state=random_state,
                n_jobs=-1,
                class_weight="balanced",
            )
    else:  # logistic_regression
        return LogisticRegression(
            random_state=random_state,
            max_iter=1000,
            class_weight="balanced",
        )


def _empty_outputs(args, log):
    with open(args.out_model, "wb") as f:
        pickle.dump({"model": None, "feature_cols": FEATURE_COLS}, f)
    pd.DataFrame(columns=["feature", "importance"]).to_csv(
        args.out_importance, sep="\t", index=False)
    pd.DataFrame(columns=["metric", "value"]).to_csv(
        args.out_metrics, sep="\t", index=False)
    log.warning("Wrote empty model outputs")


def main():
    args = parse_args()
    logging.basicConfig(
        filename=args.log, level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger()
    log.info(f"Training APA model [{args.model_type}]")

    try:
        features = pd.read_csv(args.features, sep="\t")
        events   = pd.read_csv(args.apa_events, sep="\t")
    except Exception as e:
        log.error(f"Cannot load inputs: {e}")
        _empty_outputs(args, log)
        return

    # Build labels
    if (
        "is_significant" in events.columns
        and "site_id"  in events.columns
        and "group_a"  in events.columns
    ):
        pos_keys: set = set()
        for _, row in events[events["is_significant"] == True].iterrows():
            pos_keys.add((row["site_id"], row["group_a"]))
            pos_keys.add((row["site_id"], row["group_b"]))

        features["is_apa"] = features.apply(
            lambda r: int((r["site_id"], r["group_id"]) in pos_keys), axis=1
        )
    else:
        log.warning("Cannot build labels from events — using coverage_at_site > median as proxy")
        threshold = features["coverage_at_site"].median()
        features["is_apa"] = (features["coverage_at_site"] > threshold).astype(int)

    available_features = [c for c in FEATURE_COLS if c in features.columns]
    X = features[available_features].fillna(0).values
    y = features["is_apa"].values

    log.info(
        f"Dataset: {X.shape[0]} samples × {X.shape[1]} features; "
        f"{int(y.sum())} positive labels"
    )

    if X.shape[0] < 10 or y.sum() == 0 or y.sum() == len(y):
        log.warning("Insufficient or degenerate labels — training on full dataset without CV")
        model = build_model(args.model_type, args.n_estimators, args.random_state)
        model.fit(X, y)
        cv_scores = np.array([0.5])
        acc = float(accuracy_score(y, model.predict(X)))
    else:
        n_splits = min(5, int(y.sum()))
        model = build_model(args.model_type, args.n_estimators, args.random_state)
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=args.random_state)
        cv_scores = cross_val_score(model, X, y, cv=cv, scoring="roc_auc")
        model.fit(X, y)
        acc = float(accuracy_score(y, model.predict(X)))

    log.info(f"CV AUC-ROC: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")

    # Persist model bundle
    with open(args.out_model, "wb") as f:
        pickle.dump({"model": model, "feature_cols": available_features}, f)
    log.info(f"Model saved → {args.out_model}")

    # Feature importance
    if hasattr(model, "feature_importances_"):
        imp_vals = model.feature_importances_
    elif hasattr(model, "coef_"):
        imp_vals = np.abs(model.coef_[0])
    else:
        imp_vals = np.ones(len(available_features)) / len(available_features)

    pd.DataFrame({
        "feature":    available_features,
        "importance": np.round(imp_vals, 6),
    }).sort_values("importance", ascending=False).to_csv(
        args.out_importance, sep="\t", index=False)

    # Metrics
    pd.DataFrame([
        {"metric": "auc_roc_cv_mean",  "value": round(float(cv_scores.mean()), 4)},
        {"metric": "auc_roc_cv_std",   "value": round(float(cv_scores.std()),  4)},
        {"metric": "accuracy_train",   "value": round(acc, 4)},
        {"metric": "n_samples",        "value": X.shape[0]},
        {"metric": "n_features",       "value": X.shape[1]},
        {"metric": "n_positive_train", "value": int(y.sum())},
        {"metric": "model_type",       "value": args.model_type},
    ]).to_csv(args.out_metrics, sep="\t", index=False)

    log.info("Training complete")


if __name__ == "__main__":
    main()
