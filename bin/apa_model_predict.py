#!/usr/bin/env python3
"""
apa_model_predict.py — Score candidate APA events with a trained model bundle.

The model bundle (apa_model.pkl) must be a dict produced by apa_model_train.py:
    {"model": <sklearn estimator>, "feature_cols": [str, ...]}

Outputs scored_apa_events.tsv with columns:
    site_id, group_id, group_level, apa_score, is_apa

Options
-------
  --group-level                 Filter output to this group level (default: cluster)
  --single-cell-projection      Attach per-barcode scores from cell_annotations
  --threshold                   Probability threshold for is_apa (default: 0.5)
"""
import argparse
import logging
import pickle

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--features",                required=True)
    p.add_argument("--model",                   required=True)
    p.add_argument("--group-level",             default="cluster")
    p.add_argument("--single-cell-projection",  action="store_true",
                   help="Add per-barcode APA scores via cell_annotations join")
    p.add_argument("--threshold",               type=float, default=0.5)
    p.add_argument("--out-scored",              required=True)
    p.add_argument("--log",                     required=True)
    return p.parse_args()


def _empty_output(path: str, log):
    pd.DataFrame(
        columns=["site_id", "group_id", "group_level", "apa_score", "is_apa"]
    ).to_csv(path, sep="\t", index=False)
    log.warning(f"Wrote empty scored output → {path}")


def main():
    args = parse_args()
    logging.basicConfig(
        filename=args.log, level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger()
    log.info("Starting APA model prediction")

    # Load feature table
    try:
        features = pd.read_csv(args.features, sep="\t")
    except Exception as e:
        log.error(f"Cannot load features: {e}")
        _empty_output(args.out_scored, log)
        return

    # Load model bundle
    try:
        with open(args.model, "rb") as f:
            bundle = pickle.load(f)
        model        = bundle.get("model")
        feature_cols = bundle.get("feature_cols", [])
    except Exception as e:
        log.error(f"Cannot load model: {e}")
        _empty_output(args.out_scored, log)
        return

    if model is None:
        log.warning("Model is None (empty bundle); assigning constant score 0.5")
        features["apa_score"] = 0.5
        features["is_apa"]    = False
    else:
        available = [c for c in feature_cols if c in features.columns]
        if not available:
            log.error("No matching feature columns found in feature table")
            _empty_output(args.out_scored, log)
            return

        X = features[available].fillna(0).values

        if hasattr(model, "predict_proba"):
            proba = model.predict_proba(X)
            # When only one class was seen during training, proba has shape (n, 1)
            if proba.shape[1] > 1:
                scores = proba[:, 1]
            else:
                log.warning("Model trained on single class; assigning constant score 0.5")
                scores = np.full(len(X), 0.5)
        elif hasattr(model, "decision_function"):
            raw = model.decision_function(X)
            # Min-max normalize to [0, 1]
            rng = raw.max() - raw.min()
            scores = (raw - raw.min()) / (rng + 1e-9)
        else:
            scores = model.predict(X).astype(float)

        features["apa_score"] = np.round(scores, 4)
        features["is_apa"]    = scores >= args.threshold

    # Filter to requested group_level
    if "group_level" in features.columns:
        out = features[features["group_level"] == args.group_level].copy()
    else:
        out = features.copy()
        out["group_level"] = args.group_level

    n_apa = int(out["is_apa"].sum())
    log.info(
        f"Scored {len(out)} site-group pairs in group_level='{args.group_level}'; "
        f"{n_apa} predicted APA (threshold={args.threshold})"
    )

    keep_cols = ["site_id", "group_id", "group_level", "apa_score", "is_apa"]
    missing = [c for c in keep_cols if c not in out.columns]
    for c in missing:
        out[c] = "NA"

    out[keep_cols].to_csv(args.out_scored, sep="\t", index=False)
    log.info(f"Scored events written → {args.out_scored}")


if __name__ == "__main__":
    main()
