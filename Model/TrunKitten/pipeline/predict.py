"""Predict NMD escape with TrunKitten — the reduced top-10 feature CatBoost
model derived from TrunCat.

Consumes the annotated.tsv produced by the TrunKitten annotation CLI
(`python -m minicat.cli`) and writes a prediction table with columns:
variant_id, escape_prob, escape_pred_at_youden.

Usage:
    python predict.py \
        --annotated outputs/annotated.tsv \
        --model     models/trunkitten/trunkitten.pkl \
        --metadata  results/trunkitten/top10_features.json \
        --out       outputs/predictions.tsv \
        --training-medians results/trunkitten/training_medians.json  # optional

File paths shown above reflect the TrunKitten naming. If your saved
artifacts still use the historical `reduced_top10` naming, pass those
paths — the script only cares about the file contents, not the names.

The training-medians JSON lets this script apply the SAME imputation /
zero-fill rules that the TrunCat training pipeline applied (Notebook 02).
If omitted, missing values are left as NaN and CatBoost's native NaN
handling is used. Either is defensible; using training medians reproduces
training behaviour exactly. Recommended for external-cohort scoring.

Expected training_medians.json schema:
{
  "median_impute": {"half_life_PC1": 1.1303, ...},
  "zero_fill":     ["cdsseq_AUcontentlast200",
                    "phastcons_new3utr_first200_median",
                    "phylop_ptc_to_ejc_median"]
}
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from catboost import Pool


def apply_training_imputation(
    X: pd.DataFrame, medians_path: Path | None,
) -> tuple[pd.DataFrame, dict]:
    """Apply training-time median imputation / zero-fill. Return (X_imputed, report)."""
    report = {"imputed": {}, "zero_filled": {}}
    if medians_path is None:
        return X, report

    spec = json.loads(Path(medians_path).read_text())

    for col, med in spec.get("median_impute", {}).items():
        if col in X.columns:
            n_miss = X[col].isna().sum()
            if n_miss:
                X[col] = X[col].fillna(float(med))
                report["imputed"][col] = {"n": int(n_miss), "value": float(med)}

    for col in spec.get("zero_fill", []):
        if col in X.columns:
            n_miss = X[col].isna().sum()
            if n_miss:
                X[col] = X[col].fillna(0.0)
                report["zero_filled"][col] = int(n_miss)

    return X, report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="TrunKitten scorer — predicts NMD escape using the "
                    "reduced top-10 feature model derived from TrunCat."
    )
    ap.add_argument("--annotated", required=True,
                    help="TSV from the TrunKitten annotation CLI "
                         "(python -m minicat.cli ...)")
    ap.add_argument("--model",     required=True,
                    help=".pkl of the TrunKitten CatBoost model")
    ap.add_argument("--metadata",  required=True,
                    help="top10_features.json (feature order, dtypes, threshold)")
    ap.add_argument("--out",       required=True,
                    help="output predictions TSV")
    ap.add_argument("--training-medians", default=None,
                    help="optional JSON with TrunCat training-time "
                         "median/zero-fill rules")
    args = ap.parse_args(argv)

    # --- Load feature metadata & model ---
    meta = json.loads(Path(args.metadata).read_text())
    required_features = meta["features_in_order"]
    cat_features      = meta["categorical_features"]
    threshold         = float(meta["youden_threshold"])
    model             = joblib.load(args.model)

    # --- Load annotated table ---
    df = pd.read_csv(args.annotated, sep="\t")

    missing_cols = [c for c in required_features if c not in df.columns]
    if missing_cols:
        raise KeyError(
            f"Annotated TSV is missing required feature columns: {missing_cols}. "
            f"Re-run the TrunKitten annotation CLI (python -m minicat.cli)."
        )

    # --- Feature matrix in training column order ---
    X = df[required_features].copy()

    # Apply training imputation if provided
    X, impute_report = apply_training_imputation(
        X, Path(args.training_medians) if args.training_medians else None,
    )

    # Categorical columns must be string with "NA" for missing — matches training
    for c in cat_features:
        X[c] = X[c].astype(str)
        # Pandas converts NaN → 'nan'; standardise to 'NA' to match training
        X[c] = X[c].replace({"nan": "NA", "None": "NA"})

    # --- Predict ---
    cat_idx = [X.columns.get_loc(c) for c in cat_features]
    pool = Pool(X, cat_features=cat_idx)
    probs = model.predict_proba(pool)[:, 1]

    # --- Write ---
    out = df[["variant_id"]].copy()
    if "gene" in df.columns:
        out["gene"] = df["gene"]
    if "txnames" in df.columns:
        out["txnames"] = df["txnames"]
    out["escape_prob"]             = probs
    out["escape_pred_at_youden"]   = (probs > threshold).astype(int)
    out["threshold_used"]          = threshold

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    # --- Summary ---
    n = len(out)
    n_escape = int(out["escape_pred_at_youden"].sum())
    print(f"[TrunKitten] scored {n:,} variants → {out_path}")
    print(f"             predicted escape: {n_escape:,} ({n_escape/max(n,1)*100:.1f}%)")
    print(f"             threshold: {threshold:.4f}")
    if impute_report["imputed"]:
        print(f"             median-imputed: {impute_report['imputed']}")
    if impute_report["zero_filled"]:
        print(f"             zero-filled: {impute_report['zero_filled']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
