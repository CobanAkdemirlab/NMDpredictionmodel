#!/usr/bin/env python
"""
predict.py — spooky_model v4.1 inference on new variant cohorts.

Takes a variant CSV that already has the baseline enhanced features plus the
six v4 annotation sources (codon optimality, readthrough, EJC occupancy,
PTC AUG, conservation, isoform counts — or pointers to them), and produces
per-variant NMD-escape predictions using the final CatBoost model.

Design principles
-----------------
1. The trained model's `feature_names_` is the SOURCE OF TRUTH for which
   features to use and in what order. We never re-derive the feature list.
2. We only apply DETERMINISTIC cleaning from notebook 02 (zero-fill UTR
   structural absence, gnomAD zero-fill, categorical NaN -> 'MISSING',
   loeuf_cat rename). Data-dependent imputations (median fills) use
   precomputed training medians if available; otherwise CatBoost's native
   NaN handling kicks in.
3. Verification-first: we audit every column against the model's expectations
   BEFORE predicting, and fail loudly on mismatches.

Usage
-----
    # Cohort CSV that already has all v4 features merged in:
    python predict.py \
        --input  /Users/jschmidt3/Iman_visualizations/spooky_model_v4.1/predict/gnomAD_df_stopgain_updated2026_merged.csv \
        --output /Users/jschmidt3/Iman_visualizations/spooky_model_v4.1/predict/gnomAD_predictions.csv \
        --label  gnomAD_2026

    # Cohort CSV that still needs v4 annotations merged in (point at the
    # annotation files via --codon-opt, --readthrough, etc.):
    python predict.py \
        --input  .../predict/clinvar_df_stopgain_updated2026.csv \
        --output .../predict/clinvar_predictions.csv \
        --codon-opt    .../predict/annotations/clinvar_codon_optimality.tsv \
        --readthrough  .../predict/annotations/clinvar_readthrough.csv \
        --ejc          .../predict/annotations/clinvar_ejc_occupancy.tsv \
        --ptc-aug      .../predict/annotations/clinvar_ptc_aug.tsv \
        --conservation .../predict/annotations/clinvar_conservation.csv \
        --label        clinvar_2026
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from catboost import CatBoostClassifier, Pool
except ImportError as e:  # pragma: no cover
    print("ERROR: catboost is required. pip install catboost", file=sys.stderr)
    raise

# ──────────────────────────────────────────────────────────────────────────────
# Constants mirroring notebook 02's deterministic cleaning rules
# ──────────────────────────────────────────────────────────────────────────────

# Structurally absent (zero-fill is correct — region does not exist)
UTR_ZERO_FILL = [
    "threeUTR_length", "threeUTR_AU_content", "threeUTR_UC_content",
    "ThreeUTR_AUcontentlast200", "ThreeUTR_AUcontentfirst200",
    "ThreeUTR_UCcontentlast200", "ThreeUTR_UCcontentfirst200",
    "ThreeUTR_AUcontentlast100", "ThreeUTR_AUcontentfirst100",
    "ThreeUTR_UCcontentlast100", "ThreeUTR_UCcontentfirst100",
    "fiveutr_length", "fiveUTR_AU_content", "fiveUTR_UC_content",
    "fiveUTR_AUcontentlast200", "fiveUTR_AUcontentfirst200",
    "fiveUTR_UCcontentlast200", "fiveUTR_UCcontentfirst200",
    "fiveUTR_AUcontentlast100", "fiveUTR_AUcontentfirst100",
    "fiveUTR_UCcontentlast100", "fiveUTR_UCcontentfirst100",
    "newUTR_length", "log2newUTR", "log2_3utr", "log2_5utr",
    "cdsseq_AUcontentlast200", "cdsseq_UCcontentlast200",
    "cdsseq_AUcontentfirst200", "cdsseq_UCcontentfirst200",
    "phastcons_utr5_first200_median", "phastcons_utr5_whole_median",
    "phylop_utr5_first200_median",    "phylop_utr5_whole_median",
    "phastcons_old3utr_first200_median", "phastcons_old3utr_whole_median",
    "phylop_old3utr_first200_median",    "phylop_old3utr_whole_median",
    "phastcons_ptc_to_ejc_median",
    "phylop_ptc_to_ejc_median",
]

# Median-imputed at training time (we load saved training medians for these)
MEDIAN_IMPUTE_EXPLICIT = [
    "pLI", "oe_lof_upper",
    "MedianExpression", "MedianExpression_log2", "Whole.Blood",
    "half_life_PC1",
    "CADD_phred",
    "readthrough_score_hek293t",
]


# ──────────────────────────────────────────────────────────────────────────────
# Notebook 01 merge logic — deterministic, safe to re-run on any cohort
# ──────────────────────────────────────────────────────────────────────────────

def _categorize_aug_distance(dist: float) -> str:
    if pd.isna(dist):
        return "no_inframe_AUG"
    if dist < 50:
        return "very_close"
    if dist < 100:
        return "close"
    if dist < 200:
        return "moderate"
    if dist < 500:
        return "far"
    return "very_far"


def _engineer_aug_features(df_ptcaug: pd.DataFrame) -> pd.DataFrame:
    """Mirror of notebook 01 section 7 — deterministic transforms."""
    df = df_ptcaug.copy()
    df["aug_distance_nt"] = df["nearest_inframe_aug_distance_nt"]
    df["aug_distance_category"] = df["aug_distance_nt"].apply(_categorize_aug_distance)
    df["kozak_strength"] = df["nearest_inframe_kozak_strength"].fillna("no_inframe_AUG").astype(str)
    df["has_plus1_aug"] = df["has_plus1_frame_aug"].fillna(False).astype(str)
    df["has_plus2_aug"] = df["has_plus2_frame_aug"].fillna(False).astype(str)

    def _frame_status(row):
        p1 = str(row["has_plus1_aug"]) == "True"
        p2 = str(row["has_plus2_aug"]) == "True"
        if p1 and p2:
            return "both_frames"
        if p1:
            return "plus1_only"
        if p2:
            return "plus2_only"
        return "no_frame_AUG"

    df["aug_frame_status"] = df.apply(_frame_status, axis=1)
    return df


def merge_annotations(
    df_enhanced: pd.DataFrame,
    *,
    path_codon_opt: Path | None,
    path_readthrough: Path | None,
    path_ejc: Path | None,
    path_ptc_aug: Path | None,
    path_conservation: Path | None,
) -> pd.DataFrame:
    """Mirror of notebook 01 — join the six v4 annotation sources onto
    the enhanced baseline CSV by variantID. Only merges sources that are
    actually provided; silently skips any that are None (assumes they're
    already in df_enhanced).
    """
    if "variantID" not in df_enhanced.columns:
        raise ValueError("Input CSV must contain a 'variantID' column.")

    # Drop the features notebook 01 removes before merging
    for old_col in ("AverageCodonRNAUsage", "median_half_life"):
        if old_col in df_enhanced.columns:
            df_enhanced = df_enhanced.drop(columns=[old_col])

    df_result = df_enhanced.copy()
    main_cols = set(df_enhanced.columns)

    # ── Codon optimality ─────────────────────────────────────────────────────
    if path_codon_opt is not None:
        df_codon = pd.read_csv(path_codon_opt, sep="\t")
        if "PTC_ID" in df_codon.columns and "variantID" not in df_codon.columns:
            df_codon["variantID"] = df_codon["PTC_ID"]
        codon_cols = ["variantID", "CodonOptimalityFraction_CDS", "CodonOptimalityFraction_PTCpm100nt"]
        codon_cols = [c for c in codon_cols if c in df_codon.columns]
        df_result = df_result.merge(df_codon[codon_cols], on="variantID", how="left")
        print(f"  ✓ Merged codon optimality: {df_result.shape}")

    # ── Readthrough ──────────────────────────────────────────────────────────
    if path_readthrough is not None:
        df_rt = pd.read_csv(path_readthrough)
        rt_cols = ["variantID", "readthrough_score_hek293t", "readthrough_category_hek293t"]
        rt_cols = [c for c in rt_cols if c in df_rt.columns]
        df_result = df_result.merge(df_rt[rt_cols], on="variantID", how="left")
        print(f"  ✓ Merged readthrough: {df_result.shape}")

    # ── EJC occupancy ────────────────────────────────────────────────────────
    if path_ejc is not None:
        df_ejc = pd.read_csv(path_ejc, sep="\t")
        if "PTC_ID" in df_ejc.columns:
            df_ejc = df_ejc.rename(columns={"PTC_ID": "variantID"})
        ejc_new = ["variantID"] + [c for c in df_ejc.columns
                                    if c not in main_cols and c != "variantID"]
        df_result = df_result.merge(df_ejc[ejc_new], on="variantID", how="left")
        print(f"  ✓ Merged EJC occupancy: {df_result.shape}")

    # ── PTC AUG ──────────────────────────────────────────────────────────────
    if path_ptc_aug is not None:
        df_aug = pd.read_csv(path_ptc_aug, sep="\t")
        if "variant_id" in df_aug.columns:
            df_aug = df_aug.rename(columns={"variant_id": "variantID"})
        df_aug = _engineer_aug_features(df_aug)
        aug_keep = [
            "variantID",
            "aug_distance_category", "aug_distance_nt",
            "kozak_strength", "has_plus1_aug", "has_plus2_aug", "aug_frame_status",
        ]
        aug_keep = [c for c in aug_keep if c in df_aug.columns]
        df_result = df_result.merge(df_aug[aug_keep], on="variantID", how="left")
        print(f"  ✓ Merged PTC AUG (with engineered features): {df_result.shape}")

    # ── Conservation ─────────────────────────────────────────────────────────
    if path_conservation is not None:
        df_cons = pd.read_csv(path_conservation)
        cons_new = ["variantID"] + [c for c in df_cons.columns
                                    if c.startswith(("phastcons_", "phylop_"))]
        df_result = df_result.merge(df_cons[cons_new], on="variantID", how="left")
        print(f"  ✓ Merged conservation: {df_result.shape}")

    return df_result


# ──────────────────────────────────────────────────────────────────────────────
# Feature alignment and deterministic cleaning
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class AuditReport:
    """Everything we want to know before we let CatBoost touch the data."""
    model_features: list[str]
    model_cat_features: list[str]

    present: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)   # expected by model, absent from input
    extra: list[str] = field(default_factory=list)     # in input, unused by model
    dtype_coerced: list[tuple[str, str, str]] = field(default_factory=list)

    zero_filled_structural: list[tuple[str, int]] = field(default_factory=list)
    zero_filled_gnomad: int = 0
    median_imputed: list[tuple[str, int, float]] = field(default_factory=list)
    categorical_missing_filled: list[tuple[str, int]] = field(default_factory=list)

    residual_nans: dict[str, int] = field(default_factory=dict)

    def print(self) -> None:
        print("\n" + "=" * 80)
        print("FEATURE AUDIT REPORT")
        print("=" * 80)
        print(f"\nModel expects {len(self.model_features)} features "
              f"({len(self.model_cat_features)} categorical)")
        print(f"  Present in input: {len(self.present)}")
        print(f"  Missing from input: {len(self.missing)}")
        print(f"  Extra in input (will be dropped): {len(self.extra)}")

        if self.missing:
            print(f"\n  ⚠️  Missing features (filled with NaN — CatBoost will handle natively):")
            for feat in self.missing[:20]:
                is_cat = " [CAT]" if feat in self.model_cat_features else ""
                print(f"    - {feat}{is_cat}")
            if len(self.missing) > 20:
                print(f"    ... and {len(self.missing) - 20} more")

        if self.dtype_coerced:
            print(f"\n  Dtype coercions applied: {len(self.dtype_coerced)}")
            for feat, was, now in self.dtype_coerced[:10]:
                print(f"    - {feat}: {was} → {now}")
            if len(self.dtype_coerced) > 10:
                print(f"    ... and {len(self.dtype_coerced) - 10} more")

        print(f"\nDeterministic imputation (from notebook 02):")
        print(f"  Structural zero-fills: {len(self.zero_filled_structural)} features")
        for feat, n in self.zero_filled_structural[:5]:
            print(f"    - {feat}: {n} NaNs → 0")
        if len(self.zero_filled_structural) > 5:
            print(f"    ... and {len(self.zero_filled_structural) - 5} more")
        print(f"  gnomAD_exome_ALL zero-fills: {self.zero_filled_gnomad}")
        print(f"  Categorical MISSING fills: {len(self.categorical_missing_filled)} features")

        if self.median_imputed:
            print(f"\nMedian imputation (from training medians):")
            for feat, n, med in self.median_imputed:
                print(f"    - {feat}: {n} NaNs → {med:.4f}")

        if self.residual_nans:
            print(f"\nResidual NaNs (will be handled by CatBoost natively):")
            for feat, n in sorted(self.residual_nans.items(), key=lambda x: -x[1])[:10]:
                is_cat = " [CAT]" if feat in self.model_cat_features else ""
                print(f"    - {feat}{is_cat}: {n}")
            if len(self.residual_nans) > 10:
                print(f"    ... and {len(self.residual_nans) - 10} more")
        else:
            print(f"\n✓ No residual NaNs")

    def to_dict(self) -> dict[str, Any]:
        return {
            "model_features_expected": len(self.model_features),
            "model_cat_features": len(self.model_cat_features),
            "present": len(self.present),
            "missing": self.missing,
            "extra": self.extra,
            "dtype_coerced": [list(t) for t in self.dtype_coerced],
            "zero_filled_structural": [list(t) for t in self.zero_filled_structural],
            "zero_filled_gnomad": self.zero_filled_gnomad,
            "median_imputed": [list(t) for t in self.median_imputed],
            "categorical_missing_filled": [list(t) for t in self.categorical_missing_filled],
            "residual_nans": self.residual_nans,
        }


def _get_model_feature_info(model: CatBoostClassifier) -> tuple[list[str], list[str]]:
    """Return (all_feature_names, categorical_feature_names) from the model."""
    feat_names = list(model.feature_names_)
    cat_idx = model.get_cat_feature_indices() if hasattr(model, "get_cat_feature_indices") else []
    cat_names = [feat_names[i] for i in cat_idx]
    return feat_names, cat_names


def clean_and_align(
    df: pd.DataFrame,
    model: CatBoostClassifier,
    training_medians: dict[str, float] | None,
) -> tuple[pd.DataFrame, AuditReport]:
    """
    Apply notebook 02's deterministic cleaning steps to the merged cohort CSV
    and align columns to the model's expected feature list/order.

    Returns (X_aligned, audit) where X_aligned is ready to go straight into
    a CatBoost Pool.
    """
    model_features, model_cat_features = _get_model_feature_info(model)
    audit = AuditReport(model_features=model_features, model_cat_features=model_cat_features)

    X = df.copy()

    # ── loeuf_cat naming normalisation (nb 02) ───────────────────────────────
    for variant in ("loefu_cat", "LOEUF_cat"):
        if variant in X.columns and "loeuf_cat" not in X.columns:
            X = X.rename(columns={variant: "loeuf_cat"})

    # ── Structural zero-fills (nb 02) ────────────────────────────────────────
    for col in UTR_ZERO_FILL:
        if col in X.columns:
            n = X[col].isna().sum()
            if n > 0:
                X[col] = X[col].fillna(0)
                audit.zero_filled_structural.append((col, int(n)))

    # ── gnomAD zero-fill (nb 02) ─────────────────────────────────────────────
    if "gnomAD_exome_ALL" in X.columns:
        n = X["gnomAD_exome_ALL"].isna().sum()
        if n > 0:
            X["gnomAD_exome_ALL"] = X["gnomAD_exome_ALL"].fillna(0)
            audit.zero_filled_gnomad = int(n)

    # ── Median imputation from TRAINING medians (nb 02) ──────────────────────
    if training_medians:
        for col in MEDIAN_IMPUTE_EXPLICIT:
            if col in X.columns and col in training_medians:
                n = X[col].isna().sum()
                if n > 0:
                    med = float(training_medians[col])
                    X[col] = X[col].fillna(med)
                    audit.median_imputed.append((col, int(n), med))
    else:
        warnings.warn(
            "No training_medians provided — median-impute features will retain NaN. "
            "CatBoost will handle NaNs natively, but this is NOT identical to training. "
            "Run export_training_medians.py against TOPMed_cleaned_v4.csv to fix.",
            stacklevel=2,
        )

    # ── Categorical handling (nb 02): cast to str, NaN → 'MISSING' ───────────
    for col in model_cat_features:
        if col in X.columns:
            n_missing = X[col].isna().sum()
            X[col] = X[col].astype(str)
            X[col] = X[col].replace({"nan": "MISSING", "None": "MISSING"})
            if n_missing > 0:
                audit.categorical_missing_filled.append((col, int(n_missing)))

    # ── Align to model feature list ──────────────────────────────────────────
    present = [c for c in model_features if c in X.columns]
    missing = [c for c in model_features if c not in X.columns]
    extra = [c for c in X.columns if c not in model_features]

    audit.present = present
    audit.missing = missing
    audit.extra = extra

    # Add missing columns as NaN (CatBoost handles them; categorical missings
    # become 'MISSING' string so they don't crash Pool construction)
    for col in missing:
        if col in model_cat_features:
            X[col] = "MISSING"
        else:
            X[col] = np.nan

    # Select and reorder to match the model exactly
    X_aligned = X[model_features].copy()

    # ── Final dtype coercion ─────────────────────────────────────────────────
    # CatBoost requires cat features to be string-like, numeric features numeric.
    # Note: pandas 2.x uses a distinct `str` dtype (StringDtype) separate from
    # `object`; both are acceptable to CatBoost, so we only coerce if neither.
    for col in model_features:
        if col in model_cat_features:
            is_string_like = (
                X_aligned[col].dtype == object
                or pd.api.types.is_string_dtype(X_aligned[col])
            )
            if not is_string_like:
                old = str(X_aligned[col].dtype)
                X_aligned[col] = X_aligned[col].astype(str)
                audit.dtype_coerced.append((col, old, "str"))
        else:
            if not pd.api.types.is_numeric_dtype(X_aligned[col]):
                old = str(X_aligned[col].dtype)
                X_aligned[col] = pd.to_numeric(X_aligned[col], errors="coerce")
                audit.dtype_coerced.append((col, old, "float64"))

    # ── Residual NaN tally ───────────────────────────────────────────────────
    for col in model_features:
        if col in model_cat_features:
            continue  # already handled with 'MISSING'
        n = X_aligned[col].isna().sum()
        if n > 0:
            audit.residual_nans[col] = int(n)

    return X_aligned, audit


# ──────────────────────────────────────────────────────────────────────────────
# Prediction + output
# ──────────────────────────────────────────────────────────────────────────────

def predict_cohort(
    X_aligned: pd.DataFrame,
    model: CatBoostClassifier,
    variant_ids: pd.Series,
    threshold: float,
    extra_id_cols: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Run predict_proba and assemble the output frame with threshold call."""
    _, model_cat_features = _get_model_feature_info(model)
    cat_indices = [X_aligned.columns.get_loc(c) for c in model_cat_features if c in X_aligned.columns]

    pool = Pool(X_aligned, cat_features=cat_indices)
    probs = model.predict_proba(pool)[:, 1]

    out_cols = {}
    if extra_id_cols is not None:
        for c in extra_id_cols.columns:
            out_cols[c] = extra_id_cols[c].values
    out_cols["variantID"] = variant_ids.values
    out_cols["escape_probability"] = probs
    out_cols["predicted_class"] = (probs >= threshold).astype(int)
    out_cols["predicted_label"] = np.where(probs >= threshold, "escape", "NMD")
    out_cols["threshold_used"] = threshold

    return pd.DataFrame(out_cols)


def sanity_check_predictions(pred_df: pd.DataFrame, threshold: float) -> None:
    print("\n" + "=" * 80)
    print("PREDICTION SANITY CHECK")
    print("=" * 80)
    n = len(pred_df)
    n_esc = int((pred_df["predicted_class"] == 1).sum())
    n_nmd = n - n_esc
    p = pred_df["escape_probability"]
    print(f"\nTotal variants: {n}")
    print(f"Threshold applied: {threshold:.4f}")
    print(f"\nClass distribution:")
    print(f"  escape: {n_esc:,} ({n_esc/n*100:.1f}%)")
    print(f"  NMD:    {n_nmd:,} ({n_nmd/n*100:.1f}%)")
    print(f"\nEscape probability distribution:")
    print(f"  min:  {p.min():.4f}")
    print(f"  p25:  {p.quantile(0.25):.4f}")
    print(f"  med:  {p.median():.4f}")
    print(f"  p75:  {p.quantile(0.75):.4f}")
    print(f"  max:  {p.max():.4f}")
    print(f"  mean: {p.mean():.4f}")
    # variantID uniqueness
    if pred_df["variantID"].duplicated().any():
        n_dup = int(pred_df["variantID"].duplicated().sum())
        print(f"\n⚠️  {n_dup} duplicated variantIDs in output")
    else:
        print(f"\n✓ All variantIDs unique")
    # NaN check
    if p.isna().any():
        print(f"⚠️  {p.isna().sum()} NaN probabilities — something went wrong")
    else:
        print(f"✓ No NaN probabilities")


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="spooky_model v4.1 inference on new variant cohorts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True, type=Path,
                   help="Variant CSV with baseline enhanced features (and optionally the v4 "
                        "annotations already merged in). Must contain variantID column.")
    p.add_argument("--output", required=True, type=Path,
                   help="Path to write predictions CSV.")
    p.add_argument("--model", type=Path,
                   default=Path("models/catboost_final_model.cbm"),
                   help="Path to trained CatBoost model (.cbm).")
    p.add_argument("--training-medians", type=Path,
                   default=Path("models/training_medians.json"),
                   help="JSON of medians from training cleaned CSV (see export_training_medians.py).")
    p.add_argument("--viz-summary", type=Path,
                   default=Path("results/figures/viz_summary.json"),
                   help="viz_summary.json from notebook 03, used for the Youden threshold.")
    p.add_argument("--threshold", type=float, default=None,
                   help="Override threshold instead of reading from viz_summary.json.")
    p.add_argument("--label", default=None,
                   help="Optional cohort label, used in logs and audit JSON filename.")

    # Optional annotation sources — if provided, we'll merge them in (nb 01 logic)
    p.add_argument("--codon-opt", type=Path, default=None)
    p.add_argument("--readthrough", type=Path, default=None)
    p.add_argument("--ejc", type=Path, default=None)
    p.add_argument("--ptc-aug", type=Path, default=None)
    p.add_argument("--conservation", type=Path, default=None)

    p.add_argument("--extra-id-cols", nargs="*", default=None,
                   help="Additional columns from the input to carry through to the output "
                        "(e.g. gene_symbol, clinvar_id, consequence).")
    p.add_argument("--audit-json", type=Path, default=None,
                   help="Optional path to dump the audit report as JSON.")
    return p.parse_args()


def main() -> int:
    args = _parse_args()

    label = args.label or args.input.stem
    print("=" * 80)
    print(f"spooky_model v4.1 inference: {label}")
    print("=" * 80)

    # ── Load model ───────────────────────────────────────────────────────────
    if not args.model.exists():
        print(f"ERROR: model file not found: {args.model}", file=sys.stderr)
        return 1
    model = CatBoostClassifier()
    model.load_model(str(args.model))
    feats, cat_feats = _get_model_feature_info(model)
    print(f"\n✓ Loaded model: {args.model}")
    print(f"  Features expected: {len(feats)}")
    print(f"  Categorical features: {len(cat_feats)}")

    # ── Load threshold ───────────────────────────────────────────────────────
    if args.threshold is not None:
        threshold = float(args.threshold)
        print(f"✓ Threshold (override): {threshold:.4f}")
    else:
        if not args.viz_summary.exists():
            print(f"ERROR: --viz-summary not found and --threshold not given: {args.viz_summary}",
                  file=sys.stderr)
            return 1
        with open(args.viz_summary) as f:
            summary = json.load(f)
        threshold = float(summary["threshold_used"])
        print(f"✓ Threshold from {args.viz_summary.name}: {threshold:.4f}")

    # ── Load training medians (optional but recommended) ────────────────────
    training_medians: dict[str, float] | None = None
    if args.training_medians.exists():
        with open(args.training_medians) as f:
            training_medians = json.load(f)
        print(f"✓ Loaded training medians: {len(training_medians)} features")
    else:
        print(f"⚠️  Training medians not found: {args.training_medians}")
        print(f"    CatBoost native NaN handling will be used as fallback.")

    # ── Load + merge input ──────────────────────────────────────────────────
    print(f"\nLoading input: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  Shape: {df.shape}")
    if "variantID" not in df.columns:
        print("ERROR: input CSV must contain a 'variantID' column.", file=sys.stderr)
        return 1

    any_annot = any([args.codon_opt, args.readthrough, args.ejc, args.ptc_aug, args.conservation])
    if any_annot:
        print("\nMerging v4 annotation sources...")
        df = merge_annotations(
            df,
            path_codon_opt=args.codon_opt,
            path_readthrough=args.readthrough,
            path_ejc=args.ejc,
            path_ptc_aug=args.ptc_aug,
            path_conservation=args.conservation,
        )
    else:
        print("\nNo annotation paths provided — assuming all v4 features already in --input.")

    # Stash variantID (and any extra id cols) BEFORE cleaning drops them
    variant_ids = df["variantID"].copy()
    extra_ids_df: pd.DataFrame | None = None
    if args.extra_id_cols:
        missing_extra = [c for c in args.extra_id_cols if c not in df.columns]
        if missing_extra:
            print(f"⚠️  --extra-id-cols not found in input (will skip): {missing_extra}")
        extra_ids_df = df[[c for c in args.extra_id_cols if c in df.columns]].copy()

    # ── Clean + align ────────────────────────────────────────────────────────
    X_aligned, audit = clean_and_align(df, model, training_medians)
    audit.print()

    if args.audit_json is not None:
        args.audit_json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.audit_json, "w") as f:
            json.dump(audit.to_dict(), f, indent=2)
        print(f"\n✓ Audit report written: {args.audit_json}")

    # ── Predict ──────────────────────────────────────────────────────────────
    print(f"\nRunning predictions on {len(X_aligned)} variants...")
    pred_df = predict_cohort(
        X_aligned, model, variant_ids,
        threshold=threshold, extra_id_cols=extra_ids_df,
    )

    sanity_check_predictions(pred_df, threshold)

    # ── Save ─────────────────────────────────────────────────────────────────
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pred_df.to_csv(args.output, index=False)
    print(f"\n✓ Predictions written: {args.output}")
    print(f"  Rows: {len(pred_df)}, Cols: {pred_df.shape[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
