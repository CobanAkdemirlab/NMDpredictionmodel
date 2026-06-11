#!/usr/bin/env python3
"""
Apply the NMDetective-B decision tree (Lindeboom et al., Nat Genet 2019,
Fig. 1c) to a stop-gain variant table using locally-computed features.

This is a faithful reimplementation of the published 4-test tree, not a
lookup against the figshare resource. Doing it this way uses our chosen
transcript context (matching what TrunCat sees) and avoids the
GENCODE-vs-UCSC-knownGene annotation gap that makes coordinate lookup
against the figshare GTFs unreliable.

The tree (in evaluation order):

    i.   on_last_exon                         -> 0.00  ("last_exon")
    ii.  distance_to_start < 150 nt           -> 0.12  ("start_proximal")
    iii. exon_length > 407 nt                 -> 0.41  ("long_exon")
    iv.  in_last_50nt_of_penultimate_exon     -> 0.20  ("50nt_rule")
    else                                      -> 0.65  ("trigger_NMD")

Score values are the leaf NMD efficacy means reported in Fig. 1c of the
paper (NMD=0.00, 0.12, 0.41, 0.20, 0.65 respectively). Higher = more
decay = more NMD-sensitive.

The published cutoffs (paper text) are:
    score > 0.52  -> efficiently triggers NMD ("sensitive")
    score < 0.25  -> does not trigger NMD     ("escape")
    otherwise     -> intermediate

NMDetective-A (Random Forest) is not reimplemented here — it adds only
~3% R^2 over NMDetective-B (71% vs 68% on the held-out indel set per
Lindeboom 2019), and the paper itself uses NMDetective-B as the default
in downstream analyses.

Variant table column conventions (auto-detected; falls back across these
candidates so it works on TopMed_merged_v4.csv and friends):
- on_last_exon:                last.exon, on_last_exon, last_exon, isLastExon
- distance to coding start:    PTC.2.start, distance_to_start, dist_to_start
- exon length:                 current.exon.length, exon.length, exon_length
- last 50 nt of penultimate:   penultimate.last50bp, in_last_50nt_of_penultimate_exon
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants from Lindeboom et al. 2019, Fig. 1c
# ---------------------------------------------------------------------------
SCORE_LAST_EXON      = 0.00
SCORE_START_PROXIMAL = 0.12
SCORE_LONG_EXON      = 0.41
SCORE_50NT_RULE      = 0.20
SCORE_TRIGGER_NMD    = 0.65

# Paper cutoffs for class assignment (from the Results section)
THRESH_SENSITIVE = 0.52   # > this  -> efficient NMD
THRESH_ESCAPE    = 0.25   # < this  -> no NMD

# Tree thresholds
START_PROX_NT = 150
LONG_EXON_NT  = 407

# ---------------------------------------------------------------------------
# Column auto-detection (mirrors the rest of the repo's conventions)
# ---------------------------------------------------------------------------
CANDIDATES = {
    "last_exon": [
        "last.exon", "on_last_exon", "last_exon", "isLastExon", "in_last_exon",
    ],
    "dist_to_start": [
        "PTC.2.start", "distance_to_start", "dist_to_start",
        "distance_to_coding_start", "PTC_to_start",
    ],
    "exon_length": [
        "current.exon.length", "exon.length", "exon_length",
        "currentExonLength", "PTC_exon_length",
    ],
    "penult_last50": [
        "penultimate.last50bp", "in_last_50nt_of_penultimate_exon",
        "penult_last50", "last50_penultimate",
    ],
}


def _find_col(cols, candidates, label):
    for c in candidates:
        if c in cols:
            return c
    raise ValueError(
        f"Required column for '{label}' not found. Tried: {candidates}. "
        f"Got columns: {list(cols)}"
    )


def _to_bool(s: pd.Series) -> pd.Series:
    """Coerce R-style TRUE/FALSE/T/F/0/1/yes/no into Python bools."""
    if s.dtype == bool:
        return s
    str_s = s.astype(str).str.strip().str.lower()
    truthy = {"true", "t", "1", "yes", "y"}
    falsy  = {"false", "f", "0", "no", "n", ""}
    out = pd.Series(np.nan, index=s.index, dtype=object)
    out[str_s.isin(truthy)] = True
    out[str_s.isin(falsy)] = False
    # Anything else (NA, NaN, weird strings) stays NaN -> handled downstream
    return out


# ---------------------------------------------------------------------------
# The decision tree itself
# ---------------------------------------------------------------------------
def nmdetective_b(
    last_exon: bool,
    dist_to_start: float,
    exon_length: float,
    penult_last50: bool,
):
    """
    Apply NMDetective-B to one variant. Returns (score, rule).

    NaN in any input that the active branch needs -> (NaN, "missing_feature").
    Branches not yet reached can have NaN inputs without affecting the result.
    """
    # i. last exon
    if pd.isna(last_exon):
        return np.nan, "missing_feature"
    if bool(last_exon):
        return SCORE_LAST_EXON, "last_exon"

    # ii. start-proximal
    if pd.isna(dist_to_start):
        return np.nan, "missing_feature"
    if dist_to_start < START_PROX_NT:
        return SCORE_START_PROXIMAL, "start_proximal"

    # iii. long exon
    if pd.isna(exon_length):
        return np.nan, "missing_feature"
    if exon_length > LONG_EXON_NT:
        return SCORE_LONG_EXON, "long_exon"

    # iv. last 50 nt of penultimate exon
    if pd.isna(penult_last50):
        return np.nan, "missing_feature"
    if bool(penult_last50):
        return SCORE_50NT_RULE, "50nt_rule"

    # else: triggers NMD
    return SCORE_TRIGGER_NMD, "trigger_NMD"


def classify(score: float) -> str:
    if pd.isna(score):
        return "missing"
    if score > THRESH_SENSITIVE:
        return "sensitive"
    if score < THRESH_ESCAPE:
        return "escape"
    return "intermediate"


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def annotate(variants_path: Path, output_path: Path, qc_path: Path):
    print(f"[load] {variants_path}", file=sys.stderr)
    sep = "," if variants_path.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(variants_path, sep=sep, low_memory=False)
    print(f"[load] {len(df):,} rows, {len(df.columns)} columns", file=sys.stderr)

    le_col = _find_col(df.columns, CANDIDATES["last_exon"],     "last_exon")
    ds_col = _find_col(df.columns, CANDIDATES["dist_to_start"], "dist_to_start")
    el_col = _find_col(df.columns, CANDIDATES["exon_length"],   "exon_length")
    pl_col = _find_col(df.columns, CANDIDATES["penult_last50"], "penult_last50")
    print(f"[cols] last_exon={le_col!r}  dist_to_start={ds_col!r}  "
          f"exon_length={el_col!r}  penult_last50={pl_col!r}", file=sys.stderr)

    last_exon     = _to_bool(df[le_col])
    dist_to_start = pd.to_numeric(df[ds_col], errors="coerce")
    exon_length   = pd.to_numeric(df[el_col], errors="coerce")
    penult_last50 = _to_bool(df[pl_col])

    n = len(df)
    scores = np.full(n, np.nan, dtype=np.float64)
    rules  = np.empty(n, dtype=object)

    for i in range(n):
        s, r = nmdetective_b(
            last_exon.iloc[i],
            dist_to_start.iloc[i],
            exon_length.iloc[i],
            penult_last50.iloc[i],
        )
        scores[i] = s
        rules[i] = r

    classes = pd.Series(scores).map(classify).values

    out = df.copy()
    out["NMDetective_B_score"] = scores
    out["NMDetective_B_rule"]  = rules
    out["NMDetective_B_class"] = classes

    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_sep = "," if output_path.suffix.lower() == ".csv" else "\t"
    out.to_csv(output_path, sep=out_sep, index=False)
    print(f"[write] {output_path}  ({n:,} rows)", file=sys.stderr)

    # QC summary
    rule_counts  = pd.Series(rules).value_counts(dropna=False)
    class_counts = pd.Series(classes).value_counts(dropna=False)
    rows = [("n_variants", n)]
    for k in ("last_exon", "start_proximal", "long_exon", "50nt_rule",
             "trigger_NMD", "missing_feature"):
        rows.append((f"rule__{k}", int(rule_counts.get(k, 0))))
    for k in ("sensitive", "intermediate", "escape", "missing"):
        rows.append((f"class__{k}", int(class_counts.get(k, 0))))
    rows.append(("score_median",
                 float(np.nanmedian(scores)) if np.isfinite(scores).any() else np.nan))
    rows.append(("score_mean",
                 float(np.nanmean(scores)) if np.isfinite(scores).any() else np.nan))
    qc = pd.DataFrame(rows, columns=["metric", "value"])
    qc.to_csv(qc_path, sep="\t", index=False)
    print(f"[write] {qc_path}", file=sys.stderr)
    print("[qc] rule breakdown:", file=sys.stderr)
    for k, v in rule_counts.items():
        print(f"        {k:>20s}  {v:>6d}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(
        description="Apply NMDetective-B decision tree to a variant table "
                    "using local features (no figshare lookup)."
    )
    p.add_argument("--variants", required=True, type=Path,
                   help="Variant table (csv/tsv) with the four NMDetective-B "
                        "features.")
    p.add_argument("--out", required=True, type=Path,
                   help="Output annotated table.")
    p.add_argument("--qc", required=True, type=Path,
                   help="Output QC summary tsv.")
    args = p.parse_args()

    if not args.variants.exists():
        print(f"ERROR: {args.variants} not found", file=sys.stderr)
        sys.exit(2)

    annotate(args.variants, args.out, args.qc)


if __name__ == "__main__":
    main()
