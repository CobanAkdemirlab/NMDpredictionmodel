"""Validate the TrunKitten annotation pipeline against the TrunCat training features.

Loads the first N rows of TOPMed_merged_v4.csv (or any cleaned/merged CSV),
runs the TrunKitten annotation pipeline on those variants, and compares the
10 features column-by-column against the values produced by the TrunCat
training-time feature generation.

Usage:
    python scripts/validate_against_training.py \
        --merged    /path/to/TOPMed_merged_v4.csv \
        --config    config/config.yaml \
        --n         10 \
        --out       outputs/validation_report.tsv

Requires the same inputs as a normal TrunKitten run (GTF, FASTA, BigWigs,
half-life). The underlying implementation package is `minicat`.
"""
from __future__ import annotations
import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

from minicat import REQUIRED_FEATURES
from minicat.config import PipelineConfig
from minicat.gtf_index import build_transcript_index, strip_version
from minicat.conservation import ConservationSource
from minicat.halflife import HalfLifeTable
from minicat.features import FeatureAnnotator


# Column-name map: minicat output → training table.
# All 10 feature names match training except we have to deal with txnames
# being a comma-separated list in the merged file.
MINICAT_TO_TRAINING = {
    "last.EJC": "last.EJC",
    "relativePTClocation": "relativePTClocation",
    "half_life_PC1": "half_life_PC1",
    "cdsseqs_AU_content": "cdsseqs_AU_content",
    "mut.exon": "mut.exon",
    "phastcons_new3utr_first200_median": "phastcons_new3utr_first200_median",
    "phylop_ptc_to_ejc_median": "phylop_ptc_to_ejc_median",
    "AmountExonsAfter": "AmountExonsAfter",
    "cdsseq_AUcontentlast200": "cdsseq_AUcontentlast200",
    "cdsseqs_UC_content": "cdsseqs_UC_content",
}

CATEGORICAL = {"last.EJC"}
INTEGER     = {"mut.exon", "AmountExonsAfter"}


def _setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)s  %(name)s  %(message)s",
        stream=sys.stderr,
    )


def _pick_first_txname(txnames_cell) -> str:
    """TOPMed merged CSV stores txnames as a comma-separated string like
    'ENST00000123,ENST00000456'. Use the first entry — same policy as training."""
    if pd.isna(txnames_cell):
        return ""
    s = str(txnames_cell).strip()
    if not s:
        return ""
    return s.split(",")[0].strip()


def _coerce_bool_last_ejc(v) -> str:
    """Normalise last.EJC values — both sides should be one of
    {upstream, penultimate.last50bp, last.exon}."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "NA"
    return str(v).strip()


def compare_row(our, truth, tol_numeric: float = 1e-3) -> dict:
    """Compare one minicat feature row against one training row."""
    out = {}
    for feat, truth_col in MINICAT_TO_TRAINING.items():
        a = our.get(feat)
        b = truth.get(truth_col)

        if feat in CATEGORICAL:
            a_s, b_s = _coerce_bool_last_ejc(a), _coerce_bool_last_ejc(b)
            out[f"{feat}__ours"]  = a_s
            out[f"{feat}__train"] = b_s
            out[f"{feat}__match"] = (a_s == b_s)
            continue

        a_num = float(a) if a is not None and not (isinstance(a, float) and np.isnan(a)) else np.nan
        b_num = float(b) if b is not None and not (isinstance(b, float) and np.isnan(b)) else np.nan
        out[f"{feat}__ours"]  = a_num
        out[f"{feat}__train"] = b_num

        if np.isnan(a_num) and np.isnan(b_num):
            out[f"{feat}__match"] = True
            out[f"{feat}__abs_diff"] = 0.0
            continue
        if np.isnan(a_num) or np.isnan(b_num):
            out[f"{feat}__match"] = False
            out[f"{feat}__abs_diff"] = np.inf
            continue

        diff = abs(a_num - b_num)
        if feat in INTEGER:
            out[f"{feat}__match"] = (int(round(a_num)) == int(round(b_num)))
        else:
            out[f"{feat}__match"] = (diff <= tol_numeric)
        out[f"{feat}__abs_diff"] = diff
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Validate the TrunKitten annotation pipeline against "
                    "the TrunCat training feature table."
    )
    ap.add_argument("--merged", required=True, help="TOPMed_merged_v4.csv or equivalent")
    ap.add_argument("--config", required=True, help="TrunKitten config YAML")
    ap.add_argument("--n", type=int, default=10, help="number of variants to validate")
    ap.add_argument("--out",    default="outputs/validation_report.tsv")
    ap.add_argument("--tol",    type=float, default=1e-3,
                    help="absolute tolerance for continuous features")
    args = ap.parse_args(argv)

    _setup_logging()
    log = logging.getLogger("validate")

    # 1. Load first N rows of training merged CSV
    log.info(f"Loading {args.n} rows from {args.merged}")
    truth_df = pd.read_csv(args.merged, nrows=args.n)

    # Detect the variant_id column; training CSVs commonly use 'variantID' or a composite key
    if "variant_id" in truth_df.columns:
        vid_col = "variant_id"
    elif "variantID" in truth_df.columns:
        vid_col = "variantID"
    else:
        # fall back to composite key
        truth_df["variant_id"] = (
            truth_df["contig"].astype(str) + "_"
            + truth_df["position"].astype(str) + "_"
            + truth_df["refAllele"].astype(str) + "_"
            + truth_df["altAllele"].astype(str)
        )
        vid_col = "variant_id"

    # 2. Build a variants-only DataFrame in our pipeline's schema
    needed = ["contig", "position", "refAllele", "altAllele", "gene", "txnames"]
    missing = [c for c in needed if c not in truth_df.columns]
    if missing:
        raise KeyError(f"Merged CSV missing columns: {missing}")

    variants_df = truth_df[[vid_col] + needed].copy()
    variants_df = variants_df.rename(columns={vid_col: "variant_id"})
    # Pick the first txname (training policy is to annotate with a single tx)
    variants_df["txnames"] = variants_df["txnames"].map(_pick_first_txname)
    log.info(f"Prepared {len(variants_df)} variants for annotation")
    log.info(f"  first txnames: {variants_df['txnames'].head().tolist()}")

    # 3. Run minicat on those variants — do this inline rather than writing to a file
    cfg = PipelineConfig.from_yaml(args.config)

    wanted = set(variants_df["txnames"].map(strip_version).dropna().unique())
    tx_index = build_transcript_index(cfg.gtf, wanted_tx_ids_stripped=wanted)

    fasta     = Fasta(str(cfg.fasta), as_raw=False, sequence_always_upper=True)
    phastcons = ConservationSource(cfg.phastcons, chr_style=cfg.chr_style)
    phylop    = ConservationSource(cfg.phylop,    chr_style=cfg.chr_style)
    hl        = HalfLifeTable(cfg.halflife, sheet=cfg.halflife_sheet,
                              ensg_col=cfg.halflife_ensg_col,
                              value_col=cfg.halflife_value_col,
                              symbol_col=cfg.halflife_symbol_col or None,
                              strip_versions=cfg.strip_versions)

    annot = FeatureAnnotator(
        tx_index=tx_index, fasta=fasta,
        phastcons=phastcons, phylop=phylop, halflife=hl,
        strip_versions=cfg.strip_versions,
        cds_last_window=cfg.cds_last_window,
        new3utr_window=cfg.new3utr_window,
    )

    our_rows = [annot.annotate(r).to_feature_row()
                for r in variants_df.to_dict(orient="records")]
    our_df = pd.DataFrame(our_rows)

    phastcons.close(); phylop.close()

    # 4. Row-by-row comparison (aligned by positional index)
    rows = []
    for i in range(len(variants_df)):
        ours  = our_df.iloc[i].to_dict()
        truth = truth_df.iloc[i].to_dict()
        cmp = compare_row(ours, truth, tol_numeric=args.tol)
        cmp = {
            "variant_id": ours.get("variant_id"),
            "txnames": ours.get("txnames"),
            **cmp,
        }
        rows.append(cmp)

    report = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.out, sep="\t", index=False)

    # 5. Summary
    log.info("="*72)
    log.info(f"VALIDATION SUMMARY (n={len(report)})")
    log.info("="*72)
    for feat in MINICAT_TO_TRAINING.keys():
        match_col = f"{feat}__match"
        diff_col  = f"{feat}__abs_diff"
        n_match = int(report[match_col].sum())
        pct = n_match / len(report) * 100
        line = f"  {feat:<40s}  {n_match}/{len(report)} match ({pct:.0f}%)"
        if diff_col in report.columns:
            finite = report[diff_col].replace([np.inf], np.nan).dropna()
            if len(finite):
                line += f"  max|Δ|={finite.max():.4g}  median|Δ|={finite.median():.4g}"
        log.info(line)

    log.info(f"\nFull report: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
