"""TrunKitten annotation CLI — command-line entry point.

Annotates PTC / stop-gain variants with the 10 features required by
TrunKitten, the reduced NMD-prediction model derived from TrunCat.

Usage:
    python -m minicat.cli --config config/config.yaml

Output:
    outputs/annotated.tsv   — 10-feature table (one row per input variant)
    outputs/qc_report.tsv   — QC columns per variant
    outputs/run.log         — log file
"""
from __future__ import annotations
import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from pyfaidx import Fasta

from .config import PipelineConfig
from .gtf_index import build_transcript_index, strip_version
from .conservation import ConservationSource
from .halflife import HalfLifeTable
from .features import FeatureAnnotator


def _setup_logging(logfile: Path, level: str = "INFO") -> None:
    logfile.parent.mkdir(parents=True, exist_ok=True)
    fmt = "%(asctime)s  %(levelname)s  %(name)s  %(message)s"
    handlers = [
        logging.StreamHandler(sys.stderr),
        logging.FileHandler(logfile, mode="w"),
    ]
    logging.basicConfig(level=level, format=fmt, handlers=handlers)


def _load_variants(path: Path) -> pd.DataFrame:
    sep = "\t" if str(path).endswith((".tsv", ".txt", ".tsv.gz")) else ","
    df = pd.read_csv(path, sep=sep)
    required = {"variant_id", "contig", "position", "refAllele", "altAllele", "gene", "txnames"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Input variants file missing required columns: {sorted(missing)}")
    return df


def run(cfg: PipelineConfig) -> None:
    log = logging.getLogger("minicat")
    log.info("="*72)
    log.info(f"TrunKitten annotation pipeline (package: minicat)")
    log.info("="*72)

    # 1. Variants
    variants = _load_variants(cfg.variants)
    log.info(f"Loaded {len(variants):,} variants from {cfg.variants}")

    # 2. GTF index (filter by wanted tx for speed)
    wanted = set(
        variants["txnames"].astype(str).map(strip_version).dropna().unique()
    )
    tx_index = build_transcript_index(cfg.gtf, wanted_tx_ids_stripped=wanted)

    # 3. FASTA
    fasta = Fasta(str(cfg.fasta), as_raw=False, sequence_always_upper=True)
    log.info(f"Opened FASTA: {cfg.fasta}  ({len(fasta.keys())} seqs)")

    # 4. Conservation sources
    phastcons = ConservationSource(cfg.phastcons, chr_style=cfg.chr_style)
    phylop    = ConservationSource(cfg.phylop,    chr_style=cfg.chr_style)

    # 5. Half-life table
    hl = HalfLifeTable(
        cfg.halflife,
        sheet=cfg.halflife_sheet,
        ensg_col=cfg.halflife_ensg_col,
        value_col=cfg.halflife_value_col,
        symbol_col=cfg.halflife_symbol_col or None,
        strip_versions=cfg.strip_versions,
    )

    # 6. Annotator
    annot = FeatureAnnotator(
        tx_index=tx_index,
        fasta=fasta,
        phastcons=phastcons,
        phylop=phylop,
        halflife=hl,
        strip_versions=cfg.strip_versions,
        cds_last_window=cfg.cds_last_window,
        new3utr_window=cfg.new3utr_window,
    )

    # 7. Iterate (single-threaded; BigWig handles are not process-safe, and
    #    per-variant work is I/O-dominated — for large cohorts, use a thread
    #    pool or chunk by chromosome).
    feat_rows, qc_rows = [], []
    for i, r in enumerate(variants.to_dict(orient="records"), start=1):
        res = annot.annotate(r)
        feat_rows.append(res.to_feature_row())
        qc_rows.append(res.to_qc_row())
        if i % 500 == 0:
            log.info(f"  annotated {i:,} / {len(variants):,}")

    # 8. Write outputs
    feat_df = pd.DataFrame(feat_rows)
    qc_df   = pd.DataFrame(qc_rows)
    feat_df.to_csv(cfg.out_features, sep="\t", index=False)
    qc_df.to_csv(cfg.out_qc, sep="\t", index=False)

    # 9. Summary
    log.info("="*72)
    log.info(f"Wrote features: {cfg.out_features}  ({len(feat_df):,} rows)")
    log.info(f"Wrote QC:       {cfg.out_qc}")
    log.info(f"  tx_not_in_gtf:            {qc_df['tx_not_in_gtf'].sum():,}")
    log.info(f"  version_mismatch:         {qc_df['version_mismatch'].sum():,}")
    log.info(f"  boundary_ambiguous:       {qc_df['boundary_ambiguous'].sum():,}")
    log.info(f"  cds_length_short_flag:    {qc_df['cds_length_short_flag'].sum():,}")
    log.info(f"  new3utr_empty:            {qc_df['new3utr_empty'].sum():,}")
    log.info(f"  ptc_to_ejc_empty:         {qc_df['ptc_to_ejc_empty'].sum():,}")
    log.info(f"  half_life_missing:        {qc_df['half_life_missing'].sum():,}")
    log.info(f"  any_conservation_missing: {qc_df['any_conservation_missing'].sum():,}")

    phastcons.close()
    phylop.close()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="TrunKitten PTC annotation pipeline "
                    "(reduced top-10 feature NMD model; package: minicat)"
    )
    ap.add_argument("--config", required=True, help="path to config YAML")
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args(argv)

    cfg = PipelineConfig.from_yaml(args.config)
    _setup_logging(cfg.log, level=args.log_level)
    run(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
