"""Config loading and validation."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import yaml


@dataclass(frozen=True)
class PipelineConfig:
    # paths
    variants:   Path
    gtf:        Path
    fasta:      Path
    phastcons:  Path
    halflife:   Path
    gtex:       Path   # GTEx v8 median gene-level TPM by tissue (.gct / .gct.gz / .tsv)
    halflife_sheet: Any
    halflife_ensg_col:   str
    halflife_value_col:  str
    halflife_symbol_col: str  # optional; empty string disables fallback

    out_features: Path
    out_qc:       Path
    log:          Path

    # options
    chr_style:       str  = "auto"
    strip_versions:  bool = True
    n_workers:       int  = 4
    new3utr_window:  int  = 200
    expression_symbol_fallback: bool = True   # gene-symbol fallback for GTEx ID drift

    @classmethod
    def from_yaml(cls, path: str | Path) -> "PipelineConfig":
        path = Path(path)
        with open(path) as f:
            raw = yaml.safe_load(f)
        p = raw["paths"]
        o = raw.get("options", {})

        def _pth(key: str) -> Path:
            if key not in p:
                raise KeyError(
                    f"config 'paths' is missing '{key}'. "
                    + ("TrunKitten now needs a GTEx v8 median-TPM table for MedianExpression_log2; "
                       "add `gtex: inputs/gtex_gene_median_tpm.gct.gz` (see README)." if key == "gtex" else "")
                )
            return Path(p[key]).expanduser()

        cfg = cls(
            variants   = _pth("variants"),
            gtf        = _pth("gtf"),
            fasta      = _pth("fasta"),
            phastcons  = _pth("phastcons"),
            halflife   = _pth("halflife"),
            gtex       = _pth("gtex"),
            halflife_sheet     = p.get("halflife_sheet", 0),
            halflife_ensg_col  = p.get("halflife_ensg_col", "gene_id"),
            halflife_value_col = p.get("halflife_value_col", "half_life_PC1"),
            halflife_symbol_col = p.get("halflife_symbol_col", "gene_symbol"),
            out_features = Path(p["out_features"]).expanduser(),
            out_qc       = Path(p["out_qc"]).expanduser(),
            log          = Path(p["log"]).expanduser(),
            chr_style       = o.get("chr_style", "auto"),
            strip_versions  = o.get("strip_versions", True),
            n_workers       = int(o.get("n_workers", 4)),
            new3utr_window  = int(o.get("new3utr_window", 200)),
            expression_symbol_fallback = bool(o.get("expression_symbol_fallback", True)),
        )
        cfg._validate_inputs()
        cfg.out_features.parent.mkdir(parents=True, exist_ok=True)
        cfg.log.parent.mkdir(parents=True, exist_ok=True)
        return cfg

    def _validate_inputs(self) -> None:
        missing = []
        for attr in ("variants", "gtf", "fasta", "phastcons", "halflife", "gtex"):
            pth = getattr(self, attr)
            if not pth.exists():
                missing.append(f"{attr}={pth}")
        if missing:
            raise FileNotFoundError("Missing input files: " + "; ".join(missing))
