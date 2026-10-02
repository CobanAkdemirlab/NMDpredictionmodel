"""MedianExpression_log2 from GTEx v8 median gene-level TPM.

Mirrors Features/04_gene_level_features.R exactly:

    expression$gene_id <- sub("\\..*", "", expression$Name)          # strip version
    expression$MedianExpression <- apply(tissue_cols, 1, median, na.rm = TRUE)
    MedianExpression_log2 <- log2(MedianExpression + 1)

i.e. for each gene, take the MEDIAN ACROSS TISSUES of the per-tissue median TPM
(all tissue columns, NaN ignored), then log2(x + 1).

Input file (GTEx Analysis V8 -> RNA-Seq -> "Median gene-level TPM by tissue"):
    GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
Both the original GCT (two header lines, "#1.2" and the dimensions line) and a
plain tab-separated table whose first line is the header are accepted.

GTEx v8 is annotated with GENCODE v26 gene IDs, the same release as the GTF used
for the rest of the pipeline, so lookup is by version-stripped ENSG. As with
half_life_PC1, a gene-symbol fallback (GTEx `Description` column) covers genes
whose ID differs between annotation versions; set `symbol_fallback=False` to
disable it and use IDs only.

Missing genes return None; predict.py then applies the training median
(MedianExpression_log2, from training_medians.json), as TrunCat's notebook 02
does.
"""
from __future__ import annotations
import gzip
import logging
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

from .gtf_index import strip_version

log = logging.getLogger(__name__)

_ID_COLS = ("Name", "Description")


def _header_rows_to_skip(path: Path) -> int:
    """GCT files start with '#1.2' and a dimensions line; plain TSVs do not."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        first = f.readline()
    return 2 if first.startswith("#1.2") else 0


class GTExExpressionTable:
    def __init__(
        self,
        path: Path,
        strip_versions: bool = True,
        symbol_fallback: bool = True,
    ):
        path = Path(path)
        log.info(f"Loading GTEx median-TPM table: {path}")
        df = pd.read_csv(path, sep="\t", skiprows=_header_rows_to_skip(path))
        for col in _ID_COLS:
            if col not in df.columns:
                raise KeyError(
                    f"Column '{col}' not in GTEx table. Present: {list(df.columns)[:6]} ..."
                )
        tissue_cols = [c for c in df.columns if c not in _ID_COLS]
        if len(tissue_cols) < 10:
            raise ValueError(
                f"Only {len(tissue_cols)} tissue columns found in {path}; expected the "
                f"GTEx v8 median-TPM table (54 tissues)."
            )

        tpm = df[tissue_cols].apply(pd.to_numeric, errors="coerce")
        median_tpm = tpm.median(axis=1, skipna=True)          # across tissues, NaN ignored
        log2_expr = np.log2(median_tpm + 1.0)

        keys = df["Name"].astype(str)
        if strip_versions:
            keys = keys.map(strip_version)

        self._ensg_map: Dict[str, float] = {}
        n_dup = 0
        for k, v in zip(keys, log2_expr):
            if k in self._ensg_map:
                n_dup += 1               # e.g. PAR-region genes collapse to one stripped ID
                continue
            self._ensg_map[k] = float(v) if pd.notna(v) else float("nan")
        if n_dup:
            log.warning(f"  {n_dup} GTEx rows share a version-stripped ENSG; first row kept")

        self._symbol_map: Dict[str, float] = {}
        if symbol_fallback:
            symbols = df["Description"].astype(str).str.strip()
            seen = set()
            for sym, v in zip(symbols, log2_expr):
                if sym and sym.lower() not in ("nan", "none", "") and sym not in seen:
                    seen.add(sym)
                    self._symbol_map[sym] = float(v) if pd.notna(v) else float("nan")

        n_valid = int(np.isfinite(list(self._ensg_map.values())).sum())
        log.info(
            f"  {len(self._ensg_map):,} genes ({n_valid:,} with a value) across "
            f"{len(tissue_cols)} tissues | {len(self._symbol_map):,} symbol fallback entries"
        )

    def lookup(
        self,
        ensg: Optional[str],
        strip_versions: bool = True,
        gene_symbol: Optional[str] = None,
    ) -> Optional[float]:
        """Return log2(median-across-tissues TPM + 1), or None if the gene is absent."""
        if ensg:
            key = strip_version(ensg) if strip_versions else ensg
            val = self._ensg_map.get(key)
            if val is not None and not pd.isna(val):
                return float(val)
        if gene_symbol and self._symbol_map:
            val = self._symbol_map.get(str(gene_symbol).strip())
            if val is not None and not pd.isna(val):
                return float(val)
        return None
