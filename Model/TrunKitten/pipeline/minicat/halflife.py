"""half_life_PC1 merge by version-stripped ENSG, with optional gene-symbol fallback.

Gene-ID reassignment between annotation versions can cause a transcript's
current GTF `gene_id` to differ from the ENSG used to build the half-life
table. When that happens, ENSG-based lookup fails for a genuinely-annotated
gene. We fall back to gene-symbol lookup, which handles these cases.

Example encountered in validation: ENST00000306875 (COG8). In Gencode v26
this transcript's gene_id differs from the ENSG00000213380 used in the
half-life Excel, but both annotate the COG8 gene. Symbol fallback recovers
the correct value.
"""
from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Optional
import logging

import pandas as pd
from .gtf_index import strip_version

log = logging.getLogger(__name__)


class HalfLifeTable:
    def __init__(
        self,
        path: Path,
        sheet: Any = 0,
        ensg_col: str = "gene_id",
        value_col: str = "half_life_PC1",
        symbol_col: Optional[str] = "gene_symbol",
        strip_versions: bool = True,
    ):
        log.info(f"Loading half-life table: {path}")
        df = pd.read_excel(path, sheet_name=sheet, engine="openpyxl")
        for col in (ensg_col, value_col):
            if col not in df.columns:
                raise KeyError(
                    f"Column '{col}' not in Excel. Present: {list(df.columns)}"
                )

        keys = df[ensg_col].astype(str)
        if strip_versions:
            keys = keys.map(strip_version)
        vals = pd.to_numeric(df[value_col], errors="coerce")
        self._ensg_map: Dict[str, float] = dict(zip(keys, vals))

        # Optional gene-symbol fallback map
        self._symbol_map: Dict[str, float] = {}
        if symbol_col is not None and symbol_col in df.columns:
            symbols = df[symbol_col].astype(str).str.strip()
            # For symbol collisions, first value wins; warn if duplicates exist.
            dup_mask = symbols.duplicated(keep="first")
            if dup_mask.any():
                n_dup = int(dup_mask.sum())
                log.warning(
                    f"  {n_dup} duplicate gene symbols in half-life table; "
                    f"first value per symbol will be used as fallback"
                )
            for sym, v in zip(symbols[~dup_mask], vals[~dup_mask]):
                if sym and sym.lower() not in ("nan", "none", ""):
                    self._symbol_map[sym] = v

        log.info(
            f"  loaded {len(self._ensg_map):,} ENSG → half_life_PC1 entries "
            f"({vals.notna().sum():,} non-NaN values)  |  "
            f"{len(self._symbol_map):,} symbol fallback entries"
        )

    def lookup(
        self, ensg: Optional[str], strip_versions: bool = True,
        gene_symbol: Optional[str] = None,
    ) -> Optional[float]:
        """Look up half_life_PC1 by ENSG, with gene-symbol fallback.

        Args:
            ensg: GTF-derived gene_id for the transcript.
            strip_versions: strip trailing '.N' from ENSG before lookup.
            gene_symbol: optional fallback key (from variant input or GTF).
        """
        # Primary: ENSG lookup
        if ensg:
            key = strip_version(ensg) if strip_versions else ensg
            val = self._ensg_map.get(key)
            if val is not None and not pd.isna(val):
                return float(val)

        # Fallback: gene symbol
        if gene_symbol:
            val = self._symbol_map.get(str(gene_symbol).strip())
            if val is not None and not pd.isna(val):
                return float(val)

        return None
