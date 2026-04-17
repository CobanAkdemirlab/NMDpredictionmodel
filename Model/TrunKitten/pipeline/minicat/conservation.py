"""Conservation median extraction from BigWig or BED/bedGraph input.

Medians are computed over the concatenation of valid (non-NaN) per-base
values across all genomic segments in a region — matching
`conservation_score_extraction_v2.py`.
"""
from __future__ import annotations
from pathlib import Path
from typing import List, Optional, Tuple
import logging
import numpy as np

log = logging.getLogger(__name__)


class ConservationSource:
    """Wraps a BigWig or BED-like file with a uniform .median_over_blocks API."""

    def __init__(self, path: Path, chr_style: str = "auto"):
        self.path = Path(path)
        self.chr_style = chr_style
        self._kind = self._detect_kind()
        self._bw = None
        self._bed_df = None
        self._bw_chroms = None
        self._open()

    def _detect_kind(self) -> str:
        sfx = "".join(self.path.suffixes).lower()
        if sfx.endswith(".bw") or sfx.endswith(".bigwig"):
            return "bigwig"
        if sfx.endswith(".bed") or sfx.endswith(".bed.gz") or sfx.endswith(".bedgraph") or sfx.endswith(".bedgraph.gz"):
            return "bed"
        # default: try bigwig
        return "bigwig"

    def _open(self):
        if self._kind == "bigwig":
            import pyBigWig
            self._bw = pyBigWig.open(str(self.path))
            if self._bw is None:
                raise RuntimeError(f"Failed to open BigWig: {self.path}")
            self._bw_chroms = set(self._bw.chroms().keys())
            log.info(f"  opened BigWig {self.path.name}  ({len(self._bw_chroms)} seqs)")
        else:
            import pandas as pd
            # bedGraph: chrom, start0, end0, value
            df = pd.read_csv(self.path, sep="\t", comment="#", header=None,
                             names=["chrom", "start", "end", "value"])
            df["start"] = df["start"].astype(int)
            df["end"]   = df["end"].astype(int)
            df["value"] = df["value"].astype(float)
            self._bed_df = df
            log.info(f"  loaded BED/bedGraph {self.path.name}  ({len(df):,} intervals)")

    def close(self):
        if self._bw is not None:
            self._bw.close()

    # -------- chromosome-name harmonisation --------
    def _harmonise_chrom(self, chrom: str) -> Optional[str]:
        """Return the chromosome name that matches the data source, or None."""
        if self._kind == "bigwig":
            if chrom in self._bw_chroms:
                return chrom
            # try toggling "chr"
            alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
            if alt in self._bw_chroms:
                return alt
            return None
        else:
            # BED: no precheck — pandas filtering later handles missing chroms fine
            return chrom

    # -------- per-region median --------
    def median_over_blocks(
        self, blocks_1based: List[Tuple[str, int, int]],
    ) -> Tuple[float, int, int]:
        """Compute median of valid values over the union of 1-based inclusive blocks.

        Returns (median_or_nan, total_bp, valid_bp).
        """
        if not blocks_1based:
            return (float("nan"), 0, 0)

        concat_arrays = []
        total_bp = 0
        for (chrom, s1, e1) in blocks_1based:
            if e1 < s1:
                continue
            c_harm = self._harmonise_chrom(chrom)
            if c_harm is None:
                continue
            # BigWig uses 0-based half-open
            s0, e0 = s1 - 1, e1
            total_bp += (e0 - s0)

            if self._kind == "bigwig":
                try:
                    arr = self._bw.values(c_harm, s0, e0, numpy=True)
                except (RuntimeError, Exception) as ex:
                    log.debug(f"  BigWig fetch failed {c_harm}:{s0}-{e0}: {ex}")
                    continue
                concat_arrays.append(np.asarray(arr, dtype=float))
            else:
                vec = self._bed_vector(c_harm, s0, e0)
                if vec.size:
                    concat_arrays.append(vec)

        if not concat_arrays:
            return (float("nan"), total_bp, 0)

        full = np.concatenate(concat_arrays)
        valid = full[~np.isnan(full)]
        if valid.size == 0:
            return (float("nan"), total_bp, 0)
        return (float(np.median(valid)), total_bp, int(valid.size))

    def _bed_vector(self, chrom: str, s0: int, e0: int) -> np.ndarray:
        """Build a per-base value vector over [s0, e0) from the BED/bedGraph DF."""
        df = self._bed_df
        sub = df[(df["chrom"] == chrom) & (df["end"] > s0) & (df["start"] < e0)]
        if sub.empty:
            return np.full(e0 - s0, np.nan)
        vec = np.full(e0 - s0, np.nan)
        for _, row in sub.iterrows():
            i0 = max(0, int(row["start"]) - s0)
            i1 = min(e0 - s0, int(row["end"]) - s0)
            if i1 > i0:
                vec[i0:i1] = row["value"]
        return vec
