"""Transcript sequence reconstruction and base-composition features.

Reference FASTA is DNA (A/C/G/T). We build transcript-oriented sequences by:
  + strand: concatenating exon/CDS blocks in ascending genomic order.
  - strand: reverse-complementing each block and concatenating in descending
    genomic order so the result is 5'→3' in transcript orientation.

AU / UC "content" is computed on the DNA alphabet (AU = A+T, UC = T+C);
naming is kept as-trained. N bases are excluded.
"""
from __future__ import annotations
from typing import List, Tuple, Optional
from pyfaidx import Fasta

from .gtf_index import TranscriptRecord


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def fetch_blocks(
    fasta: Fasta, chrom: str, strand: str, blocks_1based: List[Tuple[int, int]],
) -> str:
    """Fetch spliced transcript-oriented sequence from a list of 1-based inclusive blocks.

    blocks_1based must be supplied in transcript 5'→3' order (i.e., already
    sorted by gtf_index).
    """
    parts: List[str] = []
    for (s, e) in blocks_1based:
        # pyfaidx slicing is 0-based half-open
        raw = str(fasta[chrom][s - 1:e]).upper()
        if strand == "-":
            raw = revcomp(raw)
        parts.append(raw)
    return "".join(parts)


def _au_content(seq: str) -> float:
    """Match R Biostrings::alphabetFrequency(baseOnly=TRUE, as.prob=TRUE):
    denominator is the full sequence length, INCLUDING any N/other bases.
    """
    n = len(seq)
    if n == 0:
        return float("nan")
    au = sum(1 for b in seq if b == "A" or b == "T")
    return au / n


def _uc_content(seq: str) -> float:
    n = len(seq)
    if n == 0:
        return float("nan")
    uc = sum(1 for b in seq if b == "T" or b == "C")
    return uc / n


def cds_sequence(fasta: Fasta, tx: TranscriptRecord) -> str:
    """Full CDS string INCLUDING the stop codon, transcript-oriented.

    Matches Iman's R pipeline: `extractTranscriptSeqs(cdsBy(makeTxDbFromGFF(...)))`.
    Bioconductor's `cdsBy()` extends CDS ranges to include stop_codon records
    from the GTF. We replicate this by concatenating GTF CDS + stop_codon
    records in transcript order via `tx.cds_with_stop()`.
    """
    blocks = tx.cds_with_stop()
    if not blocks:
        return ""
    return fetch_blocks(fasta, tx.chrom, tx.strand, blocks)


def compute_cds_composition(
    fasta: Fasta, tx: TranscriptRecord, last_window: int = 200,
) -> dict:
    """Compute all three CDS-level composition features.

    Returns dict with:
        cdsseqs_AU_content        — AU content of full CDS (always computed)
        cdsseqs_UC_content        — UC content of full CDS (always computed)
        cdsseq_AUcontentlast200   — AU content of last `last_window` nt of CDS,
                                    or NaN if CDS shorter than window
        cds_length                — for QC
        cds_length_short_flag     — True if CDS < last_window

    Matches Iman's R:
        if (length(a)<200){ return(NA) }
    so the last-window feature is NaN when CDS < window. Consumer applies
    training zero-fill (Notebook 02 zero-filled these 5 rows).
    """
    cds = cds_sequence(fasta, tx)
    L = len(cds)
    if L == 0:
        return {
            "cdsseqs_AU_content":      float("nan"),
            "cdsseqs_UC_content":      float("nan"),
            "cdsseq_AUcontentlast200": float("nan"),
            "cds_length": 0,
            "cds_length_short_flag": True,
        }

    short_flag = L < last_window
    last200_au = _au_content(cds[-last_window:]) if not short_flag else float("nan")

    return {
        "cdsseqs_AU_content":      _au_content(cds),
        "cdsseqs_UC_content":      _uc_content(cds),
        "cdsseq_AUcontentlast200": last200_au,
        "cds_length": L,
        "cds_length_short_flag": short_flag,
    }
