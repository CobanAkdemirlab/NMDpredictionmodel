"""Transcript sequence reconstruction and base-composition features.

Reference FASTA is DNA (A/C/G/T). We build transcript-oriented sequences by:
  + strand: concatenating exon/CDS blocks in ascending genomic order.
  - strand: reverse-complementing each block and concatenating in descending
    genomic order so the result is 5'→3' in transcript orientation.

AU "content" is computed on the DNA alphabet (AU = A+T); naming is kept
as-trained. N bases are excluded.
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


def compute_cds_composition(fasta: Fasta, tx: TranscriptRecord) -> dict:
    """Compute the CDS-level composition feature kept in TrunKitten's 8-feature set.

    Returns dict with:
        cdsseqs_AU_content — AU content of the full CDS
        cds_length         — for QC only

    `cdsseqs_UC_content` and `cdsseq_AUcontentlast200` were dropped from the
    TrunKitten feature set (Sept 2026 CV-protocol correction; see
    trunkitten_features.json) and are no longer computed here. If either is
    ever needed again (e.g. for TrunCat, which still uses all three), restore
    `_uc_content` and the last-200nt slice from git history.
    """
    cds = cds_sequence(fasta, tx)
    L = len(cds)
    if L == 0:
        return {
            "cdsseqs_AU_content": float("nan"),
            "cds_length": 0,
        }
    return {
        "cdsseqs_AU_content": _au_content(cds),
        "cds_length": L,
    }