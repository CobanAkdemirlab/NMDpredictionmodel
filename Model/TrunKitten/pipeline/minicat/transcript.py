"""Transcript-coordinate math: position mapping, exon location, region building.

All public functions take a TranscriptRecord whose exons list is already
sorted in transcript 5'->3' order (gtf_index does this at build time).
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Tuple

from .gtf_index import TranscriptRecord


@dataclass
class PTCLocation:
    """Where a PTC falls within a transcript."""
    exon_rank: int            # 1-based transcript-ordered exon number
    exon_count: int
    transcript_pos: int       # 1-based within spliced transcript
    transcript_length: int
    within_exon_offset: int   # 1-based from exon's 5' end (transcript-oriented)
    exon_length: int
    boundary_ambiguous: bool  # PTC sits exactly at exon start or end


def locate_ptc_in_transcript(
    tx: TranscriptRecord, ptc_genomic_pos: int,
) -> Optional[PTCLocation]:
    """Find which transcript exon contains the PTC genomic position.

    Args:
        tx: TranscriptRecord (exons sorted in transcript order).
        ptc_genomic_pos: 1-based genomic coordinate of the mutated base.

    Returns:
        PTCLocation, or None if the position is not within any exon of this
        transcript. Boundary case: if position equals exon_start or exon_end
        of any exon, returns that exon with boundary_ambiguous=True.
    """
    if not tx.exons:
        return None

    tx_length = tx.transcript_length
    cumulative = 0

    for rank, (gstart, gend) in enumerate(tx.exons, start=1):
        # exon interval is [gstart, gend] inclusive
        if gstart <= ptc_genomic_pos <= gend:
            exon_len = gend - gstart + 1
            if tx.strand == "+":
                within = ptc_genomic_pos - gstart + 1   # 1-based from 5' end
            else:
                within = gend - ptc_genomic_pos + 1     # 5' end of minus-strand exon is gend
            boundary = (ptc_genomic_pos == gstart) or (ptc_genomic_pos == gend)
            return PTCLocation(
                exon_rank=rank,
                exon_count=len(tx.exons),
                transcript_pos=cumulative + within,
                transcript_length=tx_length,
                within_exon_offset=within,
                exon_length=exon_len,
                boundary_ambiguous=boundary,
            )
        cumulative += (gend - gstart + 1)
    return None


def amount_exons_after(tx: TranscriptRecord, ptc_exon_rank: int) -> int:
    """Number of *coding* exons strictly after the PTC-containing exon.

    "Coding" = exon overlaps any CDS record for this transcript. The
    PTC-containing exon is NOT counted, regardless of whether it is fully
    or partially coding.
    """
    coding_flags = tx.coding_exon_flags()
    # exon ranks are 1-based; count coding exons at ranks > ptc_exon_rank
    return sum(1 for i, c in enumerate(coding_flags, start=1)
               if c and i > ptc_exon_rank)


def last_ejc_category(tx: TranscriptRecord, loc: PTCLocation) -> str:
    """Return one of {'last.exon', 'penultimate.last50bp', 'upstream'}.

    Matches Iman's `Variant_annotation.R` `penultimate.exon` rule exactly:

        if (coding_pos >= fifty_pos &&
            coding_pos <= penultimate_end &&
            penultimate_length >= 50) {
          return("penultimate.last50bp")
        }

    where `fifty_pos = penultimate_end - 50`. In Python equivalents:
      - PTC must be in the penultimate coding exon (coding rank == n_coding - 1)
      - distance from PTC to 3' end of penultimate's CDS portion must be ≤ 50
      - the CDS portion of the penultimate coding exon must itself be ≥ 50 bp
        (otherwise the variant is 'upstream', even if numerically close to the
        3' end — this guards against tiny penultimate exons where the whole
        exon is effectively "the last 50 bp region")
    """
    coding_flags = tx.coding_exon_flags()
    n_coding = sum(coding_flags)
    if n_coding == 0:
        return "upstream"  # defensive

    c_rank = coding_exon_rank(tx, loc.exon_rank)
    if c_rank is None:
        return "upstream"  # PTC exon not coding (shouldn't happen for stop-gains)

    if n_coding == 1:
        return "last.exon"
    if c_rank == n_coding:
        return "last.exon"
    if c_rank == n_coding - 1:
        ptc_pos = _ptc_pos_from_loc(tx, loc)
        if ptc_pos is None:
            return "upstream"
        d = _distance_ptc_to_coding_exon_3p(tx, ptc_pos, loc)
        if d is None or d > 50:
            return "upstream"
        # Length guard: the CDS portion of the penultimate coding exon
        # must itself be ≥ 50 bp, matching Iman's `penultimate_length >= 50`.
        pen_len = _penultimate_cds_length(tx)
        if pen_len is None or pen_len < 50:
            return "upstream"
        return "penultimate.last50bp"
    return "upstream"


def _penultimate_cds_length(tx: TranscriptRecord) -> Optional[int]:
    """CDS-with-stop length of the penultimate coding exon, in transcript order.

    Matches Iman's convention: sum of CDS+stop_codon segment lengths that
    overlap the penultimate coding exon.
    """
    coding_flags = tx.coding_exon_flags()
    # Find transcript-exon indices of the coding exons
    coding_tx_indices = [i for i, c in enumerate(coding_flags) if c]
    if len(coding_tx_indices) < 2:
        return None
    pen_idx = coding_tx_indices[-2]  # penultimate coding exon's tx-exon index
    egs, ege = tx.exons[pen_idx]

    total = 0
    for (cs, ce) in tx.cds_with_stop():
        o_s = max(cs, egs)
        o_e = min(ce, ege)
        if o_e >= o_s:
            total += (o_e - o_s + 1)
    return total if total > 0 else None


def _ptc_pos_from_loc(tx: TranscriptRecord, loc: PTCLocation) -> Optional[int]:
    """Recover the genomic PTC position from a PTCLocation.

    PTCLocation stores the transcript-oriented within-exon offset but not
    the raw genomic coordinate. We reconstruct it using the PTC exon's
    genomic bounds and the strand.
    """
    if loc.exon_rank < 1 or loc.exon_rank > len(tx.exons):
        return None
    gstart, gend = tx.exons[loc.exon_rank - 1]
    if tx.strand == "+":
        return gstart + loc.within_exon_offset - 1
    else:
        return gend - loc.within_exon_offset + 1


def _distance_ptc_to_coding_exon_3p(
    tx: TranscriptRecord, ptc_genomic_pos: int, loc: PTCLocation,
) -> Optional[int]:
    """Transcript-oriented distance (in bp) from PTC to the 3' end of the
    CDS-with-stop portion of the PTC-containing exon.

    Matches Iman's R pipeline convention: she operates on `cdsBy()` blocks
    (CDS + stop_codon), so the "penultimate exon" has its 3' boundary at
    the end of the CDS+stop portion within that exon, NOT at the full
    genomic exon boundary. This matters for transcripts where the native
    stop sits mid-penultimate-exon with trailing 3'UTR — my earlier
    exon-boundary measurement over-classified these as `penultimate.last50bp`
    when Iman classified them `upstream`.

    Returns the distance in bp, or None if the PTC is outside the
    CDS-with-stop extent (shouldn't happen for a real stop-gain).
    """
    # Find the CDS-with-stop 3' end *within the PTC-containing exon*.
    cws_blocks = tx.cds_with_stop()
    if not cws_blocks:
        return None

    # Genomic bounds of the PTC exon
    egs, ege = tx.exons[loc.exon_rank - 1]

    # CDS-with-stop 3' endpoint within this exon, measured at the exon level
    # (strand-aware): the furthest CDS/stop-codon boundary that falls within
    # [egs, ege] on the 3' side of the PTC.
    cws_3p_in_exon: Optional[int] = None
    for (cs, ce) in cws_blocks:
        # overlap with the exon
        o_s = max(cs, egs)
        o_e = min(ce, ege)
        if o_e < o_s:
            continue
        if tx.strand == "+":
            # on +strand, "3' end of CDS-in-exon" is the maximum overlap end
            cws_3p_in_exon = o_e if cws_3p_in_exon is None else max(cws_3p_in_exon, o_e)
        else:
            # on -strand, "3' end" in transcript orientation is the minimum genomic start
            cws_3p_in_exon = o_s if cws_3p_in_exon is None else min(cws_3p_in_exon, o_s)

    if cws_3p_in_exon is None:
        return None

    if tx.strand == "+":
        # distance on +strand: cws_3p_in_exon - ptc (in transcript orientation)
        if cws_3p_in_exon < ptc_genomic_pos:
            return 0  # PTC is at or past the CDS-with-stop end; distance clipped to 0
        return cws_3p_in_exon - ptc_genomic_pos
    else:
        if cws_3p_in_exon > ptc_genomic_pos:
            return 0
        return ptc_genomic_pos - cws_3p_in_exon


def coding_exon_rank(tx: TranscriptRecord, transcript_exon_rank: int) -> Optional[int]:
    """Convert a whole-transcript exon rank to a coding-exon rank.

    "Coding exon" = exon that overlaps any CDS segment for this transcript.
    Returns the 1-based index of the PTC-containing exon among coding exons,
    or None if the PTC-containing exon is not itself coding.

    This matches the training convention where `mut.exon` appears to be
    numbered within the CDS (skipping 5'UTR-only first exons), which
    produces off-by-one differences vs. whole-transcript numbering when
    the 5'UTR occupies a dedicated exon.
    """
    coding_flags = tx.coding_exon_flags()
    if transcript_exon_rank < 1 or transcript_exon_rank > len(coding_flags):
        return None
    if not coding_flags[transcript_exon_rank - 1]:
        return None
    # count coding exons up to and including this one
    return sum(1 for i in range(transcript_exon_rank) if coding_flags[i])


def ptc_cds_position(tx: TranscriptRecord, ptc_genomic_pos: int) -> Optional[int]:
    """1-based CDS-coordinate position of the PTC.

    Walks CDS segments in transcript order and sums lengths until the
    PTC-containing CDS segment; returns cumulative CDS length at PTC.
    Returns None if PTC is not within any CDS segment (e.g. PTC in UTR —
    shouldn't happen for a stop-gain).
    """
    if not tx.cds:
        return None
    cumulative = 0
    for (gstart, gend) in tx.cds:
        if gstart <= ptc_genomic_pos <= gend:
            if tx.strand == "+":
                within = ptc_genomic_pos - gstart + 1
            else:
                within = gend - ptc_genomic_pos + 1
            return cumulative + within
        cumulative += (gend - gstart + 1)
    return None


def relative_ptc_location(tx: TranscriptRecord, ptc_genomic_pos: int) -> float:
    """PTC CDS-position / CDS length (stop-codon-inclusive).

    Training convention: `relativePTClocation <- coding.pos / cds_length`
    where `cds_length` comes from Bioconductor `extractTranscriptSeqs(cdsBy(...))`,
    which extends CDS ranges to include the stop_codon records. So the
    denominator is `sum(cds lengths) + 3` (for a standard 3-bp stop).

    The numerator `coding.pos` is the PTC's 1-based position within the CDS
    (not including stop codon, since the PTC is by definition before the
    native stop).
    """
    cds_pos = ptc_cds_position(tx, ptc_genomic_pos)
    cds_len_with_stop = tx.cds_length + sum(e - s + 1 for s, e in tx.stop_codon)
    if cds_pos is None or cds_len_with_stop <= 0:
        return float("nan")
    return cds_pos / cds_len_with_stop


# ---------------------------------------------------------------------------
# Region builders for conservation queries
# ---------------------------------------------------------------------------

def ptc_to_ejc_interval(
    tx: TranscriptRecord, ptc_genomic_pos: int, loc: PTCLocation,
) -> Optional[Tuple[str, int, int]]:
    """Within-exon downstream interval from PTC to the end of the current exon.

    Matches `conservation_score_extraction_v2.py` exactly. Training does NOT
    check last-exon status — it builds this region whenever PTC is not at
    the exon's 3' end, including PTCs in the last exon. So the name
    "ptc_to_ejc" is a slight misnomer: it's really "PTC to the 3' end of the
    current exon".

    Training uses 1-based start coordinates fed through `add_bed`:
      + strand: [ptc, exon_end]   (PTC base IS included)
      - strand: [exon_start, ptc]

    Returns (chrom, gstart_1based, gend_1based) or None if the interval
    would be empty (PTC exactly at exon 3' end).
    """
    # exons are in transcript order; find the current exon's genomic bounds
    gstart, gend = tx.exons[loc.exon_rank - 1]
    if tx.strand == "+":
        if gend <= ptc_genomic_pos:
            return None
        start = ptc_genomic_pos
        end   = gend
    else:
        if gstart >= ptc_genomic_pos:
            return None
        start = gstart
        end   = ptc_genomic_pos
    if end < start:
        return None
    return (tx.chrom, start, end)


def new3utr_blocks(
    tx: TranscriptRecord, ptc_genomic_pos: int, loc: PTCLocation, window: int = 200,
) -> Tuple[List[Tuple[str, int, int]], int]:
    """Build first N transcript-bases of post-PTC region, as genomic blocks.

    Walks transcript exons from the PTC-containing exon through the 3' end,
    clipping the PTC exon at (ptc±1, exon boundary) and taking the first
    `window` transcript bases total.

    Returns (list of (chrom, gstart_1based, gend_1based), actual_bases_taken).
    """
    blocks: List[Tuple[str, int, int]] = []
    remaining = window

    # index into tx.exons; tx.exons is already transcript-ordered 5'->3'
    for rank, (gstart, gend) in enumerate(tx.exons[loc.exon_rank - 1:], start=loc.exon_rank):
        if remaining <= 0:
            break

        # Clip the first block (PTC exon) to post-PTC portion
        if rank == loc.exon_rank:
            if tx.strand == "+":
                bs = ptc_genomic_pos + 1
                be = gend
            else:
                bs = gstart
                be = ptc_genomic_pos - 1
            if be < bs:
                continue
        else:
            bs, be = gstart, gend

        block_len = be - bs + 1
        if block_len <= 0:
            continue

        if block_len <= remaining:
            blocks.append((tx.chrom, bs, be))
            remaining -= block_len
        else:
            # trim to remaining length on the 5' (transcript) side of this block
            if tx.strand == "+":
                blocks.append((tx.chrom, bs, bs + remaining - 1))
            else:
                blocks.append((tx.chrom, be - remaining + 1, be))
            remaining = 0
            break

    taken = window - remaining
    return blocks, taken
