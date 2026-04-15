"""Toy-transcript validation tests — no external files required.

Run with: pytest -xvs tests/test_toy_transcripts.py

These tests validate the transcript-coordinate math (last.EJC, mut.exon,
relativePTClocation, AmountExonsAfter, region builders) using hand-built
TranscriptRecord objects. Conservation and FASTA paths are mocked.
"""
from __future__ import annotations
import math
import pytest

from minicat.gtf_index import TranscriptRecord
from minicat.transcript import (
    locate_ptc_in_transcript,
    last_ejc_category,
    amount_exons_after,
    relative_ptc_location,
    coding_exon_rank,
    ptc_to_ejc_interval,
    new3utr_blocks,
)


# ---------------------------------------------------------------------------
# Toy 1 — + strand, 3 exons, with partial 5' UTR in exon 1
# ---------------------------------------------------------------------------
#   exon1: 100-199 (len 100)   → includes 5'UTR + start of CDS
#   exon2: 300-399 (len 100)
#   exon3: 500-699 (len 200)  ← last exon; CDS ends at 602 (stop codon 600-602)
#   CDS:   150-199, 300-399, 500-599
#
def _toy1() -> TranscriptRecord:
    return TranscriptRecord(
        transcript_id_raw="ENSTTOY01.1",
        transcript_id="ENSTTOY01",
        gene_id="ENSGTOY01",
        gene_id_raw="ENSGTOY01.1",
        chrom="chr1",
        strand="+",
        exons=[(100, 199), (300, 399), (500, 699)],
        cds=[(150, 199), (300, 399), (500, 599)],
    )


def test_toy1_penultimate_boundary_in_last_50bp():
    """PTC at 350 → penultimate exon, distance-to-3'end = 49 → penultimate.last50bp."""
    tx = _toy1()
    loc = locate_ptc_in_transcript(tx, 350)
    assert loc is not None
    assert loc.exon_rank == 2
    assert loc.exon_count == 3
    # within-exon offset on +strand: 350 - 300 + 1 = 51
    assert loc.within_exon_offset == 51
    # distance to 3' end of this exon = 100 - 51 = 49 → ≤ 50 → penultimate.last50bp
    assert last_ejc_category(tx, loc) == "penultimate.last50bp"
    # mut.exon (coding rank): exon 1 has CDS (150-199) → coding; exon 2 coding → rank 2
    assert coding_exon_rank(tx, loc.exon_rank) == 2
    # relativePTClocation = CDS_pos / CDS_length
    # CDS segments (+strand, transcript order): 150-199 (50), 300-399 (100), 500-599 (100) = 250 total
    # PTC at 350 → cumulative before this segment = 50; within = 350 - 300 + 1 = 51
    # cds_pos = 50 + 51 = 101 ; relative = 101/250 = 0.404
    assert math.isclose(relative_ptc_location(tx, 350), 101 / 250, rel_tol=1e-9)
    # AmountExonsAfter: exon 3 is coding → 1
    assert amount_exons_after(tx, 2) == 1


def test_toy1_penultimate_beyond_50bp():
    """PTC at 340 → distance-to-3'end = 59 → upstream."""
    tx = _toy1()
    loc = locate_ptc_in_transcript(tx, 340)
    assert last_ejc_category(tx, loc) == "upstream"


def test_toy1_last_exon():
    tx = _toy1()
    loc = locate_ptc_in_transcript(tx, 550)
    assert loc.exon_rank == 3
    assert last_ejc_category(tx, loc) == "last.exon"
    assert amount_exons_after(tx, 3) == 0


def test_toy1_first_exon_upstream():
    tx = _toy1()
    loc = locate_ptc_in_transcript(tx, 180)
    assert loc.exon_rank == 1
    assert last_ejc_category(tx, loc) == "upstream"
    assert amount_exons_after(tx, 1) == 2


def test_toy1_ptc_to_ejc_plus_strand():
    tx = _toy1()
    loc = locate_ptc_in_transcript(tx, 350)
    iv = ptc_to_ejc_interval(tx, 350, loc)
    # Training convention (conservation_score_extraction_v2.py): includes PTC base
    assert iv == ("chr1", 350, 399)


def test_toy1_new3utr_blocks_window():
    tx = _toy1()
    loc = locate_ptc_in_transcript(tx, 350)
    blocks, taken = new3utr_blocks(tx, 350, loc, window=200)
    # post-PTC in exon 2: 351-399 (49 bp) + exon 3: 500-699 (200 bp max, but cap = 151)
    assert blocks[0] == ("chr1", 351, 399)
    assert blocks[1] == ("chr1", 500, 500 + 151 - 1)
    assert taken == 200


# ---------------------------------------------------------------------------
# Toy 2 — − strand, 3 exons (mirror of Toy1 geometry)
# ---------------------------------------------------------------------------
# Same genomic coordinates but strand = '-'.
# On '-' strand the TRANSCRIPT-ORDER is exon with largest genomic start first.
# So exon1 (tx-order) = genomic 500-699; exon2 = 300-399; exon3 = 100-199 (last).
def _toy2() -> TranscriptRecord:
    return TranscriptRecord(
        transcript_id_raw="ENSTTOY02.1",
        transcript_id="ENSTTOY02",
        gene_id="ENSGTOY02",
        gene_id_raw="ENSGTOY02.1",
        chrom="chr2",
        strand="-",
        # after gtf_index sorts: transcript 5'→3' = descending genomic start
        exons=[(500, 699), (300, 399), (100, 199)],
        cds=[(500, 599), (300, 399), (150, 199)],
    )


def test_toy2_minus_strand_last_exon():
    """PTC at genomic 150 → transcript-last exon → last.exon."""
    tx = _toy2()
    loc = locate_ptc_in_transcript(tx, 150)
    # exon_rank 3 is the genomically-leftmost exon = last in transcript order
    assert loc.exon_rank == 3
    assert loc.exon_count == 3
    assert last_ejc_category(tx, loc) == "last.exon"
    assert amount_exons_after(tx, 3) == 0


def test_toy2_minus_strand_within_exon_offset():
    tx = _toy2()
    # PTC at 350 → in exon 2 (tx-order: 300-399 on minus strand)
    # within-exon offset on -strand: gend - ptc + 1 = 399 - 350 + 1 = 50
    loc = locate_ptc_in_transcript(tx, 350)
    assert loc.exon_rank == 2
    assert loc.within_exon_offset == 50
    # distance to exon 3' end on -strand = exon_length - within_offset = 100 - 50 = 50 → ≤ 50
    assert last_ejc_category(tx, loc) == "penultimate.last50bp"


def test_toy2_minus_strand_ptc_to_ejc_interval():
    tx = _toy2()
    loc = locate_ptc_in_transcript(tx, 350)
    iv = ptc_to_ejc_interval(tx, 350, loc)
    # -strand, PTC-containing exon = (300, 399); training includes PTC base: [exon_start, ptc]
    assert iv == ("chr2", 300, 350)


def test_toy2_minus_strand_relative_ptc_location():
    tx = _toy2()
    # CDS on -strand (transcript order): 500-599 (100), 300-399 (100), 150-199 (50) = 250
    # PTC at 350 in second CDS segment (300-399) on -strand
    # cumulative before = 100; within = gend - ptc + 1 = 399 - 350 + 1 = 50
    # cds_pos = 100 + 50 = 150 ; relative = 150/250 = 0.6
    assert math.isclose(relative_ptc_location(tx, 350), 150 / 250, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# Toy 3 — single-exon transcript
# ---------------------------------------------------------------------------
def _toy3() -> TranscriptRecord:
    return TranscriptRecord(
        transcript_id_raw="ENSTTOY03.1",
        transcript_id="ENSTTOY03",
        gene_id="ENSGTOY03",
        gene_id_raw="ENSGTOY03.1",
        chrom="chr3",
        strand="+",
        exons=[(1000, 1999)],
        cds=[(1100, 1899)],
    )


def test_toy3_single_exon_always_last_exon():
    tx = _toy3()
    loc = locate_ptc_in_transcript(tx, 1500)
    assert loc.exon_rank == 1
    assert loc.exon_count == 1
    assert last_ejc_category(tx, loc) == "last.exon"
    assert amount_exons_after(tx, 1) == 0
    # ptc_to_ejc: training computes this even in last exon (PTC to exon 3' end)
    # exon is 1000-1999 on +strand; PTC at 1500 → interval [1500, 1999]
    iv = ptc_to_ejc_interval(tx, 1500, loc)
    assert iv == ("chr3", 1500, 1999)


def test_toy3_ptc_at_exon_3prime_end_returns_none():
    """PTC at exact 3' end of exon → no downstream sequence → None."""
    tx = _toy3()
    loc = locate_ptc_in_transcript(tx, 1999)  # last base of exon
    iv = ptc_to_ejc_interval(tx, 1999, loc)
    assert iv is None


# ---------------------------------------------------------------------------
# Toy 4 — + strand, 5'UTR-only first exon → validates mut.exon coding-rank
# ---------------------------------------------------------------------------
# exon 1: 100-199 (5'UTR only — no CDS overlap)
# exon 2: 300-399 (CDS starts here: 350-399)
# exon 3: 500-699 (last exon; CDS 500-599; stop 600-602)
def _toy4() -> TranscriptRecord:
    return TranscriptRecord(
        transcript_id_raw="ENSTTOY04.1",
        transcript_id="ENSTTOY04",
        gene_id="ENSGTOY04",
        gene_id_raw="ENSGTOY04.1",
        chrom="chr4",
        strand="+",
        exons=[(100, 199), (300, 399), (500, 699)],
        cds=[(350, 399), (500, 599)],
    )


def test_toy4_mut_exon_coding_rank_differs_from_transcript_rank():
    """Exon 1 is UTR-only; CDS starts in exon 2. Training uses coding rank.

    For a PTC at 370 (exon 2), transcript exon_rank=2 but coding rank=1,
    since exon 1 has no CDS overlap.
    """
    tx = _toy4()
    loc = locate_ptc_in_transcript(tx, 370)
    assert loc.exon_rank == 2
    # Coding rank: exon 1 non-coding, exon 2 coding → coding rank 1
    assert coding_exon_rank(tx, loc.exon_rank) == 1


def test_toy4_relative_ptc_cds_based():
    """Verify CDS-based relativePTClocation ignores 5' UTR.

    CDS: 350-399 (50), 500-599 (100) = 150 total
    PTC at 370 → CDS pos = 370 - 350 + 1 = 21
    relative = 21/150 = 0.14
    """
    tx = _toy4()
    assert math.isclose(relative_ptc_location(tx, 370), 21 / 150, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# Edge case — boundary PTC
# ---------------------------------------------------------------------------
def test_boundary_position_flag():
    tx = _toy1()
    loc = locate_ptc_in_transcript(tx, 199)  # exon 1 end
    assert loc.exon_rank == 1
    assert loc.boundary_ambiguous is True


# ---------------------------------------------------------------------------
# Toy 6 — + strand with native stop mid-penultimate-exon, trailing 3'UTR
# ---------------------------------------------------------------------------
# This is the structural case that caused the 10 last.EJC mismatches at n=6019.
#
# exon 1: 100-199 (CDS 150-199)
# exon 2: 300-499 (CDS 300-399, stop 400-402, rest 403-499 is 3'UTR within exon 2)
# exon 3: 600-699 (pure 3'UTR — no CDS, no stop_codon)
#
# Under Iman's cdsBy() convention, exon 3 is NOT coding (no CDS, no stop_codon).
# So the PENULTIMATE coding exon is exon 1, LAST coding exon is exon 2.
# For a PTC in exon 2: rank 2 of 2 coding exons → last.exon.
# For a PTC in exon 1: rank 1 of 2 → penultimate.
#
# BUT — the older convention that measured distance to the exon's 3' end
# (not CDS end) would give wrong answers if the penultimate had trailing
# 3'UTR. Here the penultimate is exon 1 which has no trailing UTR, so this
# test validates the correct 2-coding-exon topology.
def _toy6() -> TranscriptRecord:
    return TranscriptRecord(
        transcript_id_raw="ENSTTOY06.1",
        transcript_id="ENSTTOY06",
        gene_id="ENSGTOY06",
        gene_id_raw="ENSGTOY06.1",
        chrom="chr6",
        strand="+",
        exons=[(100, 199), (300, 499), (600, 699)],
        cds=[(150, 199), (300, 399)],    # CDS ends at 399, mid-exon-2
        stop_codon=[(400, 402)],          # stop codon within exon 2
    )


def test_toy6_trailing_utr_exon_not_coding():
    """Pure-UTR exon 3 is not coding; only exons 1 and 2 are coding."""
    tx = _toy6()
    assert tx.coding_exon_flags() == [True, True, False]


def test_toy6_ptc_in_exon2_is_last_exon():
    """Exon 2 is the LAST coding exon (exon 3 is pure UTR)."""
    tx = _toy6()
    loc = locate_ptc_in_transcript(tx, 350)
    assert coding_exon_rank(tx, loc.exon_rank) == 2
    assert last_ejc_category(tx, loc) == "last.exon"


def test_toy6_ptc_in_exon1_is_penultimate():
    """PTC in exon 1 at position 180: c_rank=1, n_coding=2 → penultimate.
    exon 1 is all CDS (no trailing UTR), within_offset=81, exon_len=100,
    distance_to_3p = 19 ≤ 50, penultimate CDS length = 100 ≥ 50 → penultimate.last50bp.
    """
    tx = _toy6()
    loc = locate_ptc_in_transcript(tx, 180)
    assert last_ejc_category(tx, loc) == "penultimate.last50bp"


# ---------------------------------------------------------------------------
# Toy 7 — + strand with SHORT penultimate coding exon (< 50 bp)
# ---------------------------------------------------------------------------
# Matches Iman's `penultimate_length >= 50` guard in Variant_annotation.R:
# if the penultimate coding exon is itself < 50 bp, a PTC inside it can't
# be classified as penultimate.last50bp — it's `upstream` instead.
#
# This is the structural case that caused the final 10 last.EJC mismatches
# at n=6019 (immunoglobulin/TCR gene segments with tiny exons).
#
# exon 1: 100-199  (CDS 150-199, length 50)
# exon 2: 300-329  (CDS 300-329, length 30 — SHORT penultimate)
# exon 3: 500-699  (CDS 500-599, stop 600-602 — last coding exon)
def _toy7() -> TranscriptRecord:
    return TranscriptRecord(
        transcript_id_raw="ENSTTOY07.1",
        transcript_id="ENSTTOY07",
        gene_id="ENSGTOY07",
        gene_id_raw="ENSGTOY07.1",
        chrom="chr7",
        strand="+",
        exons=[(100, 199), (300, 329), (500, 699)],
        cds=[(150, 199), (300, 329), (500, 599)],
        stop_codon=[(600, 602)],
    )


def test_toy7_short_penultimate_exon_is_upstream():
    """PTC inside a 30-bp penultimate coding exon → upstream, not penultimate.last50bp.

    Penultimate CDS-exon length = 30 < 50, so the length guard kicks in.
    Without the guard my code would call this penultimate.last50bp because
    the PTC is trivially within 50 bp of the exon's 3' end.
    """
    tx = _toy7()
    loc = locate_ptc_in_transcript(tx, 325)   # 4 bp from exon 2's 3' end
    assert last_ejc_category(tx, loc) == "upstream"


def test_toy7_exactly_50bp_penultimate_qualifies():
    """Penultimate exon with CDS-length == 50 bp passes the >= 50 guard."""
    # Modify toy7 so exon 2 is 50 bp exactly
    tx = TranscriptRecord(
        transcript_id_raw="ENSTTOY07b.1",
        transcript_id="ENSTTOY07b",
        gene_id="ENSGTOY07b",
        gene_id_raw="ENSGTOY07b.1",
        chrom="chr7",
        strand="+",
        exons=[(100, 199), (300, 349), (500, 699)],
        cds=[(150, 199), (300, 349), (500, 599)],
        stop_codon=[(600, 602)],
    )
    loc = locate_ptc_in_transcript(tx, 340)   # 9 bp from 3' end, within 50
    assert last_ejc_category(tx, loc) == "penultimate.last50bp"


def test_position_outside_all_exons():
    tx = _toy1()
    assert locate_ptc_in_transcript(tx, 250) is None  # intronic
    assert locate_ptc_in_transcript(tx, 50) is None   # upstream of tx


# ---------------------------------------------------------------------------
# Toy 5 — + strand with stop_codon in its own terminal exon
# ---------------------------------------------------------------------------
# This is the structural case that caused the 6019-variant AmountExonsAfter
# and last.EJC mismatches. exon 3 contains ONLY the stop codon + 3'UTR, and
# the raw GTF `CDS` records do NOT overlap it — but Bioconductor's cdsBy()
# extends CDS to include stop_codon, making exon 3 count as coding.
#
# exon 1: 100-199 (CDS 150-199)
# exon 2: 300-399 (CDS 300-399)
# exon 3: 500-699 (stop_codon 500-502 only; rest is 3'UTR)
def _toy5() -> TranscriptRecord:
    return TranscriptRecord(
        transcript_id_raw="ENSTTOY05.1",
        transcript_id="ENSTTOY05",
        gene_id="ENSGTOY05",
        gene_id_raw="ENSGTOY05.1",
        chrom="chr5",
        strand="+",
        exons=[(100, 199), (300, 399), (500, 699)],
        cds=[(150, 199), (300, 399)],          # no CDS in exon 3
        stop_codon=[(500, 502)],                # stop in exon 3
    )


def test_toy5_stop_codon_exon_counts_as_coding():
    """Exon containing only the stop codon should count as coding (cdsBy convention)."""
    tx = _toy5()
    flags = tx.coding_exon_flags()
    assert flags == [True, True, True]  # all 3 coding via cds_with_stop()


def test_toy5_amount_exons_after_counts_stop_codon_exon():
    """A PTC in exon 2 has 1 coding exon after (the stop-codon-only exon 3)."""
    tx = _toy5()
    loc = locate_ptc_in_transcript(tx, 350)
    assert amount_exons_after(tx, loc.exon_rank) == 1


def test_toy5_last_ejc_uses_coding_rank():
    """PTC in exon 2 (penultimate coding). Distance to exon 3' end = 49 → penultimate.last50bp."""
    tx = _toy5()
    # PTC at 350: within_offset = 51, exon_length = 100, distance_to_3p = 49 ≤ 50
    loc = locate_ptc_in_transcript(tx, 350)
    assert last_ejc_category(tx, loc) == "penultimate.last50bp"

    # PTC at 340: distance_to_3p = 59 > 50
    loc = locate_ptc_in_transcript(tx, 340)
    assert last_ejc_category(tx, loc) == "upstream"


def test_toy5_ptc_in_stop_codon_exon_is_last_exon():
    """A PTC in the stop-codon-only exon → last.exon (no coding exons after)."""
    tx = _toy5()
    # Use a position inside exon 3 but within the stop codon range
    loc = locate_ptc_in_transcript(tx, 501)
    assert last_ejc_category(tx, loc) == "last.exon"
    assert amount_exons_after(tx, loc.exon_rank) == 0
