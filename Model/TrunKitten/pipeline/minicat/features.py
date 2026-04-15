"""Per-variant feature orchestration.

Combines GTF index, FASTA, conservation sources, and half-life table into
a single callable that takes a variant row and returns a feature dict.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Any, Dict, Optional
import logging
import math

from pyfaidx import Fasta

from .gtf_index import TranscriptRecord, lookup_transcript, strip_version
from .transcript import (
    PTCLocation,
    locate_ptc_in_transcript,
    amount_exons_after,
    last_ejc_category,
    relative_ptc_location,
    coding_exon_rank,
    ptc_to_ejc_interval,
    new3utr_blocks,
)
from .sequence import compute_cds_composition
from .conservation import ConservationSource
from .halflife import HalfLifeTable

log = logging.getLogger(__name__)


@dataclass
class AnnotationResult:
    # identifiers / QC
    variant_id: str
    txnames: str
    transcript_id_used: Optional[str]
    gene: str
    gene_id: Optional[str]
    strand: Optional[str]
    exon_count: Optional[int]
    coding_exon_count: Optional[int]
    ptc_transcript_pos: Optional[int]
    transcript_length: Optional[int]
    cds_length: Optional[int]
    downstream_new3utr_len: Optional[int]
    boundary_ambiguous: bool
    tx_not_in_gtf: bool
    version_mismatch: bool
    cds_length_short_flag: bool
    new3utr_empty: bool
    ptc_to_ejc_empty: bool
    half_life_missing: bool
    any_conservation_missing: bool
    # features
    last_EJC: Optional[str]
    relativePTClocation: Optional[float]
    half_life_PC1: Optional[float]
    cdsseqs_AU_content: Optional[float]
    mut_exon: Optional[int]
    phastcons_new3utr_first200_median: Optional[float]
    phylop_ptc_to_ejc_median: Optional[float]
    AmountExonsAfter: Optional[int]
    cdsseq_AUcontentlast200: Optional[float]
    cdsseqs_UC_content: Optional[float]

    def to_feature_row(self) -> Dict[str, Any]:
        """10-feature row with the canonical column names (including dots)."""
        return {
            "variant_id": self.variant_id,
            "txnames": self.txnames,
            "transcript_id_used": self.transcript_id_used,
            "gene": self.gene,
            "gene_id": self.gene_id,
            "strand": self.strand,
            "last.EJC": self.last_EJC,
            "relativePTClocation": self.relativePTClocation,
            "half_life_PC1": self.half_life_PC1,
            "cdsseqs_AU_content": self.cdsseqs_AU_content,
            "mut.exon": self.mut_exon,
            "phastcons_new3utr_first200_median": self.phastcons_new3utr_first200_median,
            "phylop_ptc_to_ejc_median": self.phylop_ptc_to_ejc_median,
            "AmountExonsAfter": self.AmountExonsAfter,
            "cdsseq_AUcontentlast200": self.cdsseq_AUcontentlast200,
            "cdsseqs_UC_content": self.cdsseqs_UC_content,
        }

    def to_qc_row(self) -> Dict[str, Any]:
        return {
            "variant_id": self.variant_id,
            "txnames": self.txnames,
            "transcript_id_used": self.transcript_id_used,
            "gene_id": self.gene_id,
            "strand": self.strand,
            "exon_count": self.exon_count,
            "coding_exon_count": self.coding_exon_count,
            "ptc_transcript_pos": self.ptc_transcript_pos,
            "transcript_length": self.transcript_length,
            "cds_length": self.cds_length,
            "downstream_new3utr_len": self.downstream_new3utr_len,
            "boundary_ambiguous": self.boundary_ambiguous,
            "tx_not_in_gtf": self.tx_not_in_gtf,
            "version_mismatch": self.version_mismatch,
            "cds_length_short_flag": self.cds_length_short_flag,
            "new3utr_empty": self.new3utr_empty,
            "ptc_to_ejc_empty": self.ptc_to_ejc_empty,
            "half_life_missing": self.half_life_missing,
            "any_conservation_missing": self.any_conservation_missing,
        }


class FeatureAnnotator:
    """Stateful annotator that holds all the heavy resources."""

    def __init__(
        self,
        tx_index: Dict[str, TranscriptRecord],
        fasta: Fasta,
        phastcons: ConservationSource,
        phylop: ConservationSource,
        halflife: HalfLifeTable,
        strip_versions: bool = True,
        cds_last_window: int = 200,
        new3utr_window: int = 200,
    ):
        self.tx_index = tx_index
        self.fasta = fasta
        self.phastcons = phastcons
        self.phylop = phylop
        self.halflife = halflife
        self.strip_versions = strip_versions
        self.cds_last_window = cds_last_window
        self.new3utr_window = new3utr_window

    def annotate(self, row: Dict[str, Any]) -> AnnotationResult:
        variant_id = str(row["variant_id"])
        txname = str(row["txnames"]).strip()
        gene = str(row.get("gene", "")) if row.get("gene") is not None else ""
        contig = str(row["contig"])
        ptc_pos = int(row["position"])

        tx, version_mismatch = lookup_transcript(
            self.tx_index, txname, strip_versions=self.strip_versions,
        )

        if tx is None:
            log.warning(f"  {variant_id}: transcript '{txname}' not found in GTF")
            return self._nan_result(variant_id, txname, gene,
                                    tx_not_in_gtf=True, version_mismatch=False)

        # Harmonise contig vs GTF chrom convention
        if tx.chrom != contig:
            # try toggling "chr" prefix
            alt = contig[3:] if contig.startswith("chr") else f"chr{contig}"
            if tx.chrom != alt:
                log.warning(
                    f"  {variant_id}: chromosome mismatch variant={contig} "
                    f"transcript={tx.chrom}"
                )

        # --- Locate PTC in transcript ---
        loc = locate_ptc_in_transcript(tx, ptc_pos)
        if loc is None:
            log.warning(
                f"  {variant_id}: position {ptc_pos} not within any exon of "
                f"transcript {tx.transcript_id} (strand={tx.strand})"
            )
            return self._nan_result(variant_id, txname, gene,
                                    tx_not_in_gtf=False,
                                    version_mismatch=version_mismatch,
                                    transcript_id_used=tx.transcript_id,
                                    gene_id=tx.gene_id,
                                    strand=tx.strand)

        # --- Categorical + positional features ---
        # NOTE on training conventions (verified against TOPMed_merged_v4.csv):
        #   - last.EJC          → all transcript exons (matches 100%)
        #   - mut.exon          → coding-exon rank (falls back to transcript rank
        #                         if PTC exon is non-coding — shouldn't happen for stop-gains)
        #   - relativePTClocation → PTC_CDS_pos / CDS_length (CDS-internal; NOT tx-spliced)
        #   - AmountExonsAfter  → coding exons strictly after (matches 100%)
        last_ejc = last_ejc_category(tx, loc)
        rel_ptc  = relative_ptc_location(tx, ptc_pos)
        coding_rank = coding_exon_rank(tx, loc.exon_rank)
        mut_exon = coding_rank if coding_rank is not None else loc.exon_rank
        n_after  = amount_exons_after(tx, loc.exon_rank)
        coding_exon_count = sum(tx.coding_exon_flags())

        # --- CDS sequence composition features ---
        cds_comp = compute_cds_composition(
            self.fasta, tx, last_window=self.cds_last_window,
        )

        # --- half_life_PC1 ---
        # Primary key is GTF-derived ENSG; gene-symbol fallback handles cases
        # where the transcript's current gene_id differs from the ENSG used to
        # build the half-life table (annotation-version drift).
        hl = self.halflife.lookup(
            tx.gene_id,
            strip_versions=self.strip_versions,
            gene_symbol=gene if gene else None,
        )
        half_life_missing = hl is None

        # --- phastcons_new3utr_first200_median ---
        new3_blocks, taken = new3utr_blocks(
            tx, ptc_pos, loc, window=self.new3utr_window,
        )
        new3utr_empty = (len(new3_blocks) == 0)
        if new3utr_empty:
            phc_new3 = float("nan")
        else:
            phc_new3, _bp, valid = self.phastcons.median_over_blocks(new3_blocks)

        # --- phylop_ptc_to_ejc_median ---
        p2e = ptc_to_ejc_interval(tx, ptc_pos, loc)
        ptc_to_ejc_empty = (p2e is None)
        if ptc_to_ejc_empty:
            phy_p2e = float("nan")
        else:
            phy_p2e, _bp, _valid = self.phylop.median_over_blocks([p2e])

        any_cons_miss = _is_nan(phc_new3) or _is_nan(phy_p2e)

        return AnnotationResult(
            variant_id=variant_id,
            txnames=txname,
            transcript_id_used=tx.transcript_id,
            gene=gene,
            gene_id=tx.gene_id,
            strand=tx.strand,
            exon_count=tx.exon_count,
            coding_exon_count=coding_exon_count,
            ptc_transcript_pos=loc.transcript_pos,
            transcript_length=loc.transcript_length,
            cds_length=cds_comp["cds_length"],
            downstream_new3utr_len=taken,
            boundary_ambiguous=loc.boundary_ambiguous,
            tx_not_in_gtf=False,
            version_mismatch=version_mismatch,
            cds_length_short_flag=cds_comp["cds_length_short_flag"],
            new3utr_empty=new3utr_empty,
            ptc_to_ejc_empty=ptc_to_ejc_empty,
            half_life_missing=half_life_missing,
            any_conservation_missing=any_cons_miss,
            last_EJC=last_ejc,
            relativePTClocation=rel_ptc,
            half_life_PC1=hl,
            cdsseqs_AU_content=cds_comp["cdsseqs_AU_content"],
            mut_exon=mut_exon,
            phastcons_new3utr_first200_median=phc_new3,
            phylop_ptc_to_ejc_median=phy_p2e,
            AmountExonsAfter=n_after,
            cdsseq_AUcontentlast200=cds_comp["cdsseq_AUcontentlast200"],
            cdsseqs_UC_content=cds_comp["cdsseqs_UC_content"],
        )

    def _nan_result(self, variant_id, txname, gene, **kw) -> AnnotationResult:
        base = dict(
            variant_id=variant_id, txnames=txname, gene=gene,
            transcript_id_used=None, gene_id=None, strand=None,
            exon_count=None, coding_exon_count=None,
            ptc_transcript_pos=None, transcript_length=None, cds_length=None,
            downstream_new3utr_len=None,
            boundary_ambiguous=False, tx_not_in_gtf=False, version_mismatch=False,
            cds_length_short_flag=False, new3utr_empty=False, ptc_to_ejc_empty=False,
            half_life_missing=True, any_conservation_missing=True,
            last_EJC=None, relativePTClocation=None, half_life_PC1=None,
            cdsseqs_AU_content=None, mut_exon=None,
            phastcons_new3utr_first200_median=None,
            phylop_ptc_to_ejc_median=None,
            AmountExonsAfter=None, cdsseq_AUcontentlast200=None,
            cdsseqs_UC_content=None,
        )
        base.update(kw)
        return AnnotationResult(**base)


def _is_nan(x) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x))
