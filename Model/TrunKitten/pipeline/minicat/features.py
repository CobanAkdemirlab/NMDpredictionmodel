"""Per-variant feature orchestration.

Combines GTF index, FASTA, phastCons, half-life and GTEx expression tables into
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
    last_ejc_category,
    relative_ptc_location,
    coding_exon_rank,
    new3utr_blocks,
)
from .sequence import compute_cds_composition
from .conservation import ConservationSource
from .halflife import HalfLifeTable
from .expression import GTExExpressionTable

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
    new3utr_empty: bool
    half_life_missing: bool
    expression_missing: bool
    any_conservation_missing: bool
    # features (8 — TrunKitten final feature set; order = trunkitten_features.json)
    last_EJC: Optional[str]
    relativePTClocation: Optional[float]
    half_life_PC1: Optional[float]
    cdsseqs_AU_content: Optional[float]
    mut_exon: Optional[int]
    cdsseqs_UC_content: Optional[float]
    phastcons_new3utr_first200_median: Optional[float]
    MedianExpression_log2: Optional[float]

    def to_feature_row(self) -> Dict[str, Any]:
        """8-feature row (plus identifiers) with the canonical column names (including dots)."""
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
            "cdsseqs_UC_content": self.cdsseqs_UC_content,
            "phastcons_new3utr_first200_median": self.phastcons_new3utr_first200_median,
            "MedianExpression_log2": self.MedianExpression_log2,
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
            "new3utr_empty": self.new3utr_empty,
            "half_life_missing": self.half_life_missing,
            "expression_missing": self.expression_missing,
            "any_conservation_missing": self.any_conservation_missing,
        }


class FeatureAnnotator:
    """Stateful annotator that holds all the heavy resources."""

    def __init__(
        self,
        tx_index: Dict[str, TranscriptRecord],
        fasta: Fasta,
        phastcons: ConservationSource,
        halflife: HalfLifeTable,
        expression: GTExExpressionTable,
        strip_versions: bool = True,
        new3utr_window: int = 200,
    ):
        self.tx_index = tx_index
        self.fasta = fasta
        self.phastcons = phastcons
        self.halflife = halflife
        self.expression = expression
        self.strip_versions = strip_versions
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
        last_ejc = last_ejc_category(tx, loc)
        rel_ptc  = relative_ptc_location(tx, ptc_pos)
        coding_rank = coding_exon_rank(tx, loc.exon_rank)
        mut_exon = coding_rank if coding_rank is not None else loc.exon_rank
        coding_exon_count = sum(tx.coding_exon_flags())

        # --- CDS sequence composition features (AU and UC content) ---
        cds_comp = compute_cds_composition(self.fasta, tx)

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

        # --- MedianExpression_log2 (GTEx v8, same ID/symbol lookup as half-life) ---
        expr = self.expression.lookup(
            tx.gene_id,
            strip_versions=self.strip_versions,
            gene_symbol=gene if gene else None,
        )
        expression_missing = expr is None

        # --- phastcons_new3utr_first200_median ---
        new3_blocks, taken = new3utr_blocks(
            tx, ptc_pos, loc, window=self.new3utr_window,
        )
        new3utr_empty = (len(new3_blocks) == 0)
        if new3utr_empty:
            phc_new3 = float("nan")
        else:
            phc_new3, _bp, valid = self.phastcons.median_over_blocks(new3_blocks)

        any_cons_miss = _is_nan(phc_new3)

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
            new3utr_empty=new3utr_empty,
            half_life_missing=half_life_missing,
            expression_missing=expression_missing,
            any_conservation_missing=any_cons_miss,
            last_EJC=last_ejc,
            relativePTClocation=rel_ptc,
            half_life_PC1=hl,
            cdsseqs_AU_content=cds_comp["cdsseqs_AU_content"],
            mut_exon=mut_exon,
            cdsseqs_UC_content=cds_comp["cdsseqs_UC_content"],
            phastcons_new3utr_first200_median=phc_new3,
            MedianExpression_log2=expr,
        )

    def _nan_result(self, variant_id, txname, gene, **kw) -> AnnotationResult:
        base = dict(
            variant_id=variant_id, txnames=txname, gene=gene,
            transcript_id_used=None, gene_id=None, strand=None,
            exon_count=None, coding_exon_count=None,
            ptc_transcript_pos=None, transcript_length=None, cds_length=None,
            downstream_new3utr_len=None,
            boundary_ambiguous=False, tx_not_in_gtf=False, version_mismatch=False,
            new3utr_empty=False,
            half_life_missing=True, expression_missing=True,
            any_conservation_missing=True,
            last_EJC=None, relativePTClocation=None, half_life_PC1=None,
            cdsseqs_AU_content=None, mut_exon=None,
            cdsseqs_UC_content=None,
            phastcons_new3utr_first200_median=None,
            MedianExpression_log2=None,
        )
        base.update(kw)
        return AnnotationResult(**base)


def _is_nan(x) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x))