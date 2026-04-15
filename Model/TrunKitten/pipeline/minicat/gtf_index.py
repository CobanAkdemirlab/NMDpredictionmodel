"""GTF parsing → per-transcript exon/CDS index.

Builds a compact in-memory structure:

    tx_index[tx_id_no_version] = TranscriptRecord(
        transcript_id_raw,     # original ID with version, if any
        transcript_id,         # version-stripped
        gene_id,               # version-stripped ENSG
        chrom, strand,
        exons = [(gstart, gend), ...] sorted in transcript 5'->3' order,
        cds   = [(gstart, gend), ...] sorted in transcript 5'->3' order,
    )

All coordinates are 1-based inclusive (GTF convention).
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import gzip
import re
import logging

log = logging.getLogger(__name__)

VERSION_RE = re.compile(r"\.\d+$")


def strip_version(tid: Optional[str]) -> Optional[str]:
    """Strip '.N' version suffix from Ensembl/Gencode IDs."""
    if tid is None:
        return None
    return VERSION_RE.sub("", str(tid))


@dataclass
class TranscriptRecord:
    transcript_id_raw: str
    transcript_id: str            # version-stripped
    gene_id: str                  # version-stripped
    gene_id_raw: str
    chrom: str
    strand: str                   # '+' or '-'
    exons: List[Tuple[int, int]] = field(default_factory=list)  # sorted tx 5'->3'
    cds:   List[Tuple[int, int]] = field(default_factory=list)  # sorted tx 5'->3'
    # stop_codon records from GTF, sorted tx 5'->3'. Typically exactly one
    # 3-bp range, but split across exons in rare cases (two ranges summing to 3 bp).
    stop_codon: List[Tuple[int, int]] = field(default_factory=list)

    @property
    def exon_count(self) -> int:
        return len(self.exons)

    @property
    def transcript_length(self) -> int:
        return sum(e - s + 1 for s, e in self.exons)

    @property
    def cds_length(self) -> int:
        return sum(e - s + 1 for s, e in self.cds)

    def cds_with_stop(self) -> List[Tuple[int, int]]:
        """Return CDS blocks extended to include the stop codon, in tx order.

        Matches Bioconductor `cdsBy(makeTxDbFromGFF(gencode))` output, which
        auto-appends stop_codon records to CDS ranges. This is the convention
        used by Iman's R pipeline for all *cdsseqs_* / *cdsseq_* sequence
        composition features.
        """
        if not self.stop_codon:
            return list(self.cds)
        # Concatenate in transcript order (cds already ordered; stop_codon already ordered).
        # stop_codon comes after the last CDS segment in transcript orientation.
        return list(self.cds) + list(self.stop_codon)

    def coding_exon_flags(self) -> List[bool]:
        """For each exon in transcript order, True if it overlaps any CDS-with-stop segment.

        Matches Iman's R pipeline convention: coding exons are identified via
        `cdsBy(makeTxDbFromGFF(...))` which extends CDS ranges to include
        stop_codon records. So an exon containing ONLY the stop codon (and
        no sense codons) still counts as coding — matters for transcripts
        where the stop codon sits in its own tiny terminal exon.
        """
        blocks = self.cds_with_stop()
        flags = []
        for (es, ee) in self.exons:
            overlaps = any(not (ee < cs or es > ce) for (cs, ce) in blocks)
            flags.append(overlaps)
        return flags


_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def _parse_attrs(attr_field: str) -> Dict[str, str]:
    return dict(_ATTR_RE.findall(attr_field))


def _open_gtf(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def build_transcript_index(
    gtf_path: Path,
    wanted_tx_ids_stripped: Optional[set[str]] = None,
) -> Dict[str, TranscriptRecord]:
    """Parse GTF and return {version-stripped transcript_id: TranscriptRecord}.

    Args:
        gtf_path: path to .gtf or .gtf.gz (Gencode convention).
        wanted_tx_ids_stripped: optional allowlist for early filtering
            (all version-stripped). Dramatically speeds up annotation when
            the input is a small cohort.

    Returns:
        dict keyed by version-stripped transcript id.
    """
    log.info(f"Parsing GTF: {gtf_path}")
    records: Dict[str, TranscriptRecord] = {}
    n_lines = 0
    n_kept  = 0

    with _open_gtf(gtf_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            n_lines += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _src, feature, start, end, _score, strand, _frame, attrs = parts
            if feature not in ("exon", "CDS", "transcript", "stop_codon"):
                continue

            a = _parse_attrs(attrs)
            tx_raw = a.get("transcript_id")
            if tx_raw is None:
                continue
            tx_stripped = strip_version(tx_raw)
            if wanted_tx_ids_stripped is not None and tx_stripped not in wanted_tx_ids_stripped:
                continue

            gene_raw = a.get("gene_id", "")
            gene_stripped = strip_version(gene_raw) or ""

            rec = records.get(tx_stripped)
            if rec is None:
                rec = TranscriptRecord(
                    transcript_id_raw=tx_raw,
                    transcript_id=tx_stripped,
                    gene_id=gene_stripped,
                    gene_id_raw=gene_raw,
                    chrom=chrom,
                    strand=strand,
                )
                records[tx_stripped] = rec

            if feature == "exon":
                rec.exons.append((int(start), int(end)))
            elif feature == "CDS":
                rec.cds.append((int(start), int(end)))
            elif feature == "stop_codon":
                rec.stop_codon.append((int(start), int(end)))
            # transcript feature is only used to populate gene_id if not seen

    # Sort exons and CDS in transcript (5'->3') order
    for rec in records.values():
        if rec.strand == "+":
            rec.exons.sort(key=lambda t: t[0])
            rec.cds.sort(key=lambda t: t[0])
            rec.stop_codon.sort(key=lambda t: t[0])
        else:
            rec.exons.sort(key=lambda t: -t[0])
            rec.cds.sort(key=lambda t: -t[0])
            rec.stop_codon.sort(key=lambda t: -t[0])
        n_kept += 1

    log.info(f"  lines scanned: {n_lines:,}  transcripts indexed: {n_kept:,}")
    return records


def lookup_transcript(
    tx_index: Dict[str, TranscriptRecord],
    txname: str,
    strip_versions: bool = True,
) -> Tuple[Optional[TranscriptRecord], bool]:
    """Look up a transcript by txname; return (record, version_mismatch_flag).

    Match order: exact → version-stripped (if strip_versions).
    """
    if txname in tx_index:
        return tx_index[txname], False
    stripped = strip_version(txname) if strip_versions else txname
    if stripped and stripped in tx_index:
        return tx_index[stripped], (stripped != txname)
    return None, False
