"""Contract + unit tests for TrunKitten's 8-feature set — no external data needed.

Run with: pytest -v tests/test_feature_contract.py

Covers:
  * REQUIRED_FEATURES matches the trained model's trunkitten_features.json (same names, same order)
  * AU / UC content (incl. minus strand and N handling)
  * GTEx MedianExpression_log2 (median across tissues, log2(x+1), versions, symbol fallback)
  * FeatureAnnotator output has exactly the 8 features and no retired ones
  * predict.py imputation (median-impute expression/half-life, no zero-fill)
"""
from __future__ import annotations
import gzip
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from minicat import REQUIRED_FEATURES
from minicat.expression import GTExExpressionTable
from minicat.features import FeatureAnnotator
from minicat.gtf_index import TranscriptRecord
from minicat.sequence import _au_content, _uc_content, compute_cds_composition

TRUNKITTEN_DIR = Path(__file__).resolve().parents[2]
RETIRED = {"phylop_ptc_to_ejc_median", "AmountExonsAfter", "cdsseq_AUcontentlast200"}


# ---------------------------------------------------------------------------
# Contract with the trained model
# ---------------------------------------------------------------------------
def _features_json() -> Path:
    override = os.environ.get("TRUNKITTEN_FEATURES_JSON")
    return Path(override) if override else TRUNKITTEN_DIR / "model" / "trunkitten_features.json"


def test_required_features_match_trained_model():
    path = _features_json()
    if not path.exists():
        pytest.skip(f"{path} not found")
    meta = json.loads(path.read_text())
    assert list(REQUIRED_FEATURES) == meta["features_in_order"], (
        "minicat.REQUIRED_FEATURES is out of sync with the trained model's feature list/order. "
        "Column order changes CatBoost results, so both must match exactly."
    )
    assert meta["n_features"] == len(REQUIRED_FEATURES) == 8


def test_no_retired_features_required():
    assert not (set(REQUIRED_FEATURES) & RETIRED)


# ---------------------------------------------------------------------------
# AU / UC content
# ---------------------------------------------------------------------------
def test_au_uc_content_basic():
    assert _au_content("AATTCCGG") == pytest.approx(0.5)
    assert _uc_content("AATTCCGG") == pytest.approx(0.5)
    assert _au_content("AAAA") == 1.0 and _uc_content("AAAA") == 0.0
    assert _uc_content("TTCC") == 1.0 and _au_content("TTCC") == 0.5


def test_au_uc_denominator_includes_n_and_empty_is_nan():
    # Matches R Biostrings::alphabetFrequency(as.prob=TRUE): N stays in the denominator
    assert _au_content("ANNN") == pytest.approx(0.25)
    assert _uc_content("TNNN") == pytest.approx(0.25)
    assert math.isnan(_au_content("")) and math.isnan(_uc_content(""))


def _write_fasta(tmp_path, seq: str, name="chr1") -> Path:
    p = tmp_path / "toy.fa"
    p.write_text(f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    return p


def _toy_genome(n=800) -> str:
    # deterministic, non-trivial composition
    rng = np.random.default_rng(7)
    return "".join(rng.choice(list("ACGT"), size=n, p=[0.35, 0.2, 0.15, 0.3]))


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _expected_cds(genome: str, blocks, strand: str) -> str:
    parts = [genome[s - 1:e] for s, e in blocks]
    if strand == "-":
        parts = [_rc(p) for p in parts]
    return "".join(parts)


def _toy_tx(strand: str) -> TranscriptRecord:
    if strand == "+":
        return TranscriptRecord("ENSTTOY01.1", "ENSTTOY01", "ENSGTOY01", "ENSGTOY01.1", "chr1", "+",
                                exons=[(100, 199), (300, 399), (500, 699)],
                                cds=[(150, 199), (300, 399), (500, 599)],
                                stop_codon=[(600, 602)])
    # minus strand: transcript order is descending genomic
    return TranscriptRecord("ENSTTOY02.1", "ENSTTOY02", "ENSGTOY02", "ENSGTOY02.1", "chr1", "-",
                            exons=[(500, 699), (300, 399), (100, 199)],
                            cds=[(500, 599), (300, 399), (150, 199)],
                            stop_codon=[(147, 149)])


@pytest.mark.parametrize("strand", ["+", "-"])
def test_compute_cds_composition_matches_hand_computation(tmp_path, strand):
    pyfaidx = pytest.importorskip("pyfaidx")
    genome = _toy_genome()
    fa = pyfaidx.Fasta(str(_write_fasta(tmp_path, genome)), sequence_always_upper=True)
    tx = _toy_tx(strand)
    comp = compute_cds_composition(fa, tx)
    cds = _expected_cds(genome, tx.cds_with_stop(), strand)
    assert comp["cds_length"] == len(cds)
    assert comp["cdsseqs_AU_content"] == pytest.approx(_au_content(cds))
    assert comp["cdsseqs_UC_content"] == pytest.approx(_uc_content(cds))
    assert set(comp) == {"cdsseqs_AU_content", "cdsseqs_UC_content", "cds_length"}


# ---------------------------------------------------------------------------
# GTEx MedianExpression_log2
# ---------------------------------------------------------------------------
TISSUES = [f"Tissue_{i}" for i in range(12)]


def _write_gct(path: Path, rows, gz=False, gct_header=True):
    """rows: list of (Name, Description, [12 tissue values])"""
    lines = []
    if gct_header:
        lines += ["#1.2", f"{len(rows)}\t{len(TISSUES)}"]
    lines.append("\t".join(["Name", "Description"] + TISSUES))
    for name, desc, vals in rows:
        lines.append("\t".join([name, desc] + ["" if v is None else str(v) for v in vals]))
    text = "\n".join(lines) + "\n"
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)


ROWS = [
    ("ENSG00000000001.5", "GENEA", [0.0] * 12),                       # median 0   -> log2(1) = 0
    ("ENSG00000000002.7", "GENEB", list(range(1, 13))),               # median 6.5 -> log2(7.5)
    ("ENSG00000000003.1", "GENEC", [10.0] * 6 + [None] * 6),          # NaN ignored -> median 10
    ("ENSG00000000004.2", "GENED", [None] * 12),                      # all NaN -> missing
]


@pytest.mark.parametrize("gz,gct_header", [(True, True), (False, True), (False, False)])
def test_expression_values_match_r_logic(tmp_path, gz, gct_header):
    p = tmp_path / ("expr.gct.gz" if gz else "expr.gct")
    _write_gct(p, ROWS, gz=gz, gct_header=gct_header)
    t = GTExExpressionTable(p)
    assert t.lookup("ENSG00000000001.9") == pytest.approx(0.0)                 # version-stripped
    assert t.lookup("ENSG00000000002") == pytest.approx(math.log2(6.5 + 1))
    assert t.lookup("ENSG00000000003") == pytest.approx(math.log2(10.0 + 1))
    assert t.lookup("ENSG00000000004") is None                                  # all NaN
    assert t.lookup("ENSG99999999999") is None                                  # unknown gene


def test_expression_symbol_fallback_and_switch(tmp_path):
    p = tmp_path / "expr.gct"
    _write_gct(p, ROWS)
    on = GTExExpressionTable(p, symbol_fallback=True)
    off = GTExExpressionTable(p, symbol_fallback=False)
    assert on.lookup("ENSGDRIFTED", gene_symbol="GENEB") == pytest.approx(math.log2(7.5))
    assert off.lookup("ENSGDRIFTED", gene_symbol="GENEB") is None
    # ID match wins over symbol
    assert on.lookup("ENSG00000000001", gene_symbol="GENEB") == pytest.approx(0.0)


def test_expression_duplicate_stripped_ids_keep_first(tmp_path):
    p = tmp_path / "expr.gct"
    rows = [("ENSG00000000009.1", "X", [1.0] * 12), ("ENSG00000000009.1_PAR_Y", "Y", [100.0] * 12)]
    _write_gct(p, rows)
    assert GTExExpressionTable(p).lookup("ENSG00000000009") == pytest.approx(1.0)


def test_expression_rejects_non_gtex_table(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("Name\tDescription\tT1\tT2\nENSG1\tA\t1\t2\n")
    with pytest.raises(ValueError):
        GTExExpressionTable(p)


# ---------------------------------------------------------------------------
# FeatureAnnotator end to end (toy transcript, fake conservation / half-life)
# ---------------------------------------------------------------------------
class _FakePhastcons:
    def median_over_blocks(self, blocks):
        # real signature: blocks are (chrom, start_1based, end_1based); returns (median, total_bp, valid_bp)
        total = sum(e - s + 1 for _chrom, s, e in blocks)
        return 0.42, total, total


class _FakeHalfLife:
    def lookup(self, ensg, strip_versions=True, gene_symbol=None):
        return 1.25


def test_annotator_returns_exactly_the_eight_features(tmp_path):
    pyfaidx = pytest.importorskip("pyfaidx")
    genome = _toy_genome()
    fa = pyfaidx.Fasta(str(_write_fasta(tmp_path, genome)), sequence_always_upper=True)
    expr_path = tmp_path / "expr.gct"
    _write_gct(expr_path, [("ENSGTOY01.1", "TOYGENE", list(range(1, 13)))])
    tx = _toy_tx("+")
    annot = FeatureAnnotator(
        tx_index={tx.transcript_id: tx}, fasta=fa, phastcons=_FakePhastcons(),
        halflife=_FakeHalfLife(), expression=GTExExpressionTable(expr_path),
    )
    res = annot.annotate({"variant_id": "v1", "txnames": "ENSTTOY01.1", "gene": "TOYGENE",
                          "contig": "chr1", "position": 350})
    row = res.to_feature_row()

    ids = {"variant_id", "txnames", "transcript_id_used", "gene", "gene_id", "strand"}
    assert set(row) - ids == set(REQUIRED_FEATURES)
    assert not (set(row) & RETIRED)

    cds = _expected_cds(genome, tx.cds_with_stop(), "+")
    assert row["cdsseqs_AU_content"] == pytest.approx(_au_content(cds))
    assert row["cdsseqs_UC_content"] == pytest.approx(_uc_content(cds))
    assert row["MedianExpression_log2"] == pytest.approx(math.log2(6.5 + 1))
    assert row["half_life_PC1"] == pytest.approx(1.25)
    assert row["phastcons_new3utr_first200_median"] == pytest.approx(0.42)
    assert row["last.EJC"] == "penultimate.last50bp" and row["mut.exon"] == 2

    qc = res.to_qc_row()
    assert qc["expression_missing"] is False and "ptc_to_ejc_empty" not in qc


def test_annotator_flags_missing_expression(tmp_path):
    pyfaidx = pytest.importorskip("pyfaidx")
    fa = pyfaidx.Fasta(str(_write_fasta(tmp_path, _toy_genome())), sequence_always_upper=True)
    expr_path = tmp_path / "expr.gct"
    _write_gct(expr_path, [("ENSGOTHER.1", "OTHER", [1.0] * 12)])
    tx = _toy_tx("+")
    annot = FeatureAnnotator({tx.transcript_id: tx}, fa, _FakePhastcons(), _FakeHalfLife(),
                             GTExExpressionTable(expr_path))
    res = annot.annotate({"variant_id": "v1", "txnames": "ENSTTOY01", "gene": "NOTHERE",
                          "contig": "chr1", "position": 350})
    assert res.MedianExpression_log2 is None and res.expression_missing is True


# ---------------------------------------------------------------------------
# predict.py imputation (needs catboost installed)
# ---------------------------------------------------------------------------
def test_predict_imputation_uses_training_medians_and_does_not_zero_fill(tmp_path):
    pytest.importorskip("catboost")
    pytest.importorskip("joblib")
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import predict as tk_predict

    medians = tmp_path / "medians.json"
    medians.write_text(json.dumps({"MedianExpression_log2": 3.83, "half_life_PC1": 1.13, "CADD_phred": 37.0}))
    X = pd.DataFrame({
        "half_life_PC1": [np.nan, 0.5],
        "MedianExpression_log2": [np.nan, 4.0],
        "phastcons_new3utr_first200_median": [np.nan, 0.3],
    })
    out, report = tk_predict.apply_training_imputation(X.copy(), medians)
    assert out["MedianExpression_log2"].tolist() == [3.83, 4.0]
    assert out["half_life_PC1"].tolist() == [1.13, 0.5]
    assert math.isnan(out.loc[0, "phastcons_new3utr_first200_median"])   # not zero-filled (matches training)
    assert set(report["imputed"]) == {"half_life_PC1", "MedianExpression_log2"}
    assert tk_predict.ZERO_FILL_COLUMNS == []
