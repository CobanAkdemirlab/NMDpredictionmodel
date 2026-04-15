"""TrunKitten PTC annotation pipeline (package name: `minicat`).

Produces the 10 features required by TrunKitten — the reduced top-10 feature
NMD-prediction model derived from TrunCat (TRUNcation-aware Classifier using
Annotated Transcripts) — for externally-called stop-gain variants.

Conventions match the TrunCat training pipeline in
CobanAkdemirLab/NMDpredictionmodel.

The package is named `minicat` for historical reasons; `python -m minicat.cli`
and all imports remain unchanged. "TrunKitten" is the user-facing pipeline
name.
"""

__version__ = "0.1.0"

REQUIRED_FEATURES = [
    "last.EJC",
    "relativePTClocation",
    "half_life_PC1",
    "cdsseqs_AU_content",
    "mut.exon",
    "phastcons_new3utr_first200_median",
    "phylop_ptc_to_ejc_median",
    "AmountExonsAfter",
    "cdsseq_AUcontentlast200",
    "cdsseqs_UC_content",
]
