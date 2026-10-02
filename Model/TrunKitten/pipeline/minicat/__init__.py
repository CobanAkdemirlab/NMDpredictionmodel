"""TrunKitten PTC annotation pipeline (package name: `minicat`).

Produces the 8 features required by TrunKitten — the reduced 8-feature
(leave-one-fold-out-stable) NMD-prediction model derived from TrunCat (TRUNcation-aware Classifier using
Annotated Transcripts) — for externally-called stop-gain variants.

Conventions match the TrunCat training pipeline in
CobanAkdemirLab/NMDpredictionmodel.

The package is named `minicat` for historical reasons; `python -m minicat.cli`
and all imports remain unchanged. "TrunKitten" is the user-facing pipeline
name.
"""

__version__ = "0.2.0"

# Order must match Model/TrunKitten/model/trunkitten_features.json["features_in_order"]
# (enforced by tests/test_feature_contract.py). Column order affects CatBoost results.
REQUIRED_FEATURES = [
    "last.EJC",
    "relativePTClocation",
    "half_life_PC1",
    "cdsseqs_AU_content",
    "mut.exon",
    "cdsseqs_UC_content",
    "phastcons_new3utr_first200_median",
    "MedianExpression_log2",
]
