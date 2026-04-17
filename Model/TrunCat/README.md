# TrunCat

**TRUN**cation-aware **C**lassifier using **A**nnotated **T**ranscripts.

TrunCat is the full NMD-escape prediction pipeline: a CatBoost classifier
trained on ~5,700 premature termination codon (PTC) variants from TOPMed,
using ~730 transcript-level, genomic, and sequence-based features to predict
whether stopgain variants will escape nonsense-mediated decay surveillance.

**Out-of-fold ROC-AUC ≈ 0.78** | **Youden-optimal threshold ≈ 0.38–0.41**

## Layout

Model/TrunCat/
├── config/
│   └── config.yaml          ← paths, hyperparameters, feature lists
├── data/                    ← input variant and annotation files (not distributed)
│   └── README.md            ← sourcing instructions
├── notebooks/               ← three-stage pipeline (see notebooks/README.md)
│   ├── 01_data_loading_and_merging.ipynb
│   ├── 02_feature_cleaning_and_selection.ipynb
│   └── 03_model_training.ipynb
├── scripts/                 ← standalone Python equivalents (see scripts/README.md)
├── model/                   ← trained TrunCat artifacts (tracked)
│   ├── truncat.cbm
│   └── truncat.pkl
└── results/                 ← CV outputs, figures (regenerated; not tracked)

## Pipeline

Notebooks and scripts share `config/config.yaml` and produce identical
outputs. Use whichever interface fits your workflow.

1. **Notebook 01** — merge TOPMed variants with six annotation sources into `data/TOPMed_merged.csv`
2. **Notebook 02** — feature cleaning, leakage removal, correlation filtering → `data/TOPMed_cleaned.csv`
3. **Notebook 03** — 5-fold stratified CV training, SHAP analysis, final model → `model/truncat.{cbm,pkl}` + `results/`

See [`notebooks/README.md`](notebooks/README.md) for a per-notebook walkthrough.

## Top predictive features (mean |SHAP|)

`last.EJC`, `relativePTClocation`, `half_life_PC1`, `cdsseqs_AU_content`,
`mut.exon`, `phastcons_new3utr_first200_median`, `phylop_ptc_to_ejc_median`,
`AmountExonsAfter`, `cdsseq_AUcontentlast200`, `cdsseqs_UC_content`.

## Related

For a lightweight version using only the top 10 features — suitable for
scoring external cohorts without reproducing the full annotation pipeline —
see [`../TrunKitten/`](../TrunKitten/).