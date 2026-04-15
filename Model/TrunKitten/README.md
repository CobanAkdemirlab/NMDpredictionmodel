# TrunKitten

A reduced version of [TrunCat](../TrunCat/) using only the 10 most informative
transcript features (by mean |SHAP|). Designed for prediction on external
cohorts where reproducing the full ~730-feature annotation pipeline is
impractical.

## Layout

- `config/config.yaml` — TrunKitten configuration (paths, hyperparameters, feature lists)
- `notebooks/04_reduced_model_top10_shap.ipynb` — training and evaluation of the reduced model
- `inputs/shap_feature_importance_rankings.csv` — snapshot of TrunCat's SHAP rankings used to define the top-10 selection
- `model/reduced_top10/` — trained model artifacts
  - `catboost_reduced_top10.cbm`, `catboost_reduced_top10.pkl`
  - `top10_features.json` — feature order, categorical specification, Youden threshold
- `pipeline/` — standalone annotation + prediction pipeline for scoring new variants
  - See `pipeline/README.md` for the annotation CLI and `predict.py` invocation

## Top 10 features

`last.EJC`, `relativePTClocation`, `half_life_PC1`, `cdsseqs_AU_content`,
`mut.exon`, `phastcons_new3utr_first200_median`, `phylop_ptc_to_ejc_median`,
`AmountExonsAfter`, `cdsseq_AUcontentlast200`, `cdsseqs_UC_content`.

## Relationship to TrunCat

TrunKitten is trained on the same TOPMed variant set as TrunCat
(`../TrunCat/data/TOPMed_cleaned.csv`), restricted to the top 10 features by
SHAP importance from TrunCat's CV-averaged ranking. Hyperparameters match
TrunCat v4.1. It trades some predictive performance for substantially lower
annotation burden on new inputs.

For the full model and feature engineering pipeline, see [`../TrunCat/`](../TrunCat/).
