# TrunKitten

A reduced version of [TrunCat](../TrunCat/) using only the 10 most informative
features (by mean |SHAP|). Designed for prediction on external cohorts where
reproducing the full ~850-feature annotation pipeline is impractical.

## Layout

- `config/config.yaml` — paths, hyperparameters, feature lists
- `notebooks/04_reduced_model_top10_shap.ipynb` — training and evaluation of the reduced model
- `inputs/shap_feature_importance_rankings.csv` — snapshot of TrunCat's SHAP rankings used to define the top-10 selection
- `model/` — trained model artifacts
  - `trunkitten.cbm`, `trunkitten.pkl`
  - `trunkitten_features.json` — feature order, categorical specification, Youden threshold
- `pipeline/` — standalone annotation + prediction pipeline for scoring new variants

## Top 10 features

`last.EJC`, `relativePTClocation`, `half_life_PC1`, `cdsseqs_AU_content`,
`mut.exon`, `phastcons_new3utr_first200_median`, `phylop_ptc_to_ejc_median`,
`AmountExonsAfter`, `cdsseq_AUcontentlast200`, `cdsseqs_UC_content`.

## Scoring new variants

The `pipeline/` subdirectory contains a standalone annotation and prediction
workflow: annotate a variant TSV with the 10 features TrunKitten requires,
then score with the trained model.

```bash
cd pipeline/

# Annotate variants (produces annotated.tsv)
python -m minicat.cli --input inputs/your_variants.tsv --out outputs/annotated.tsv

# Score with TrunKitten
python predict.py \
    --annotated outputs/annotated.tsv \
    --model     ../model/trunkitten.pkl \
    --metadata  ../model/trunkitten_features.json \
    --out       outputs/predictions.tsv
```

See [`pipeline/README.md`](pipeline/README.md) for full usage, input TSV format,
optional training-medians imputation, and example outputs.

## Relationship to TrunCat

TrunKitten is trained on the same TOPMed variant set as TrunCat
(`../TrunCat/data/TOPMed_cleaned.csv`), restricted to the top 10 features by
SHAP importance from TrunCat's CV-averaged ranking. Hyperparameters match
TrunCat's current configuration.

Despite using only ~1.2% of the features, TrunKitten retains 99% of TrunCat's
out-of-fold ROC-AUC (0.779 vs. 0.785 on the full 853-feature model).
This makes TrunKitten suitable for scoring external cohorts where reproducing
the full annotation pipeline is impractical.

For the full model and feature engineering pipeline, see [`../TrunCat/`](../TrunCat/).
