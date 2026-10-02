# TrunKitten

A reduced version of [TrunCat](../TrunCat/) that uses only 8 features. Designed for
prediction on external cohorts where reproducing the full 619-feature annotation
pipeline is impractical.

## Layout

- `config/config.yaml` — paths, CatBoost hyperparameters (TrunKitten's own), feature lists
- `notebooks/04_reduced_model_top10_shap.ipynb` — selection of the feature set, training and evaluation of the reduced model
- `inputs/` — SHAP-rankings snapshot, half-life table, example variants, and download instructions for the large reference files (`inputs/README.md`)
- `model/` — trained model artifacts
  - `trunkitten.cbm`, `trunkitten.pkl`
  - `trunkitten_features.json` — feature order, categorical specification, Youden threshold
- `results/` — cross-validated predictions and performance summary
- `pipeline/` — standalone annotation + prediction pipeline for scoring new variants

## The 8 features

`last.EJC`, `relativePTClocation`, `half_life_PC1`, `cdsseqs_AU_content`, `mut.exon`,
`cdsseqs_UC_content`, `phastcons_new3utr_first200_median`, `MedianExpression_log2`.

**How they were chosen.** The candidates are TrunCat's top 10 features by mean |SHAP|.
A candidate is kept only if it also appears in the top 10 of all five
leave-one-fold-out SHAP rankings (the table is
`../TrunCat/results/lofo_stability/lofo_feature_stability_table.csv`). Eight features
meet that rule. The feature order above is the model's column order and is fixed in
`trunkitten_features.json`.

## Scoring new variants

The `pipeline/` subdirectory contains a standalone annotation and prediction
workflow: annotate a variant TSV with the 8 features TrunKitten requires, then score
with the trained model.

```bash
cd pipeline/

# Annotate variants (paths and reference files are set in config/config.yaml)
python -m minicat.cli --config config/config.yaml

# Score with TrunKitten
python predict.py \
    --annotated outputs/annotated.tsv \
    --model     ../model/trunkitten.pkl \
    --metadata  ../model/trunkitten_features.json \
    --out       outputs/predictions.tsv \
    --training-medians ../../TrunCat/predict/training_medians.json
```

Besides a GTF, genome FASTA, a phastCons BigWig and the mRNA half-life table, the
annotation step needs the GTEx v8 median-TPM table for `MedianExpression_log2`.
See [`pipeline/README.md`](pipeline/README.md) for inputs, the optional
training-medians imputation, and validation against the training features.

## Relationship to TrunCat

TrunKitten is trained on the same TOPMed variant set as TrunCat
(`../TrunCat/data/TOPMed_cleaned.csv`) with gene-grouped 5-fold cross-validation and
a fixed number of trees (no early stopping). It uses its own tuned CatBoost
hyperparameters (`config/config.yaml`).

With 8 of TrunCat's 619 features, TrunKitten reaches an out-of-fold ROC-AUC of 0.7754,
against 0.7773 for the full model (99.8% retained; paired gene-clustered bootstrap
difference −0.0019, 95% CI −0.0072 to +0.0033). This makes TrunKitten suitable for
scoring external cohorts where reproducing the full annotation pipeline is
impractical.

For the full model and feature engineering pipeline, see [`../TrunCat/`](../TrunCat/).
