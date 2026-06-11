# Prediction Pipeline

Reusable inference for `spooky_model_v4.1` on new variant cohorts (ClinVar,
gnomAD, or any other TOPMed-style stopgain cohort).

## Files

| File                                       | Purpose |
|--------------------------------------------|---------|
| `scripts/predict.py`                       | Main inference script. Takes a cohort CSV and produces per-variant escape probability + threshold call. |
| `scripts/export_training_medians.py`       | One-off helper. Extracts the medians notebook 02 uses for imputation and persists them as JSON so `predict.py` can reuse them identically. |
| `notebooks/04_predict_new_variants.ipynb`  | Interactive walk-through. Loops over cohorts, runs `predict.py`'s helpers, plots distributions, verifies output columns. |

## Design principle

The trained model's `feature_names_` is the **source of truth** for which
features to use and in what order. `predict.py` never re-runs the
data-dependent cleaning from notebook 02 (correlation filtering,
zero-variance detection, cohort-wise medians) because those decisions are
already baked into the model. Instead it:

1. Applies the **deterministic** rules from notebook 02:
   - `loefu_cat` / `LOEUF_cat` → `loeuf_cat` rename
   - Structural zero-fills (UTR windows, conservation in nonexistent regions)
   - gnomAD zero-fill (absent = ultra-rare)
   - Categorical cast to string, NaN → `'MISSING'`
2. Uses **precomputed training medians** for the explicit median-impute block
   (`pLI`, `MedianExpression_log2`, `half_life_PC1`, `CADD_phred`,
   `readthrough_score_hek293t`, etc.)
3. Aligns columns to `model.feature_names_`, adding missing columns as NaN
   (CatBoost handles those natively) and dropping extras.

## Workflow

### One-time setup

After notebook 02 finishes and `TOPMed_cleaned_v4.csv` exists:

```bash
python scripts/export_training_medians.py \
    --cleaned-csv /Users/jschmidt3/Iman_visualizations/spooky_model_v4.1/TOPMed_cleaned_v4.csv \
    --out models/training_medians.json
```

This writes one number per feature that notebook 02 median-imputes. Re-run
only if the training set changes.

### Per cohort (command-line)

If the cohort CSV is already fully merged (has all v4 annotations):

```bash
python scripts/predict.py \
    --input  /Users/jschmidt3/Iman_visualizations/spooky_model_v4.1/predict/clinvar_df_stopgain_updated2026.csv \
    --output /Users/jschmidt3/Iman_visualizations/spooky_model_v4.1/predict/clinvar_predictions_v4.1.csv \
    --label  clinvar_2026 \
    --extra-id-cols GENE_ID hgnc_symbol clinvar_significance \
    --audit-json /Users/jschmidt3/Iman_visualizations/spooky_model_v4.1/predict/clinvar_audit_v4.1.json
```

If the cohort CSV only has the baseline features and the v4 annotation files
are separate (TSV/CSV outputs of the feature-generation scripts), pass them
explicitly and `predict.py` will merge them in using notebook 01's logic:

```bash
python scripts/predict.py \
    --input         .../predict/clinvar_df_stopgain_updated2026.csv \
    --output        .../predict/clinvar_predictions_v4.1.csv \
    --codon-opt     .../predict/annotations/clinvar_codon_optimality.tsv \
    --readthrough   .../predict/annotations/clinvar_readthrough.csv \
    --ejc           .../predict/annotations/clinvar_ejc_occupancy.tsv \
    --ptc-aug       .../predict/annotations/clinvar_ptc_aug.tsv \
    --conservation  .../predict/annotations/clinvar_conservation.csv \
    --label         clinvar_2026
```

### Per cohort (notebook)

Open `notebooks/04_predict_new_variants.ipynb`, edit the `COHORTS` dict at the
top (step 1) to point at your cohort files, and run all cells. The notebook
walks through:

1. Config + path setup
2. Load model, threshold, training medians
3. Inspect model expectations (feature count, categorical list, top importances)
4. Loop over cohorts: merge → clean + align → audit → predict → save
5. Column verification summary table
6. Categorical value spot check (detects upstream feature-generation bugs)
7. Cohort prediction distribution plots
8. Sample head/tail predictions
9. Final column checklist (pass/fail per cohort)

## Output format

Each prediction CSV contains:

| Column                | Description |
|-----------------------|-------------|
| `variantID`           | Primary key (chr_pos_ref_alt) |
| `escape_probability`  | Model-predicted probability of NMD escape |
| `predicted_class`     | 0 = NMD-sensitive, 1 = escape (at threshold) |
| `predicted_label`     | `'NMD'` or `'escape'` |
| `threshold_used`      | Youden threshold from `viz_summary.json` |
| *(extra id cols)*     | Any columns passed via `--extra-id-cols` |

## Audit report

`predict.py` prints (and optionally writes as JSON via `--audit-json`) a
report covering:

- Feature count present / missing / extra vs the model's expectation
- Dtype coercions applied during alignment
- Structural zero-fills (feature name + NaN count)
- gnomAD zero-fill count
- Median imputations (feature, NaN count, fill value)
- Residual NaNs per feature after all cleaning (these go to CatBoost's native
  NaN handling — which works, but is worth looking at to catch upstream bugs)

**Rule of thumb**: if `missing` is non-zero, either (a) your feature-generation
pipeline didn't run to completion, or (b) notebook 02 dropped that feature at
training time and it shouldn't even be in the model's `feature_names_` — check
which case before trusting the predictions.

## Categorical value spot check (important)

The engineered v4 categoricals (`aug_distance_category`, `kozak_strength`,
`aug_frame_status`, `has_plus1_aug`, `has_plus2_aug`) have explicit domain
categories encoding real biology — `no_inframe_AUG`, `no_frame_AUG`, and the
distance bins. You should **never** see the generic `'MISSING'` string in
these columns. If the notebook's section 6 flags any, re-check the PTC AUG
annotation run for that cohort.
