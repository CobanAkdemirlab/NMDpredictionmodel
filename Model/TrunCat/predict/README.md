# Prediction Pipeline

Reusable inference for TrunCat on new variant cohorts (ClinVar, gnomAD,
GREGoR, or any other TOPMed-style stop-gain cohort).

## Files

| File | Purpose |
|------|---------|
| `predict_cohorts.ipynb` | Main entry point. Annotates and scores the ClinVar, gnomAD and GREGoR tables, with audits and distribution plots. |
| `predict.py` | Inference module and command-line script. Takes a cohort table and produces a per-variant escape probability and threshold call. |
| `training_medians.json` | Training-set medians used for the explicit median-imputation block, so new cohorts are imputed exactly as the training data were. |
| `clinvar_predictions.csv`, `gnomAD_predictions.csv`, `gregor_predictions.csv` | Current cohort predictions from `predict_cohorts.ipynb`. |

Older notebooks (`predict_new_variants`, `predict_gregor`) are superseded by
`predict_cohorts.ipynb`.

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
2. Uses **precomputed training medians** (`training_medians.json`) for the
   explicit median-impute block (`pLI`, `MedianExpression_log2`,
   `half_life_PC1`, `CADD_phred`, `readthrough_score_hek293t`, etc.)
3. Aligns columns to `model.feature_names_`, adding missing columns as NaN
   (CatBoost handles those natively) and dropping extras.

Column order follows the model, not the input table. CatBoost results shift
slightly with column order, so do not reorder after alignment.

## Workflow

### Cohorts (notebook)

Open `predict_cohorts.ipynb` and run all cells. Annotation outputs for each
cohort are written to `predict/annotations/`. Cohort annotation uses the same
feature-generation scripts as training (see `../../../Features/`), run
per cohort.

### Single cohort (command line)

```bash
python predict.py --help
```

lists the current options. A typical call on a cohort table that already has
all annotations:

```bash
python predict.py \
    --input  path/to/cohort_table.csv \
    --output path/to/predict/cohort_predictions.csv \
    --label  cohort_name \
    --threshold 0.4747 \
    --audit-json path/to/predict/cohort_audit.json
```

`--threshold` is required unless a `viz_summary.json` is available at the
`--viz-summary` path (default `results/figures/viz_summary.json`, relative to
the working directory; the repository does not ship this file). 0.4747 is the
Youden-optimal threshold from the TrunCat out-of-fold predictions.

## Cohort input notes

Cohort tables built before the corrected-features rerun can differ from the
training table. `predict_cohorts.ipynb` handles the known cases:

- `EJC.downstream.cut` bins unseen at training time are mapped to the training levels.
- `.` placeholder strings in numeric gnomAD columns are zero-filled, as in training.
- Missing intron flags are filled with the training level (`3UTR-No`).

If you add a new cohort, check these first.

## Output format

Each prediction CSV contains:

| Column | Description |
|--------|-------------|
| `variantID` | Primary key (chr_pos_ref_alt) |
| `escape_probability` | Model-predicted probability of NMD escape |
| `predicted_class` | 0 = NMD-sensitive, 1 = escape (at threshold) |
| `predicted_label` | `'NMD'` or `'escape'` |
| `threshold_used` | Youden threshold used for the call (0.4747 for the current model; set with `--threshold`) |
| `GENE_ID`, `hgnc_symbol`, `CHROM`, `POS`, `REF_ALLELE`, `ALT_ALLELE` | Identifier columns carried through from the input (the cohort files pass these via `--extra-id-cols`) |

Predicted probabilities are on the TOPMed out-of-fold scale; see
`../results/calibration/` for how well that scale matches observed escape
rates, and `../notebooks/tier_analysis_clinvar_combined.ipynb` for confidence
tiers.

## Audit report

`predict.py` prints (and optionally writes as JSON via `--audit-json`) a
report covering:

- Feature count present / missing / extra vs the model's expectation
- Dtype coercions applied during alignment
- Structural zero-fills (feature name + NaN count)
- gnomAD zero-fill count
- Median imputations (feature, NaN count, fill value)
- Residual NaNs per feature after all cleaning (these go to CatBoost's native
  NaN handling, which works, but is worth looking at to catch upstream bugs)

**Rule of thumb**: if `missing` is non-zero, either (a) your feature-generation
pipeline didn't run to completion, or (b) notebook 02 dropped that feature at
training time and it shouldn't even be in the model's `feature_names_`. Check
which case before trusting the predictions.

## Categorical value spot check (important)

The engineered categoricals (`aug_distance_category`, `kozak_strength`,
`aug_frame_status`, `has_plus1_aug`, `has_plus2_aug`) have explicit domain
categories encoding real biology: `no_inframe_AUG`, `no_frame_AUG`, and the
distance bins. You should **never** see the generic `'MISSING'` string in
these columns. If `predict_cohorts.ipynb` flags any, re-check the PTC AUG
annotation run for that cohort.
