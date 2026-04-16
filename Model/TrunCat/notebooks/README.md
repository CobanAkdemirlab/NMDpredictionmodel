# TrunCat Notebooks

Three-stage pipeline: data merging → feature cleaning → model training.
All three share `../config/config.yaml`. Paths resolve via `BASE_DIR`
anchored to `Model/TrunCat/`, so notebooks work regardless of Jupyter
launch directory.

Run in order:

```bash
jupyter notebook 01_data_loading_and_merging.ipynb
jupyter notebook 02_feature_cleaning_and_selection.ipynb
jupyter notebook 03_model_training.ipynb
```

---

## 01 — Data loading and merging

Loads six input files, filters indel variants, engineers AUG features, and
merges into a single training table.

**Inputs** (from `config.data.*`):

| Key | File | Description |
|-----|------|-------------|
| `enhanced_features` | `data/TOPMed_stopgain.csv` | Main variant CSV with NMD features and expression |
| `codon_optimality` | `data/annotations/codon_optimality.tsv` | CodonOptimalityFraction scores |
| `readthrough` | `data/annotations/readthrough.csv` | HEK293T readthrough scores and categories |
| `ejc` | `data/annotations/ejc_occupancy.tsv` | EJC occupancy counts and overlap flags |
| `ptc_aug` | `data/annotations/ptc_aug.tsv` | Downstream in-frame AUG features |
| `conservation` | `data/annotations/conservation_medians.csv` | PhastCons and PhyloP regional medians |

**Output:** `data/TOPMed_merged.csv`

**Key operations:** indel removal via V7/V8 allele scanning; `variantID`
key standardization; deprecated feature drop (`AverageCodonRNAUsage`,
`median_half_life`); AUG distance/Kozak/frame engineering; left-join merge.

---

## 02 — Feature cleaning and selection

Domain drop rules, automated quality checks, imputation.

**Input:** `data/TOPMed_merged.csv`
**Outputs:** `data/TOPMed_cleaned.csv`, `data/final_feature_list.csv`

| Rule | Decision |
|------|----------|
| AU vs GC content | Keep AU (biologically relevant), drop GC |
| RBP motif features | Protect all ~700 from correlation removal |
| PTC distance | Keep `relativePTClocation` (normalized), drop raw distances |
| Leaky features | Drop `ALLELE.RAT`, `refCount`, `altCount`, `Freq`, `Whole.Blood` |
| Expression | Keep `MedianExpression_log2`, drop raw `MedianExpression` |
| UTR composition | Keep whole + first100 + last100 windows; drop 200nt windows |
| Identifiers | Drop `variantID`, `GENE_ID`, coordinates, allele strings |
| Zero-importance RBP | Drop 109 features confirmed zero importance across all CV folds |

Automated checks then remove duplicate columns, zero-variance features, and
correlated features (Pearson |r| > 0.95, with RBP protection). Missing values
imputed by structural zero (absent UTR regions, gnomAD absence) or median.

---

## 03 — Model training

Trains the CatBoost classifier, computes SHAP-based feature importances,
generates publication figures.

**Input:** `data/TOPMed_cleaned.csv`

**Outputs:**

model/
├── truncat.cbm           ← final model (all data)
└── truncat.pkl

results/
├── models/cv_folds/                      ← per-fold models
├── cv_predictions.csv                    ← out-of-fold predicted probabilities
├── feature_importances_cv_averaged.csv
├── model_performance_summary.txt
└── figures/visualizations_cv/
├── 1_roc_curve.{png,pdf}
├── 2_precision_recall_curve.{png,pdf}
├── 3_confusion_matrix.{png,pdf}
├── 4_probability_distribution.{png,pdf}
├── 5_performance_dashboard.{png,pdf}
├── 7a_feature_importance_shap_top20.{png,pdf}
├── 7b_feature_importance_native_top20.{png,pdf}
├── 7c_feature_importance_comparison.{png,pdf}
└── shap_manuscript/
├── 6_shap_summary_cv_averaged.{png,pdf}
├── 8a_shap_waterfall_escape_example.{png,pdf}
├── 8b_shap_waterfall_sensitive_example.{png,pdf}
├── shap_feature_importance_rankings.csv
└── feature_importance_comparison.csv

**Key operations:**
- 5-fold stratified CV with out-of-fold predictions for unbiased metrics
- Native CatBoost importances averaged across all folds
- SHAP values computed per-fold (500-sample subsets) then aggregated
- Youden-optimal threshold derived from OOF predictions
- Final model trained on complete dataset for deployment

---

## Model design notes

| Choice | Rationale |
|--------|-----------|
| CatBoost classifier | Native categorical support, robust to missing values, ordered boosting |
| 5-fold stratified CV | Unbiased AUC, consistent with published NMD predictor benchmarks |
| OOF predictions | Unbiased probability estimates for threshold selection |
| Mean \|SHAP\| for importance | Fairer than native importance for correlated feature sets |
| Youden-optimal threshold | Balances sensitivity and specificity |
| Final model on all data | Maximizes training signal for genome-wide deployment |