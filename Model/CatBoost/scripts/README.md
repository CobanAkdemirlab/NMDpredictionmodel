# NMD Escape Prediction Model

Machine learning pipeline for predicting nonsense-mediated decay (NMD) escape in stopgain variants carrying premature termination codons (PTCs).

## Overview

This repository contains the complete workflow for training and evaluating a CatBoost classifier that predicts whether stopgain variants will escape NMD surveillance. The model incorporates positional information (last exon rule, relative PTC location), conservation scores, RNA-binding protein motifs, codon optimality, mRNA half-life, and transcript characteristics.

**Out-of-fold ROC-AUC: ~0.78 | Optimal threshold: ~0.38–0.41 (Youden's J)**

---

## Repository Structure

```
repo_root/
├── config/
│   └── config.yaml                          ← shared config for notebooks and scripts
├── data/
│   ├── README.md                            ← data sourcing instructions
│   ├── TOPMed_stopgain.csv                  ← main variant CSV (not distributed)
│   └── annotations/
│       ├── codon_optimality.tsv
│       ├── readthrough.csv
│       ├── ejc_occupancy.tsv
│       ├── ptc_aug.tsv
│       └── conservation_medians.csv
├── notebooks/
│   ├── 01_data_loading_and_merging.ipynb
│   ├── 02_feature_cleaning_and_selection.ipynb
│   └── 03_model_training.ipynb
├── scripts/                                 ← standalone Python scripts (this folder)
│   ├── 01_data_loading_and_merging.py
│   ├── 02_feature_cleaning_and_selection.py
│   ├── 03_model_training.py
│   ├── setup.py
│   ├── requirements.txt
│   └── README.md
└── results/                                 ← created automatically on first run
    ├── models/
    │   ├── catboost_final_model.cbm
    │   ├── catboost_final_model.pkl
    │   └── cv_folds/
    ├── cv_predictions.csv
    ├── cv_metrics.csv
    ├── feature_importances_cv_averaged.csv
    ├── model_performance_summary.txt
    └── figures/
        └── visualizations_cv/
            └── shap_manuscript/
```

The scripts and notebooks share the same `config/config.yaml` and produce identical outputs — use whichever interface suits your workflow.

---

## Requirements

- Python 3.12+
- See `requirements.txt` for full package list

### Installation

```bash
# Create a virtual environment (recommended)
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Verify setup

```bash
python setup.py
```

This checks that the config is valid, all input data files are present, output directories exist, and all dependencies are installed.

---

## Configuration

All paths and hyperparameters are controlled by `config/config.yaml` at the repo root. Paths in the config are relative to the repo root and are resolved automatically — no edits needed unless you change the directory layout or update your data files.

Key sections:

```yaml
data:
  enhanced_features: "data/TOPMed_stopgain.csv"
  codon_optimality:  "data/annotations/codon_optimality.tsv"
  readthrough:       "data/annotations/readthrough.csv"
  ejc:               "data/annotations/ejc_occupancy.tsv"
  ptc_aug:           "data/annotations/ptc_aug.tsv"
  conservation:      "data/annotations/conservation_medians.csv"
  merged:            "data/TOPMed_merged.csv"
  cleaned:           "data/TOPMed_cleaned.csv"

model:
  target: "NMD.ESCAPEE"
  random_seed: 42
  n_folds: 5
  catboost:
    iterations: 1294
    learning_rate: 0.03072886163898516
    depth: 8
    # ... (Optuna-optimized — do not change for exact replication)
```

See `data/README.md` for instructions on obtaining the input files.

---

## Usage

Run the three scripts in order. Each accepts an optional `--config` flag if you need to point to a non-default config location.

### Step 1 — Data loading and merging

```bash
python 01_data_loading_and_merging.py
```

Loads six input files, filters indel variants, engineers AUG features, and merges everything into a single dataset.

**Inputs** (from `config.data.*`):
| Key | File | Description |
|-----|------|-------------|
| `enhanced_features` | `TOPMed_stopgain.csv` | Main variant CSV with NMD features and expression |
| `codon_optimality` | `annotations/codon_optimality.tsv` | CodonOptimalityFraction scores |
| `readthrough` | `annotations/readthrough.csv` | HEK293T readthrough scores and categories |
| `ejc` | `annotations/ejc_occupancy.tsv` | EJC occupancy counts and overlap flags |
| `ptc_aug` | `annotations/ptc_aug.tsv` | Downstream in-frame AUG features |
| `conservation` | `annotations/conservation_medians.csv` | PhastCons and PhyloP regional medians |

**Output:** `data/TOPMed_merged.csv`

**Key operations:**
- Identifies and removes indel variants via V7/V8 allele scanning
- Standardizes join keys (`variantID`) across all annotation files
- Drops deprecated features (`AverageCodonRNAUsage`, `median_half_life`)
- Engineers AUG distance categories, Kozak strength, and frame status
- Merges all sources on `variantID` (left join from main variant CSV)

---

### Step 2 — Feature cleaning and selection

```bash
python 02_feature_cleaning_and_selection.py
```

Applies domain-knowledge drop rules, automated quality checks, and imputation.

**Input:** `data/TOPMed_merged.csv`

**Outputs:** `data/TOPMed_cleaned.csv`, `data/final_feature_list.csv`

**Key operations:**

| Rule | Decision |
|------|----------|
| AU vs GC content | Keep AU (biologically relevant), drop GC |
| RBP motif features | Protect all ~700 features from correlation removal |
| PTC distance | Keep `relativePTClocation` (normalized), drop raw distances |
| Leaky features | Drop `ALLELE.RAT`, `refCount`, `altCount`, `Freq`, `Whole.Blood` |
| Expression | Keep `MedianExpression_log2`, drop raw `MedianExpression` |
| UTR composition | Keep whole + first100 + last100 windows; drop 200nt windows |
| Identifiers | Drop `variantID`, `GENE_ID`, coordinates, allele strings |
| Zero-importance RBP | Drop 109 features confirmed zero importance across all CV folds |

Automated checks then remove duplicate columns, zero-variance features, and features with Pearson |r| > 0.95 (with RBP protection). Missing values are imputed by structural zeros (absent UTR regions, gnomAD absence) or median (annotation gaps).

---

### Step 3 — Model training

```bash
python 03_model_training.py
```

Trains the CatBoost classifier, computes SHAP-based feature importances, and generates all publication-quality figures.

**Input:** `data/TOPMed_cleaned.csv`

**Outputs:**

```
results/
├── models/
│   ├── catboost_final_model.cbm         ← final model (all data)
│   ├── catboost_final_model.pkl
│   └── cv_folds/
│       ├── fold_1_auc_*.cbm             ← per-fold models
│       └── ...
├── cv_predictions.csv                   ← out-of-fold predicted probabilities
├── feature_importances_cv_averaged.csv  ← native CatBoost importances (CV-averaged)
├── model_performance_summary.txt
└── figures/
    └── visualizations_cv/
        ├── 1_roc_curve.png/pdf
        ├── 2_precision_recall_curve.png/pdf
        ├── 3_confusion_matrix.png/pdf
        ├── 4_probability_distribution.png/pdf
        ├── 5_performance_dashboard.png/pdf
        ├── 7a_feature_importance_shap_top20.png/pdf
        ├── 7b_feature_importance_native_top20.png/pdf
        ├── 7c_feature_importance_comparison.png/pdf
        └── shap_manuscript/
            ├── 6_shap_summary_cv_averaged.png/pdf
            ├── 8a_shap_waterfall_escape_example.png/pdf
            ├── 8b_shap_waterfall_sensitive_example.png/pdf
            ├── shap_feature_importance_rankings.csv
            └── feature_importance_comparison.csv
```

**Key operations:**
- 5-fold stratified CV with out-of-fold (OOF) predictions for unbiased metrics
- Native CatBoost feature importances averaged across all 5 folds
- SHAP values computed on 500-sample subsets per fold then aggregated (CV-averaged OOF)
- Optimal classification threshold via Youden's J index (~0.38–0.41)
- Final model trained on complete dataset for deployment/prediction

---

## Model Design

| Choice | Rationale |
|--------|-----------|
| CatBoost classifier | Native categorical feature support; robust to missing values |
| 5-fold stratified CV | Unbiased AUC; consistent with published NMD predictor benchmarks |
| OOF predictions | Unbiased probability estimates for threshold selection and plotting |
| Mean \|SHAP\| for importance | Fairer than native importance for correlated feature sets |
| Youden's J threshold | Balances sensitivity and specificity; preferred over default 0.5 |
| Final model on all data | Maximises training signal for genome-wide deployment |

### Top predictive features (by mean |SHAP|)

1. `last.EJC` — whether the PTC is in the last exon (canonical NMD rule)
2. `relativePTClocation` — normalized PTC position within the CDS
3. `half_life_PC1` — mRNA stability principal component
4. `cdsseqs_AU_content` — AU nucleotide content of the CDS
5. `phylop_ptc_to_ejc_median` — evolutionary conservation in the PTC–EJC window
6. `MedianExpression_log2` — log-transformed median expression (read-depth leakage control)
7. `mut.exon` — exon number of the PTC-bearing exon

---

## Citation

If you use this pipeline, please cite: *(manuscript in preparation)*

---

## Contact

Jacob Schmidt — computational biology, NMD prediction pipeline
