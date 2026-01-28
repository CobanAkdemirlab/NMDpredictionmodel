# NMD Escape Prediction Model

Machine learning pipeline for predicting nonsense-mediated decay (NMD) escape in genetic variants with premature termination codons (PTCs).

## Overview

This repository contains the complete workflow for training and evaluating a CatBoost classifier that predicts whether stop-gain variants will escape NMD surveillance. The model incorporates positional information, RNA-binding protein motifs, conservation scores, and transcript characteristics to achieve high prediction accuracy.

## Repository Structure

```
NMDpredictionmodel/Model/CatBoost/
├── config/
│   └── config.yaml              # Shared configuration file
├── notebooks/
│   ├── 01_data_loading_and_merging.ipynb
│   ├── 02_feature_cleaning_and_selection.ipynb
│   └── 03_model_training.ipynb
└── scripts/                     # Python scripts (this folder)
    ├── 01_data_loading_and_merging.py
    ├── 02_feature_cleaning_and_selection.py
    ├── 03_model_training.py
    ├── requirements.txt
    ├── setup.py
    └── README.md
```

**Note:** These scripts use the same `config/config.yaml` file as the notebooks. No separate configuration needed!

## Requirements

- Python 3.12+
- See `requirements.txt` for full dependencies

### Installation

```bash
# Create a virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Configuration

These scripts use the **same configuration file** as the notebooks: `../config/config.yaml`

If you've already set up the notebooks, you're good to go! Otherwise, configure paths in `../config/config.yaml`:

```yaml
data:
  enhanced_features: "path/to/enhanced_features.csv"
  codon_optimality: "path/to/codon_optimality.tsv"
  halflife: "path/to/halflife_features.csv"
  readthrough: "path/to/readthrough_scores.csv"
  merged: "data/processed/merged_data.csv"
  cleaned: "data/processed/cleaned_data.csv"
  
model:
  target: "escapee"
  random_seed: 42
  n_folds: 5
  
  catboost:
    iterations: 1000
    learning_rate: 0.05
    depth: 6
    l2_leaf_reg: 3
    # ... additional parameters
```

## Usage

### Step 1: Data Loading and Merging

Loads multiple data sources, filters indel variants, and merges features:

```bash
python 01_data_loading_and_merging.py
```

**Inputs:**
- Enhanced features CSV
- Codon optimality TSV
- Half-life features CSV
- Readthrough scores CSV

**Output:**
- `data/processed/merged_data.csv`

**Key operations:**
- Identifies and removes indel variants
- Replaces outdated features (e.g., AverageCodonRNAUsage → CodonOptimalityFraction)
- Merges all data sources by variant key

### Step 2: Feature Cleaning and Selection

Applies domain knowledge-based feature engineering and quality checks:

```bash
python 02_feature_cleaning_and_selection.py
```

**Input:**
- `data/processed/merged_data.csv`

**Output:**
- `data/processed/cleaned_data.csv`
- `data/processed/feature_list.csv`

**Key operations:**
- Removes identifier and leaky features
- Applies biological domain rules:
  - Keeps AU content, drops GC content (AU more relevant)
  - Protects RNA-binding protein features from correlation removal
  - Keeps PTC.2.EJC over PTC_dist_exon_end_0b
- Filters highly correlated features (threshold: 0.95)
- Handles categorical features and imputes missing values

### Step 3: Model Training

Trains CatBoost classifier with 5-fold cross-validation:

```bash
python 03_model_training.py
```

**Input:**
- `data/processed/cleaned_data.csv`

**Outputs:**
- Final model trained on all data
- All 5 CV fold models
- Feature importances (CV-averaged)
- Out-of-fold predictions
- Performance visualizations
- SHAP analysis plots

**Key operations:**
- 5-fold stratified cross-validation
- Feature importance analysis
- SHAP interpretability analysis
- Comprehensive performance evaluation

## Model Performance

Expected performance metrics (approximate):
- **ROC-AUC:** ~0.75 (out-of-fold)
- **Optimal threshold:** ~0.38-0.41 (Youden's index)

## Output Files

### Models
```
output/models/
├── final_model.cbm           # Final model (CatBoost format)
├── final_model.pkl           # Final model (pickle format)
└── cv_folds/                 # Individual CV fold models
    ├── fold_1_auc_*.cbm
    ├── fold_2_auc_*.cbm
    └── ...
```

### Results
```
output/results/
├── cv_predictions.csv                    # Out-of-fold predictions
├── feature_importances_cv_averaged.csv   # Feature importance rankings
└── model_performance_summary.txt         # Performance summary
```

### Visualizations
```
output/figures/visualizations_cv/
├── 1_roc_curve.png
├── 2_precision_recall_curve.png
├── 3_confusion_matrix.png
├── 4_feature_importance_top20.png
├── 5_probability_distribution.png
└── shap_manuscript/
    ├── Fig_SHAP_summary_cv_averaged.png
    ├── Fig_SHAP_summary_cv_averaged.pdf
    └── shap_feature_importance_rankings.csv
```



### Model Design Principles
1. **CatBoost classifier:** Native handling of categorical features and missing values
2. **Cross-validation:** 5-fold stratified CV for unbiased performance estimates
3. **Feature protection:** Biologically meaningful features retained despite correlation
4. **Interpretability:** SHAP analysis for understanding predictions
