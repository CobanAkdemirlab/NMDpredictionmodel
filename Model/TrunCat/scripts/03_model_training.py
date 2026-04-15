#!/usr/bin/env python3
"""
03_model_training.py
====================
Trains the CatBoost NMD escape prediction model using 5-fold stratified CV,
generates SHAP-based feature importances, and saves all outputs.

Steps:
  1. Load cleaned data (from script 02)
  2. Identify categorical features
  3. 5-fold stratified cross-validation
  4. Feature importance (CV-averaged, native + SHAP)
  5. Train final model on all data
  6. Save models and results
  7. Generate publication-quality visualizations

Output: results/models/, results/figures/, results/cv_predictions.csv

Usage:
    python 03_model_training.py
    python 03_model_training.py --config /path/to/config.yaml
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for script mode
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import joblib
import json
import warnings
import argparse
import sys
warnings.filterwarnings('ignore')

from catboost import CatBoostClassifier, Pool
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    roc_auc_score,
    roc_curve,
    auc,
    precision_recall_curve,
    average_precision_score,
    confusion_matrix,
    classification_report,
    precision_score,
    recall_score,
    f1_score
)
from tqdm import tqdm
import yaml
import shap


# ==============================================================================
# CONFIG
# ==============================================================================

def load_config(config_path=None):
    """Load config and resolve all paths relative to repo root."""
    if config_path is None:
        script_dir = Path(__file__).resolve().parent
        config_path = script_dir.parent / "config" / "config.yaml"
    else:
        config_path = Path(config_path).resolve()

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path) as f:
        config = yaml.safe_load(f)

    BASE_DIR = config_path.parent.parent
    for key, val in config['data'].items():
        config['data'][key] = str(BASE_DIR / val)
    for key, val in config['output'].items():
        if isinstance(val, str):
            config['output'][key] = str(BASE_DIR / val)

    print(f"✓ Configuration loaded from: {config_path}")
    return config


# ==============================================================================
# SETUP
# ==============================================================================

def setup(config):
    print("=" * 80)
    print("CONFIGURATION LOADED")
    print("=" * 80)

    PATH_INPUT = config['data']['cleaned']
    TARGET     = config['model']['target']
    RANDOM_SEED = config['model']['random_seed']
    N_FOLDS    = config['model']['n_folds']
    CATBOOST_PARAMS = config['model']['catboost'].copy()
    CATEGORICAL_FEATURES = config['features']['categorical']

    MODELS_DIR    = Path(config['output']['models_dir'])
    CV_MODELS_DIR = Path(config['output']['cv_models_dir'])
    RESULTS_DIR   = Path(config['output']['results_dir'])
    FIGURES_DIR   = Path(config['output']['figures_dir'])

    print(f"\nInput data: {PATH_INPUT}")
    print(f"Target: {TARGET}")
    print(f"Random seed: {RANDOM_SEED}")
    print(f"CV folds: {N_FOLDS}")
    print(f"\nCatBoost hyperparameters:")
    for key, value in CATBOOST_PARAMS.items():
        print(f"  {key}: {value}")
    print(f"\nOutput directories:")
    print(f"  Models:   {MODELS_DIR}")
    print(f"  CV folds: {CV_MODELS_DIR}")
    print(f"  Results:  {RESULTS_DIR}")
    print(f"  Figures:  {FIGURES_DIR}")
    print(f"\nCategorical features: {len(CATEGORICAL_FEATURES)}")

    for directory in [MODELS_DIR, CV_MODELS_DIR, RESULTS_DIR, FIGURES_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
    print("\n✓ All directories created")

    return (PATH_INPUT, TARGET, RANDOM_SEED, N_FOLDS, CATBOOST_PARAMS,
            CATEGORICAL_FEATURES, MODELS_DIR, CV_MODELS_DIR, RESULTS_DIR, FIGURES_DIR)


# ==============================================================================
# LOAD DATA
# ==============================================================================

def load_data(PATH_INPUT, TARGET, CATEGORICAL_FEATURES):
    print("\n" + "=" * 80)
    print("LOADING DATA")
    print("=" * 80)

    df = pd.read_csv(PATH_INPUT)
    print(f"\n✓ Loaded: {df.shape}")

    y = df[TARGET].astype(int)
    X = df.drop(columns=[TARGET])

    print(f"\nDataset Info:")
    print(f"  Samples: {len(y)}")
    print(f"  Features: {X.shape[1]}")
    print(f"  Escapees: {y.sum()} ({y.mean()*100:.1f}%)")
    print(f"  NMD: {(~y.astype(bool)).sum()} ({(~y.astype(bool)).sum()/len(y)*100:.1f}%)")

    # Identify categorical features
    print("\n" + "=" * 80)
    print("IDENTIFYING CATEGORICAL FEATURES")
    print("=" * 80)
    print(f"Categorical features from config: {len(CATEGORICAL_FEATURES)}")

    cat_features = [c for c in CATEGORICAL_FEATURES if c in X.columns]
    other_objs = [c for c in X.columns if X[c].dtype == "object" and c not in cat_features]
    cat_features.extend(other_objs)

    cat_indices = [X.columns.get_loc(c) for c in cat_features]

    print(f"\nCategorical features: {len(cat_features)}")
    print(f"Numeric features: {X.shape[1] - len(cat_features)}")
    print(f"\nCategorical feature names:")
    for cat in sorted(cat_features):
        print(f"  - {cat}")

    return X, y, cat_features, cat_indices


# ==============================================================================
# CROSS-VALIDATION
# ==============================================================================

def run_cross_validation(X, y, cat_indices, N_FOLDS, RANDOM_SEED, CATBOOST_PARAMS, CV_MODELS_DIR):
    print("\n" + "=" * 80)
    print(f"{N_FOLDS}-FOLD CROSS-VALIDATION")
    print("=" * 80)

    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_SEED)
    cv_preds = np.zeros(len(y))
    cv_models = []
    fold_aucs = []
    feature_importances_folds = []

    best_fold_idx = None
    best_fold_auc = -np.inf

    print(f"\nTraining {N_FOLDS} models...\n")

    for fold, (train_idx, val_idx) in enumerate(tqdm(skf.split(X, y), total=N_FOLDS, desc="CV Folds")):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

        train_pool = Pool(X_train, y_train, cat_features=cat_indices)
        val_pool   = Pool(X_val,   y_val,   cat_features=cat_indices)

        model = CatBoostClassifier(
            random_seed=RANDOM_SEED,
            cat_features=cat_indices,
            **CATBOOST_PARAMS
        )
        model.fit(train_pool, eval_set=val_pool, verbose=False)

        val_preds = model.predict_proba(val_pool)[:, 1]
        cv_preds[val_idx] = val_preds

        fold_auc = roc_auc_score(y_val, val_preds)
        fold_aucs.append(fold_auc)

        cv_models.append(model)
        feature_importances_folds.append(model.feature_importances_)

        if fold_auc > best_fold_auc:
            best_fold_auc = fold_auc
            best_fold_idx = fold

        print(f"  Fold {fold+1}: AUC = {fold_auc:.4f}")

    cv_auc  = roc_auc_score(y, cv_preds)
    mean_auc = np.mean(fold_aucs)
    std_auc  = np.std(fold_aucs)

    print(f"\nOut-of-Fold AUC: {cv_auc:.4f}")
    print(f"Mean Fold AUC: {mean_auc:.4f} ± {std_auc:.4f}")
    print(f"Fold AUCs: {[f'{a:.4f}' for a in fold_aucs]}")
    print(f"Best Fold: {best_fold_idx+1} (AUC = {best_fold_auc:.4f})")

    return cv_preds, cv_models, fold_aucs, feature_importances_folds, best_fold_idx, best_fold_auc, cv_auc, mean_auc, std_auc


# ==============================================================================
# FEATURE IMPORTANCE
# ==============================================================================

def compute_feature_importance(X, cat_features, feature_importances_folds, RESULTS_DIR):
    print("\n" + "=" * 80)
    print("FEATURE IMPORTANCE (AVERAGED ACROSS FOLDS)")
    print("=" * 80)

    avg_importances = np.mean(feature_importances_folds, axis=0)
    std_importances = np.std(feature_importances_folds,  axis=0)

    importance_df = pd.DataFrame({
        'feature': X.columns,
        'importance_mean': avg_importances,
        'importance_std':  std_importances,
        'is_categorical':  [col in cat_features for col in X.columns]
    }).sort_values('importance_mean', ascending=False)

    print("\nTop 20 Most Important Features:")
    print("=" * 80)
    for idx, row in importance_df.head(20).iterrows():
        cat_marker = "[CAT]" if row['is_categorical'] else "[NUM]"
        print(f"  {row['feature']:50s} {cat_marker:6s} {row['importance_mean']:8.4f} ± {row['importance_std']:.4f}")

    importance_path = Path(RESULTS_DIR) / "feature_importances_cv_averaged.csv"
    importance_df.to_csv(importance_path, index=False)
    print(f"\n✓ Saved: {importance_path}")

    return importance_df


# ==============================================================================
# FINAL MODEL
# ==============================================================================

def train_final_model(X, y, cat_indices, RANDOM_SEED, CATBOOST_PARAMS):
    print("\n" + "=" * 80)
    print("TRAINING FINAL MODEL ON ALL DATA")
    print("=" * 80)

    full_pool = Pool(X, y, cat_features=cat_indices)
    final_model = CatBoostClassifier(
        random_seed=RANDOM_SEED,
        cat_features=cat_indices,
        **CATBOOST_PARAMS
    )
    print("\nTraining...")
    final_model.fit(full_pool, verbose=100)
    print(f"\n✓ Final model trained on all {len(y)} samples")

    return final_model


# ==============================================================================
# SAVE MODELS AND RESULTS
# ==============================================================================

def save_models_and_results(
    final_model, cv_models, fold_aucs, best_fold_idx, best_fold_auc,
    cv_preds, y, X, cat_features, N_FOLDS, cv_auc, mean_auc, std_auc,
    CATBOOST_PARAMS, config, MODELS_DIR, CV_MODELS_DIR, RESULTS_DIR
):
    print("\n" + "=" * 80)
    print("SAVING MODELS AND RESULTS")
    print("=" * 80)

    # Final model
    print("\n📦 Saving final model (trained on all data)...")
    final_model_cbm = Path(config['output']['final_model_cbm'])
    final_model_pkl = Path(config['output']['final_model_pkl'])
    final_model_cbm.parent.mkdir(parents=True, exist_ok=True)
    final_model.save_model(str(final_model_cbm))
    joblib.dump(final_model, final_model_pkl)
    print(f"✓ Final model saved (trained on all {len(y)} samples):")
    print(f"  {final_model_cbm}")
    print(f"  {final_model_pkl}")

    # CV fold models
    print(f"\n📦 Saving all {N_FOLDS} CV fold models...")
    cv_models_dir = Path(CV_MODELS_DIR)
    cv_models_dir.mkdir(parents=True, exist_ok=True)
    for fold_idx, model in enumerate(cv_models):
        fold_num = fold_idx + 1
        fold_auc = fold_aucs[fold_idx]
        is_best  = (fold_idx == best_fold_idx)
        cbm_path = cv_models_dir / f"fold_{fold_num}_auc_{fold_auc:.4f}.cbm"
        pkl_path = cv_models_dir / f"fold_{fold_num}_auc_{fold_auc:.4f}.pkl"
        model.save_model(str(cbm_path))
        joblib.dump(model, pkl_path)
        best_marker = " ← BEST" if is_best else ""
        print(f"  Fold {fold_num} (AUC={fold_auc:.4f}){best_marker}")
    print(f"  All CV models saved to: {cv_models_dir}/")
    print(f"  Best fold: {best_fold_idx+1} (AUC={best_fold_auc:.4f})")

    # CV predictions
    print("\n💾 Saving cross-validation predictions...")
    cv_pred_df = pd.DataFrame({
        'true_label':      y,
        'predicted_prob':  cv_preds,
        'predicted_class': (cv_preds > 0.5).astype(int)
    })
    cv_pred_path = Path(config['output']['cv_predictions'])
    cv_pred_path.parent.mkdir(parents=True, exist_ok=True)
    cv_pred_df.to_csv(cv_pred_path, index=False)
    print(f"✓ CV predictions saved: {cv_pred_path}")

    # Performance summary
    print("\n📊 Saving performance summary...")
    summary = {
        'n_samples':     len(y),
        'n_features':    X.shape[1],
        'n_categorical': len(cat_features),
        'n_numeric':     X.shape[1] - len(cat_features),
        'escapee_rate':  y.mean(),
        'cv_auc_oof':    cv_auc,
        'cv_auc_mean':   mean_auc,
        'cv_auc_std':    std_auc,
        'best_fold':     best_fold_idx + 1,
        'best_fold_auc': best_fold_auc,
        'fold_aucs':     fold_aucs
    }
    summary_path = Path(RESULTS_DIR) / "model_performance_summary.txt"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("NMD ESCAPEE MODEL - PERFORMANCE SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        f.write("Dataset Info:\n")
        f.write(f"  Samples: {summary['n_samples']}\n")
        f.write(f"  Features: {summary['n_features']}\n")
        f.write(f"    - Categorical: {summary['n_categorical']}\n")
        f.write(f"    - Numeric: {summary['n_numeric']}\n")
        f.write(f"  Escapee Rate: {summary['escapee_rate']*100:.1f}%\n\n")
        f.write("Cross-Validation Performance:\n")
        f.write(f"  Out-of-Fold AUC: {summary['cv_auc_oof']:.4f}\n")
        f.write(f"  Mean Fold AUC: {summary['cv_auc_mean']:.4f} ± {summary['cv_auc_std']:.4f}\n")
        f.write(f"  Fold AUCs: {[f'{a:.4f}' for a in summary['fold_aucs']]}\n")
        f.write(f"  Best Fold: {summary['best_fold']} (AUC={summary['best_fold_auc']:.4f})\n\n")
        f.write("Models Saved:\n")
        f.write(f"  - Final model: {final_model_cbm.name}\n")
        f.write(f"  - All CV folds: {cv_models_dir.name}/fold_[1-{N_FOLDS}]_auc_*.cbm/.pkl\n")
        f.write(f"  - Best CV fold: Fold {summary['best_fold']} (AUC={summary['best_fold_auc']:.4f})\n\n")
        f.write("Hyperparameters (from config):\n")
        for key, value in CATBOOST_PARAMS.items():
            if key not in ['verbose']:
                f.write(f"  {key}: {value}\n")
    print(f"✓ Summary saved: {summary_path}")

    print("\n" + "=" * 80)
    print("✅ ALL MODELS AND RESULTS SAVED")
    print("=" * 80)
    print(f"\n📁 Output locations:")
    print(f"  Models:   {MODELS_DIR}")
    print(f"  CV folds: {CV_MODELS_DIR}")
    print(f"  Results:  {RESULTS_DIR}")
    print(f"\n🎯 Key files created:")
    print(f"  - {final_model_cbm}")
    print(f"  - {cv_pred_path}")
    print(f"  - {summary_path}")

    return summary, cv_pred_path


# ==============================================================================
# VISUALIZATIONS
# ==============================================================================

def generate_visualizations(
    X, y, cv_preds, cv_models, cat_indices, importance_df,
    N_FOLDS, RANDOM_SEED, config, FIGURES_DIR
):
    print("=" * 80)
    print("CREATING PUBLICATION-QUALITY VISUALIZATIONS")
    print("=" * 80)

    VIZ_DIR  = Path(FIGURES_DIR) / "visualizations_cv"
    SHAP_DIR = VIZ_DIR / "shap_manuscript"
    VIZ_DIR.mkdir(parents=True, exist_ok=True)
    SHAP_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {VIZ_DIR}")

    # Youden threshold
    fpr_tmp, tpr_tmp, thr_tmp = roc_curve(y, cv_preds)
    youden_index = tpr_tmp - fpr_tmp
    THRESHOLD = thr_tmp[np.argmax(youden_index)]
    print(f"\nOptimal threshold (Youden): {THRESHOLD:.3f}")

    fig_dpi = config['visualization']['figure_dpi']
    top_n_features = config['visualization'].get('top_n_features', 20)
    print(f"Figure DPI: {fig_dpi}")

    # ------------------------------------------------------------------
    # SHAP (runs first — needed for all downstream plots)
    # ------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("SHAP ANALYSIS (CV-AVERAGED)")
    print("=" * 80)

    print("\nReconstructing CV splits...")
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_SEED)
    fold_test_indices = [va for tr, va in skf.split(X, y)]

    print("Computing SHAP values across folds...")
    all_shap_values = []
    all_samples     = []
    all_base_values = []
    SHAP_SAMPLE_SIZE = 500

    for fold_idx in range(N_FOLDS):
        print(f"  Fold {fold_idx + 1}/{N_FOLDS}...")
        fold_model = cv_models[fold_idx]
        test_idx   = fold_test_indices[fold_idx]
        X_test_fold = X.iloc[test_idx]

        if len(X_test_fold) > SHAP_SAMPLE_SIZE:
            X_test_fold = X_test_fold.sample(n=SHAP_SAMPLE_SIZE, random_state=42)

        explainer = shap.TreeExplainer(fold_model)
        shap_vals = explainer.shap_values(X_test_fold)

        if isinstance(shap_vals, list):
            shap_vals = shap_vals[1]
            base_val  = explainer.expected_value[1]
        else:
            base_val  = explainer.expected_value

        all_shap_values.append(shap_vals)
        all_samples.append(X_test_fold)
        all_base_values.append(base_val)

    shap_values_oof = np.vstack(all_shap_values)
    X_shap_oof      = pd.concat(all_samples, ignore_index=True)
    base_value_avg  = np.mean(all_base_values)

    print(f"\nTotal SHAP samples: {len(X_shap_oof)}")
    print(f"Average base value: {base_value_avg:.4f}")

    mean_abs_shap = np.abs(shap_values_oof).mean(axis=0)
    shap_importance_df = pd.DataFrame({
        'feature':       X_shap_oof.columns,
        'mean_abs_shap': mean_abs_shap
    }).sort_values('mean_abs_shap', ascending=False)

    shap_importance_df.to_csv(SHAP_DIR / "shap_feature_importance_rankings.csv", index=False)
    print("✓ Saved: shap_feature_importance_rankings.csv")

    comparison_df = importance_df.merge(shap_importance_df, on='feature', how='outer').fillna(0)
    comparison_df['shap_rank']   = comparison_df['mean_abs_shap'].rank(ascending=False)
    comparison_df['native_rank'] = comparison_df['importance_mean'].rank(ascending=False)
    comparison_df.to_csv(SHAP_DIR / "feature_importance_comparison.csv", index=False)
    print("✓ Saved: feature_importance_comparison.csv")

    # ------------------------------------------------------------------
    # 1. ROC Curve
    # ------------------------------------------------------------------
    print("\n1. ROC Curve...")
    fpr, tpr, _ = roc_curve(y, cv_preds)
    roc_auc_val = auc(fpr, tpr)

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(fpr, tpr, color='#2E86AB', lw=3, label=f'ROC (AUC = {roc_auc_val:.4f})')
    ax.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--', label='Random')
    ax.set_xlim([0.0, 1.0]); ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate', fontsize=14, fontweight='bold')
    ax.set_ylabel('True Positive Rate', fontsize=14, fontweight='bold')
    ax.set_title('ROC Curve - OOF', fontsize=16, fontweight='bold', pad=20)
    ax.legend(loc="lower right", fontsize=12); ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(VIZ_DIR / "1_roc_curve.png", dpi=fig_dpi, bbox_inches='tight')
    plt.savefig(VIZ_DIR / "1_roc_curve.pdf", bbox_inches='tight')
    plt.close()
    print("✓ Saved: 1_roc_curve.png/pdf")

    # ------------------------------------------------------------------
    # 2. Precision-Recall Curve
    # ------------------------------------------------------------------
    print("2. Precision-Recall Curve...")
    precision, recall, _ = precision_recall_curve(y, cv_preds)
    pr_auc_val = average_precision_score(y, cv_preds)
    baseline   = y.mean()

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(recall, precision, color='#A23B72', lw=3, label=f'PR (AUC = {pr_auc_val:.4f})')
    ax.plot([0, 1], [baseline, baseline], color='gray', lw=2, linestyle='--',
            label=f'Baseline ({baseline:.3f})')
    ax.set_xlim([0.0, 1.0]); ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('Recall', fontsize=14, fontweight='bold')
    ax.set_ylabel('Precision', fontsize=14, fontweight='bold')
    ax.set_title('Precision-Recall Curve - OOF', fontsize=16, fontweight='bold', pad=20)
    ax.legend(loc="upper right", fontsize=12); ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(VIZ_DIR / "2_precision_recall_curve.png", dpi=fig_dpi, bbox_inches='tight')
    plt.savefig(VIZ_DIR / "2_precision_recall_curve.pdf", bbox_inches='tight')
    plt.close()
    print("✓ Saved: 2_precision_recall_curve.png/pdf")

    # ------------------------------------------------------------------
    # 3. Confusion Matrices
    # ------------------------------------------------------------------
    print("3. Confusion Matrices...")
    oof_labels = (cv_preds >= THRESHOLD).astype(int)
    cm      = confusion_matrix(y, oof_labels)
    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax1,
                cbar_kws={'label': 'Count'}, annot_kws={'size': 16, 'weight': 'bold'})
    ax1.set_xlabel('Predicted', fontsize=14, fontweight='bold')
    ax1.set_ylabel('True', fontsize=14, fontweight='bold')
    ax1.set_title(f'Confusion Matrix (Counts) @ {THRESHOLD:.2f}', fontsize=16, fontweight='bold', pad=20)
    ax1.set_xticklabels(['Sensitive', 'Escape'], fontsize=12)
    ax1.set_yticklabels(['Sensitive', 'Escape'], fontsize=12, rotation=0)

    sns.heatmap(cm_norm, annot=True, fmt='.2%', cmap='Oranges', ax=ax2,
                cbar_kws={'label': 'Percentage'}, annot_kws={'size': 16, 'weight': 'bold'})
    ax2.set_xlabel('Predicted', fontsize=14, fontweight='bold')
    ax2.set_ylabel('True', fontsize=14, fontweight='bold')
    ax2.set_title(f'Confusion Matrix (Normalized) @ {THRESHOLD:.2f}', fontsize=16, fontweight='bold', pad=20)
    ax2.set_xticklabels(['Sensitive', 'Escape'], fontsize=12)
    ax2.set_yticklabels(['Sensitive', 'Escape'], fontsize=12, rotation=0)

    plt.tight_layout()
    plt.savefig(VIZ_DIR / "3_confusion_matrix.png", dpi=fig_dpi, bbox_inches='tight')
    plt.savefig(VIZ_DIR / "3_confusion_matrix.pdf", bbox_inches='tight')
    plt.close()
    print("✓ Saved: 3_confusion_matrix.png/pdf")

    # ------------------------------------------------------------------
    # 4. Probability Distribution
    # ------------------------------------------------------------------
    print("4. Probability Distribution...")
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.hist(cv_preds[y == 0], bins=50, alpha=0.6, label='NMD Sensitive (True)', color='#2E86AB', edgecolor='black')
    ax.hist(cv_preds[y == 1], bins=50, alpha=0.6, label='NMD Escape (True)',    color='#C73E1D', edgecolor='black')
    ax.axvline(x=THRESHOLD, color='black', linestyle='--', linewidth=2,
               label=f'Decision Threshold ({THRESHOLD:.2f})')
    ax.set_xlabel('Predicted Probability', fontsize=14, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=14, fontweight='bold')
    ax.set_title('Distribution of OOF Predicted Probabilities by True Class',
                 fontsize=16, fontweight='bold', pad=20)
    ax.legend(fontsize=12); ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(VIZ_DIR / "4_probability_distribution.png", dpi=fig_dpi, bbox_inches='tight')
    plt.savefig(VIZ_DIR / "4_probability_distribution.pdf", bbox_inches='tight')
    plt.close()
    print("✓ Saved: 4_probability_distribution.png/pdf")

    # ------------------------------------------------------------------
    # 5. Performance Dashboard
    # ------------------------------------------------------------------
    print("5. Performance Dashboard...")
    fig = plt.figure(figsize=(16, 10))
    gs  = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    ax1 = fig.add_subplot(gs[0:2, 0])
    ax1.plot(fpr, tpr, color='#2E86AB', lw=2.5, label=f'AUC = {roc_auc_val:.4f}')
    ax1.plot([0, 1], [0, 1], 'k--', lw=1.5, alpha=0.5)
    ax1.set_xlabel('FPR', fontsize=10, fontweight='bold')
    ax1.set_ylabel('TPR', fontsize=10, fontweight='bold')
    ax1.set_title('ROC (OOF)', fontsize=12, fontweight='bold')
    ax1.legend(loc="lower right", fontsize=9); ax1.grid(True, alpha=0.3)

    ax2 = fig.add_subplot(gs[0:2, 1])
    ax2.plot(recall, precision, color='#A23B72', lw=2.5, label=f'AUC = {pr_auc_val:.4f}')
    ax2.axhline(y=baseline, color='k', linestyle='--', lw=1.5, alpha=0.5)
    ax2.set_xlabel('Recall', fontsize=10, fontweight='bold')
    ax2.set_ylabel('Precision', fontsize=10, fontweight='bold')
    ax2.set_title('Precision-Recall (OOF)', fontsize=12, fontweight='bold')
    ax2.legend(loc="upper right", fontsize=9); ax2.grid(True, alpha=0.3)

    ax3 = fig.add_subplot(gs[0:2, 2])
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax3, cbar=False,
                annot_kws={'size': 12, 'weight': 'bold'})
    ax3.set_xlabel('Predicted', fontsize=10, fontweight='bold')
    ax3.set_ylabel('True', fontsize=10, fontweight='bold')
    ax3.set_title(f'Confusion Matrix @ {THRESHOLD:.2f}', fontsize=12, fontweight='bold')
    ax3.set_xticklabels(['Sensitive', 'Escape'], fontsize=9)
    ax3.set_yticklabels(['Sensitive', 'Escape'], fontsize=9, rotation=0)

    ax4 = fig.add_subplot(gs[2, :])
    top_10_shap = shap_importance_df.head(10)
    ax4.barh(range(len(top_10_shap)), top_10_shap['mean_abs_shap'].values,
             color='#2E86AB', edgecolor='black', linewidth=0.5)
    ax4.set_yticks(range(len(top_10_shap)))
    ax4.set_yticklabels(top_10_shap['feature'].values, fontsize=9)
    ax4.set_xlabel('Mean |SHAP Value|', fontsize=10, fontweight='bold')
    ax4.set_title('Top 10 Features by Mean |SHAP|', fontsize=12, fontweight='bold')
    ax4.invert_yaxis(); ax4.grid(axis='x', alpha=0.3)

    fig.suptitle('NMD Escapee Model - CV/OOF Performance Summary',
                 fontsize=18, fontweight='bold', y=0.98)
    plt.savefig(VIZ_DIR / "5_performance_dashboard.png", dpi=fig_dpi, bbox_inches='tight')
    plt.savefig(VIZ_DIR / "5_performance_dashboard.pdf", bbox_inches='tight')
    plt.close()
    print("✓ Saved: 5_performance_dashboard.png/pdf")

    # ------------------------------------------------------------------
    # 6. SHAP Summary Plot
    # ------------------------------------------------------------------
    print("\n6. SHAP Summary Plot...")
    top_features = shap_importance_df.head(top_n_features)['feature'].tolist()
    top_features = [f for f in top_features if f in X_shap_oof.columns]
    top_idx = [X_shap_oof.columns.get_loc(f) for f in top_features]

    plt.figure(figsize=(12, 10))
    shap.summary_plot(shap_values_oof[:, top_idx], X_shap_oof.iloc[:, top_idx],
                      show=False, plot_size=(12, 10))
    plt.title('SHAP Feature Importance - CV-Averaged (OOF)', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('SHAP Value (Impact on Escape Prediction)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(SHAP_DIR / "6_shap_summary_cv_averaged.png", dpi=fig_dpi, bbox_inches='tight')
    plt.savefig(SHAP_DIR / "6_shap_summary_cv_averaged.pdf", bbox_inches='tight')
    plt.close()
    print("✓ Saved: 6_shap_summary_cv_averaged.png/pdf")

    # ------------------------------------------------------------------
    # 7a. Top 20 by Mean |SHAP|
    # ------------------------------------------------------------------
    print("7a. Feature Importances — Mean |SHAP| (primary)...")
    top_20_shap = shap_importance_df.head(20)
    fig, ax = plt.subplots(figsize=(12, 10))
    bars = ax.barh(range(len(top_20_shap)), top_20_shap['mean_abs_shap'].values,
                   color='#2E86AB', edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(top_20_shap)))
    ax.set_yticklabels(top_20_shap['feature'].values, fontsize=11)
    ax.set_xlabel('Mean |SHAP Value| (OOF)', fontsize=14, fontweight='bold')
    ax.set_title('Top 20 Feature Importances by Mean |SHAP| (CV-Averaged OOF)',
                 fontsize=16, fontweight='bold', pad=20)
    ax.invert_yaxis(); ax.grid(axis='x', alpha=0.3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    for bar, val in zip(bars, top_20_shap['mean_abs_shap'].values):
        ax.text(val + 0.002, bar.get_y() + bar.get_height()/2,
                f'{val:.3f}', va='center', fontsize=9, fontweight='bold')
    plt.tight_layout()
    plt.savefig(VIZ_DIR / "7a_feature_importance_shap_top20.png", dpi=fig_dpi, bbox_inches='tight')
    plt.savefig(VIZ_DIR / "7a_feature_importance_shap_top20.pdf", bbox_inches='tight')
    plt.close()
    print("✓ Saved: 7a_feature_importance_shap_top20.png/pdf")

    # ------------------------------------------------------------------
    # 7b. Top 20 Native CatBoost importance
    # ------------------------------------------------------------------
    print("7b. Feature Importances — Native CatBoost (supplementary)...")
    top_20_native = importance_df.head(20)
    fig, ax = plt.subplots(figsize=(12, 10))
    bars = ax.barh(range(len(top_20_native)), top_20_native['importance_mean'].values,
                   xerr=top_20_native['importance_std'].values,
                   color='#5C6BC0', edgecolor='black', linewidth=0.5, capsize=3)
    ax.set_yticks(range(len(top_20_native)))
    ax.set_yticklabels(top_20_native['feature'].values, fontsize=11)
    ax.set_xlabel('Importance Score (Mean ± Std)', fontsize=14, fontweight='bold')
    ax.set_title('Top 20 Feature Importances — Native CatBoost (CV-Averaged)',
                 fontsize=16, fontweight='bold', pad=20)
    ax.invert_yaxis(); ax.grid(axis='x', alpha=0.3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    for bar, val, std in zip(bars, top_20_native['importance_mean'].values, top_20_native['importance_std'].values):
        ax.text(val + std + 0.15, bar.get_y() + bar.get_height()/2,
                f'{val:.2f}±{std:.2f}', va='center', fontsize=9, fontweight='bold')
    plt.tight_layout()
    plt.savefig(VIZ_DIR / "7b_feature_importance_native_top20.png", dpi=fig_dpi, bbox_inches='tight')
    plt.savefig(VIZ_DIR / "7b_feature_importance_native_top20.pdf", bbox_inches='tight')
    plt.close()
    print("✓ Saved: 7b_feature_importance_native_top20.png/pdf")

    # ------------------------------------------------------------------
    # 7c. Rank comparison: SHAP vs Native (normalized)
    # ------------------------------------------------------------------
    print("7c. Feature Importance Rank Comparison...")
    top_20_comp  = comparison_df.nlargest(20, 'mean_abs_shap').copy()
    top_20_comp  = top_20_comp.sort_values('shap_rank')
    shap_norm    = top_20_comp['mean_abs_shap'].values / top_20_comp['mean_abs_shap'].max()
    native_norm  = top_20_comp['importance_mean'].values / top_20_comp['importance_mean'].max()

    fig, ax = plt.subplots(figsize=(12, 10))
    y_pos = np.arange(len(top_20_comp))
    ax.barh(y_pos,       shap_norm,   color='#2E86AB', alpha=0.8, label='Mean |SHAP| (normalized)', height=0.4)
    ax.barh(y_pos + 0.4, native_norm, color='#5C6BC0', alpha=0.8, label='Native Importance (normalized)', height=0.4)
    ax.set_yticks(y_pos + 0.2)
    ax.set_yticklabels(top_20_comp['feature'].values, fontsize=11)
    ax.set_xlabel('Normalized Score (max = 1)', fontsize=14, fontweight='bold')
    ax.set_title('Feature Importance: Mean |SHAP| vs Native CatBoost\n(Top 20 by SHAP, both normalized to max = 1)',
                 fontsize=14, fontweight='bold', pad=20)
    ax.invert_yaxis()
    ax.legend(fontsize=11, loc='lower right'); ax.grid(axis='x', alpha=0.3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(VIZ_DIR / "7c_feature_importance_comparison.png", dpi=fig_dpi, bbox_inches='tight')
    plt.savefig(VIZ_DIR / "7c_feature_importance_comparison.pdf", bbox_inches='tight')
    plt.close()
    print("✓ Saved: 7c_feature_importance_comparison.png/pdf")

    # ------------------------------------------------------------------
    # 8. SHAP Waterfall Examples
    # ------------------------------------------------------------------
    print("\n8. SHAP Waterfall Examples...")
    best_fold_model    = cv_models[np.argmax([roc_auc_score(y.iloc[fold_test_indices[i]],
                                               cv_models[i].predict_proba(
                                                   Pool(X.iloc[fold_test_indices[i]], cat_features=cat_indices))[:, 1])
                                              for i in range(N_FOLDS)])]
    explainer_waterfall = shap.TreeExplainer(best_fold_model)
    probs_full = best_fold_model.predict_proba(Pool(X, cat_features=cat_indices))[:, 1]

    # Example 1: high-confidence escape
    escape_candidates = np.where((y.values == 1) & (probs_full > 0.85))[0]
    if len(escape_candidates) > 0:
        example_escape_idx  = int(escape_candidates[0])
        shap_vals_escape    = explainer_waterfall.shap_values(X.iloc[[example_escape_idx]])
        if isinstance(shap_vals_escape, list):
            shap_vals_escape = shap_vals_escape[1]
            base_val_escape  = explainer_waterfall.expected_value[1]
        else:
            base_val_escape  = explainer_waterfall.expected_value

        exp_escape = shap.Explanation(
            values=shap_vals_escape[0],
            base_values=base_val_escape if np.ndim(base_val_escape) == 0 else base_val_escape.ravel()[0],
            data=X.iloc[example_escape_idx].values,
            feature_names=X.columns.tolist()
        )
        plt.figure(figsize=(12, 8))
        shap.waterfall_plot(exp_escape, show=False, max_display=15)
        plt.title(f"SHAP Explanation - NMD Escape Example\n(True=Escape, Pred={probs_full[example_escape_idx]:.3f})",
                  fontsize=14, fontweight='bold', pad=20)
        plt.tight_layout()
        plt.savefig(SHAP_DIR / "8a_shap_waterfall_escape_example.png", dpi=fig_dpi, bbox_inches='tight')
        plt.savefig(SHAP_DIR / "8a_shap_waterfall_escape_example.pdf", bbox_inches='tight')
        plt.close()
        print("✓ Saved: 8a_shap_waterfall_escape_example.png/pdf")

    # Example 2: high-confidence sensitive
    sensitive_candidates = np.where((y.values == 0) & (probs_full < 0.15))[0]
    if len(sensitive_candidates) > 0:
        example_sensitive_idx = int(sensitive_candidates[0])
        shap_vals_sens        = explainer_waterfall.shap_values(X.iloc[[example_sensitive_idx]])
        if isinstance(shap_vals_sens, list):
            shap_vals_sens = shap_vals_sens[1]
            base_val_sens  = explainer_waterfall.expected_value[1]
        else:
            base_val_sens  = explainer_waterfall.expected_value

        exp_sens = shap.Explanation(
            values=shap_vals_sens[0],
            base_values=base_val_sens if np.ndim(base_val_sens) == 0 else base_val_sens.ravel()[0],
            data=X.iloc[example_sensitive_idx].values,
            feature_names=X.columns.tolist()
        )
        plt.figure(figsize=(12, 8))
        shap.waterfall_plot(exp_sens, show=False, max_display=15)
        plt.title(f"SHAP Explanation - NMD Sensitive Example\n(True=Sensitive, Pred={probs_full[example_sensitive_idx]:.3f})",
                  fontsize=14, fontweight='bold', pad=20)
        plt.tight_layout()
        plt.savefig(SHAP_DIR / "8b_shap_waterfall_sensitive_example.png", dpi=fig_dpi, bbox_inches='tight')
        plt.savefig(SHAP_DIR / "8b_shap_waterfall_sensitive_example.pdf", bbox_inches='tight')
        plt.close()
        print("✓ Saved: 8b_shap_waterfall_sensitive_example.png/pdf")

    # Summary metrics + json
    prec = precision_score(y, oof_labels)
    rec  = recall_score(y, oof_labels)
    f1   = f1_score(y, oof_labels)

    viz_summary = {
        "oof_auc":              float(roc_auc_val),
        "oof_pr_auc":           float(pr_auc_val),
        "baseline_rate":        float(baseline),
        "threshold_used":       float(THRESHOLD),
        "precision_at_threshold": float(prec),
        "recall_at_threshold":  float(rec),
        "f1_at_threshold":      float(f1),
        "n_folds_used":         N_FOLDS,
        "shap_samples_total":   len(X_shap_oof),
    }
    with open(VIZ_DIR / "viz_summary.json", "w") as f:
        json.dump(viz_summary, f, indent=2)

    print("\n" + "=" * 80)
    print("VISUALIZATIONS COMPLETE! 🎉")
    print("=" * 80)
    print(f"\nAll visualizations saved to: {VIZ_DIR}")
    print(f"SHAP visualizations saved to: {SHAP_DIR}")

    # Also save the config-path confusion matrix for compatibility
    cm_path = Path(config['output']['confusion_matrix'])
    cm_path.parent.mkdir(parents=True, exist_ok=True)
    y_pred_class = (cv_preds > THRESHOLD).astype(int)
    cm2 = confusion_matrix(y, y_pred_class)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(cm2, annot=True, fmt='d', cmap='Blues',
                xticklabels=['NMD', 'Escapee'],
                yticklabels=['NMD', 'Escapee'],
                cbar_kws={'label': 'Count'}, ax=ax)
    ax.set_xlabel('Predicted Label', fontsize=12)
    ax.set_ylabel('True Label', fontsize=12)
    ax.set_title(f'Confusion Matrix @ Threshold={THRESHOLD:.3f} (Out-of-Fold)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(cm_path, dpi=config['visualization']['figure_dpi'], bbox_inches='tight')
    plt.close()
    print(f"\n✓ Confusion matrix plot saved: {cm_path}")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="Model training for NMD escape prediction")
    parser.add_argument('--config', default=None, help='Path to config.yaml (default: ../config/config.yaml)')
    args = parser.parse_args()

    try:
        config = load_config(args.config)

        (PATH_INPUT, TARGET, RANDOM_SEED, N_FOLDS, CATBOOST_PARAMS,
         CATEGORICAL_FEATURES, MODELS_DIR, CV_MODELS_DIR, RESULTS_DIR, FIGURES_DIR) = setup(config)

        X, y, cat_features, cat_indices = load_data(PATH_INPUT, TARGET, CATEGORICAL_FEATURES)

        (cv_preds, cv_models, fold_aucs, feature_importances_folds,
         best_fold_idx, best_fold_auc, cv_auc, mean_auc, std_auc) = run_cross_validation(
            X, y, cat_indices, N_FOLDS, RANDOM_SEED, CATBOOST_PARAMS, CV_MODELS_DIR
        )

        importance_df = compute_feature_importance(X, cat_features, feature_importances_folds, RESULTS_DIR)

        final_model = train_final_model(X, y, cat_indices, RANDOM_SEED, CATBOOST_PARAMS)

        save_models_and_results(
            final_model, cv_models, fold_aucs, best_fold_idx, best_fold_auc,
            cv_preds, y, X, cat_features, N_FOLDS, cv_auc, mean_auc, std_auc,
            CATBOOST_PARAMS, config, MODELS_DIR, CV_MODELS_DIR, RESULTS_DIR
        )

        generate_visualizations(
            X, y, cv_preds, cv_models, cat_indices, importance_df,
            N_FOLDS, RANDOM_SEED, config, FIGURES_DIR
        )

        print("\n" + "=" * 80)
        print("MODEL TRAINING COMPLETE! 🎉")
        print("=" * 80)
        print(f"""
Performance Summary:
  Out-of-Fold AUC: {cv_auc:.4f}
  Mean Fold AUC: {mean_auc:.4f} ± {std_auc:.4f}
  Best Fold: {best_fold_idx+1} (AUC = {best_fold_auc:.4f})

Models Saved:
  ✓ Final model (trained on all data)
  ✓ All {N_FOLDS} CV fold models

Outputs Saved:
  ✓ Feature importances (CV-averaged, native + SHAP)
  ✓ CV predictions
  ✓ Performance plots
  ✓ Confusion matrix
  ✓ Summary report

All outputs saved to: {RESULTS_DIR}
""")
        return 0

    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
