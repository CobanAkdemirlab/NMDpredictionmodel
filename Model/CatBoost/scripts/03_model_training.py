#!/usr/bin/env python3
"""
Model Training Script
=====================
This script trains a CatBoost classifier for NMD escape prediction:
- 5-fold stratified cross-validation
- Feature importance analysis
- SHAP analysis for interpretability
- Performance visualizations
- Model persistence

Requirements:
    - Python 3.12+
    - pandas
    - numpy
    - catboost
    - scikit-learn
    - matplotlib
    - seaborn
    - shap
    - pyyaml
    - joblib
    - tqdm

Usage:
    python 03_model_training.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import joblib
import json
import warnings
import sys

warnings.filterwarnings('ignore')

from catboost import CatBoostClassifier, Pool
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    roc_auc_score, 
    roc_curve,
    auc,
    precision_recall_curve,
    confusion_matrix,
    classification_report,
    precision_score,
    recall_score,
    f1_score
)
from tqdm import tqdm
import yaml
import shap


def load_config(config_path="../config/config.yaml"):
    """Load configuration from YAML file."""
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_data(config):
    """Load cleaned dataset."""
    print("=" * 80)
    print("LOADING DATA")
    print("=" * 80)
    
    df = pd.read_csv(config['data']['cleaned'])
    print(f"\n✓ Loaded: {df.shape}")
    
    target = config['model']['target']
    y = df[target].astype(int)
    X = df.drop(columns=[target])
    
    print(f"\nDataset Info:")
    print(f"  Samples: {len(y)}")
    print(f"  Features: {X.shape[1]}")
    print(f"  Escapees: {y.sum()} ({y.mean()*100:.1f}%)")
    print(f"  NMD: {(~y.astype(bool)).sum()} ({(~y.astype(bool)).sum()/len(y)*100:.1f}%)")
    
    return X, y


def identify_categorical_features(X, config):
    """Identify categorical features and their indices."""
    print("\n" + "=" * 80)
    print("IDENTIFYING CATEGORICAL FEATURES")
    print("=" * 80)
    
    declared_cat = config['features']['categorical']
    print(f"Categorical features from config: {len(declared_cat)}")
    
    cat_features = [c for c in declared_cat if c in X.columns]
    other_objs = [c for c in X.columns if X[c].dtype == "object" and c not in cat_features]
    cat_features.extend(other_objs)
    
    cat_indices = [X.columns.get_loc(c) for c in cat_features]
    
    print(f"\nCategorical features: {len(cat_features)}")
    print(f"Numeric features: {X.shape[1] - len(cat_features)}")
    
    return cat_features, cat_indices


def cross_validation_training(X, y, cat_indices, config):
    """Perform 5-fold cross-validation training."""
    print("\n" + "=" * 80)
    print("5-FOLD CROSS-VALIDATION")
    print("=" * 80)
    
    n_folds = config['model']['n_folds']
    random_seed = config['model']['random_seed']
    catboost_params = config['model']['catboost'].copy()
    
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=random_seed)
    cv_preds = np.zeros(len(y))
    cv_models = []
    fold_aucs = []
    feature_importances_folds = []
    
    best_fold_idx = None
    best_fold_auc = -np.inf
    
    print(f"\nTraining {n_folds} models...\n")
    
    for fold, (train_idx, val_idx) in enumerate(tqdm(skf.split(X, y), total=n_folds, desc="CV Folds")):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]
        
        train_pool = Pool(X_train, y_train, cat_features=cat_indices)
        val_pool = Pool(X_val, y_val, cat_features=cat_indices)
        
        model = CatBoostClassifier(
            random_seed=random_seed,
            cat_features=cat_indices,
            **catboost_params
        )
        model.fit(train_pool, eval_set=val_pool, verbose=False)
        
        val_preds = model.predict_proba(val_pool)[:, 1]
        cv_preds[val_idx] = val_preds
        
        auc_score = roc_auc_score(y_val, val_preds)
        fold_aucs.append(auc_score)
        
        cv_models.append(model)
        feature_importances_folds.append(model.feature_importances_)
        
        if auc_score > best_fold_auc:
            best_fold_auc = auc_score
            best_fold_idx = fold
        
        print(f"  Fold {fold+1}: AUC = {auc_score:.4f}")
    
    cv_auc = roc_auc_score(y, cv_preds)
    mean_auc = np.mean(fold_aucs)
    std_auc = np.std(fold_aucs)
    
    print("\n" + "=" * 80)
    print("CROSS-VALIDATION RESULTS")
    print("=" * 80)
    print(f"\nOut-of-Fold AUC: {cv_auc:.4f}")
    print(f"Mean Fold AUC: {mean_auc:.4f} ± {std_auc:.4f}")
    print(f"\nFold AUCs: {[f'{a:.4f}' for a in fold_aucs]}")
    print(f"Best Fold: {best_fold_idx+1} (AUC = {best_fold_auc:.4f})")
    
    return cv_models, cv_preds, fold_aucs, feature_importances_folds, best_fold_idx


def analyze_feature_importance(X, feature_importances_folds, cat_features, config):
    """Analyze and save feature importance."""
    print("\n" + "=" * 80)
    print("FEATURE IMPORTANCE (AVERAGED ACROSS FOLDS)")
    print("=" * 80)
    
    avg_importances = np.mean(feature_importances_folds, axis=0)
    std_importances = np.std(feature_importances_folds, axis=0)
    
    importance_df = pd.DataFrame({
        'feature': X.columns,
        'importance_mean': avg_importances,
        'importance_std': std_importances,
        'is_categorical': [col in cat_features for col in X.columns]
    }).sort_values('importance_mean', ascending=False)
    
    print("\nTop 20 Most Important Features:")
    print("=" * 80)
    for idx, row in importance_df.head(20).iterrows():
        cat_marker = "[CAT]" if row['is_categorical'] else "[NUM]"
        print(f"  {row['feature']:50s} {cat_marker:6s} {row['importance_mean']:8.4f} ± {row['importance_std']:.4f}")
    
    results_dir = Path(config['output']['results_dir'])
    importance_path = results_dir / "feature_importances_cv_averaged.csv"
    importance_path.parent.mkdir(parents=True, exist_ok=True)
    importance_df.to_csv(importance_path, index=False)
    print(f"\n✓ Saved: {importance_path}")
    
    return importance_df


def train_final_model(X, y, cat_indices, config):
    """Train final model on all data."""
    print("\n" + "=" * 80)
    print("TRAINING FINAL MODEL ON ALL DATA")
    print("=" * 80)
    
    random_seed = config['model']['random_seed']
    catboost_params = config['model']['catboost'].copy()
    
    full_pool = Pool(X, y, cat_features=cat_indices)
    
    final_model = CatBoostClassifier(
        random_seed=random_seed,
        cat_features=cat_indices,
        **catboost_params
    )
    
    print("\nTraining...")
    final_model.fit(full_pool, verbose=100)
    
    print(f"\n✓ Final model trained on all {len(y)} samples")
    
    return final_model


def save_models(final_model, cv_models, cv_preds, y, fold_aucs, best_fold_idx, config):
    """Save all models and results."""
    print("\n" + "=" * 80)
    print("SAVING MODELS AND RESULTS")
    print("=" * 80)
    
    n_folds = config['model']['n_folds']
    
    # Save final model
    print("\n📦 Saving final model (trained on all data)...")
    final_model_cbm = Path(config['output']['final_model_cbm'])
    final_model_pkl = Path(config['output']['final_model_pkl'])
    final_model_cbm.parent.mkdir(parents=True, exist_ok=True)
    
    final_model.save_model(str(final_model_cbm))
    joblib.dump(final_model, final_model_pkl)
    
    print(f"✓ Final model saved:")
    print(f"  {final_model_cbm}")
    print(f"  {final_model_pkl}")
    
    # Save CV fold models
    print(f"\n📦 Saving all {n_folds} CV fold models...")
    cv_models_dir = Path(config['output']['cv_models_dir'])
    cv_models_dir.mkdir(parents=True, exist_ok=True)
    
    for fold_idx, model in enumerate(cv_models):
        fold_num = fold_idx + 1
        auc = fold_aucs[fold_idx]
        is_best = (fold_idx == best_fold_idx)
        
        cbm_path = cv_models_dir / f"fold_{fold_num}_auc_{auc:.4f}.cbm"
        pkl_path = cv_models_dir / f"fold_{fold_num}_auc_{auc:.4f}.pkl"
        
        model.save_model(str(cbm_path))
        joblib.dump(model, pkl_path)
        
        best_marker = " ← BEST" if is_best else ""
        print(f"  Fold {fold_num} (AUC={auc:.4f}){best_marker}")
    
    # Save CV predictions
    print("\n💾 Saving cross-validation predictions...")
    cv_pred_df = pd.DataFrame({
        'true_label': y,
        'predicted_prob': cv_preds,
        'predicted_class': (cv_preds > 0.5).astype(int)
    })
    
    cv_pred_path = Path(config['output']['cv_predictions'])
    cv_pred_path.parent.mkdir(parents=True, exist_ok=True)
    cv_pred_df.to_csv(cv_pred_path, index=False)
    print(f"✓ CV predictions saved: {cv_pred_path}")
    
    # Save performance summary
    print("\n📊 Saving performance summary...")
    cv_auc = roc_auc_score(y, cv_preds)
    mean_auc = np.mean(fold_aucs)
    std_auc = np.std(fold_aucs)
    
    summary_path = Path(config['output']['results_dir']) / "model_performance_summary.txt"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(summary_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("NMD ESCAPEE MODEL - PERFORMANCE SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Dataset Info:\n")
        f.write(f"  Samples: {len(y)}\n")
        f.write(f"  Escapee Rate: {y.mean()*100:.1f}%\n\n")
        f.write(f"Cross-Validation Performance:\n")
        f.write(f"  Out-of-Fold AUC: {cv_auc:.4f}\n")
        f.write(f"  Mean Fold AUC: {mean_auc:.4f} ± {std_auc:.4f}\n")
        f.write(f"  Fold AUCs: {[f'{a:.4f}' for a in fold_aucs]}\n")
        f.write(f"  Best Fold: {best_fold_idx+1} (AUC={fold_aucs[best_fold_idx]:.4f})\n")
    
    print(f"✓ Summary saved: {summary_path}")


def create_visualizations(X, y, cv_preds, cv_models, cat_indices, fold_aucs, 
                         best_fold_idx, importance_df, config):
    """Create comprehensive visualizations."""
    print("\n" + "=" * 80)
    print("CREATING PUBLICATION-QUALITY VISUALIZATIONS")
    print("=" * 80)
    
    viz_dir = Path(config['output']['figures_dir']) / "visualizations_cv"
    shap_dir = viz_dir / "shap_manuscript"
    viz_dir.mkdir(parents=True, exist_ok=True)
    shap_dir.mkdir(parents=True, exist_ok=True)
    
    fig_dpi = config['visualization']['figure_dpi']
    
    # Compute optimal threshold
    fpr_tmp, tpr_tmp, thr_tmp = roc_curve(y, cv_preds)
    youden_index = tpr_tmp - fpr_tmp
    optimal_idx = np.argmax(youden_index)
    threshold = thr_tmp[optimal_idx]
    print(f"\nOptimal threshold (Youden): {threshold:.3f}")
    
    # 1. ROC Curve
    print("\n1. ROC Curve...")
    fpr, tpr, _ = roc_curve(y, cv_preds)
    roc_auc_val = auc(fpr, tpr)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(fpr, tpr, color='#2E86AB', lw=3, label=f'ROC (AUC = {roc_auc_val:.4f})')
    ax.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--', label='Random')
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate', fontsize=14, fontweight='bold')
    ax.set_ylabel('True Positive Rate', fontsize=14, fontweight='bold')
    ax.set_title('ROC Curve - OOF', fontsize=16, fontweight='bold', pad=20)
    ax.legend(loc="lower right", fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(viz_dir / "1_roc_curve.png", dpi=fig_dpi, bbox_inches='tight')
    plt.close()
    print("✓ Saved: 1_roc_curve.png")
    
    # 2. Precision-Recall Curve
    print("2. Precision-Recall Curve...")
    precision, recall, _ = precision_recall_curve(y, cv_preds)
    pr_auc_val = auc(recall, precision)
    baseline = y.mean()
    
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(recall, precision, color='#A23B72', lw=3, label=f'PR (AUC = {pr_auc_val:.4f})')
    ax.plot([0, 1], [baseline, baseline], color='gray', lw=2, linestyle='--',
            label=f'Baseline ({baseline:.3f})')
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('Recall', fontsize=14, fontweight='bold')
    ax.set_ylabel('Precision', fontsize=14, fontweight='bold')
    ax.set_title('Precision-Recall Curve - OOF', fontsize=16, fontweight='bold', pad=20)
    ax.legend(loc="upper right", fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(viz_dir / "2_precision_recall_curve.png", dpi=fig_dpi, bbox_inches='tight')
    plt.close()
    print("✓ Saved: 2_precision_recall_curve.png")
    
    # 3. Confusion Matrices
    print("3. Confusion Matrices...")
    oof_labels = (cv_preds >= threshold).astype(int)
    cm = confusion_matrix(y, oof_labels)
    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax1,
                cbar_kws={'label': 'Count'},
                annot_kws={'size': 16, 'weight': 'bold'})
    ax1.set_xlabel('Predicted', fontsize=14, fontweight='bold')
    ax1.set_ylabel('True', fontsize=14, fontweight='bold')
    ax1.set_title(f'Confusion Matrix (Counts) @ {threshold:.2f}',
                 fontsize=16, fontweight='bold', pad=20)
    ax1.set_xticklabels(['Sensitive', 'Escape'], fontsize=12)
    ax1.set_yticklabels(['Sensitive', 'Escape'], fontsize=12, rotation=0)
    
    sns.heatmap(cm_norm, annot=True, fmt='.2%', cmap='Oranges', ax=ax2,
                cbar_kws={'label': 'Percentage'},
                annot_kws={'size': 16, 'weight': 'bold'})
    ax2.set_xlabel('Predicted', fontsize=14, fontweight='bold')
    ax2.set_ylabel('True', fontsize=14, fontweight='bold')
    ax2.set_title(f'Confusion Matrix (Normalized) @ {threshold:.2f}',
                 fontsize=16, fontweight='bold', pad=20)
    ax2.set_xticklabels(['Sensitive', 'Escape'], fontsize=12)
    ax2.set_yticklabels(['Sensitive', 'Escape'], fontsize=12, rotation=0)
    
    plt.tight_layout()
    plt.savefig(viz_dir / "3_confusion_matrix.png", dpi=fig_dpi, bbox_inches='tight')
    plt.close()
    print("✓ Saved: 3_confusion_matrix.png")
    
    # 4. Feature Importance
    print("4. Feature Importances...")
    top_20 = importance_df.head(20)
    
    fig, ax = plt.subplots(figsize=(12, 10))
    bars = ax.barh(range(len(top_20)), top_20['importance_mean'].values,
                   xerr=top_20['importance_std'].values,
                   color='#2E86AB', edgecolor='black', linewidth=0.5,
                   capsize=3)
    ax.set_yticks(range(len(top_20)))
    ax.set_yticklabels(top_20['feature'].values, fontsize=11)
    ax.set_xlabel('Importance Score (Mean ± Std)', fontsize=14, fontweight='bold')
    ax.set_title('Top 20 Feature Importances (CV-Averaged)',
                 fontsize=16, fontweight='bold', pad=20)
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(viz_dir / "4_feature_importance_top20.png", dpi=fig_dpi, bbox_inches='tight')
    plt.close()
    print("✓ Saved: 4_feature_importance_top20.png")
    
    # 5. Probability Distribution
    print("5. Probability Distribution...")
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.hist(cv_preds[y == 0], bins=50, alpha=0.6,
            label='NMD Sensitive (True)', color='#2E86AB', edgecolor='black')
    ax.hist(cv_preds[y == 1], bins=50, alpha=0.6,
            label='NMD Escape (True)', color='#C73E1D', edgecolor='black')
    ax.axvline(x=threshold, color='black', linestyle='--', linewidth=2,
               label=f'Decision Threshold ({threshold:.2f})')
    ax.set_xlabel('Predicted Probability', fontsize=14, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=14, fontweight='bold')
    ax.set_title('Distribution of OOF Predicted Probabilities by True Class',
                 fontsize=16, fontweight='bold', pad=20)
    ax.legend(fontsize=12)
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(viz_dir / "5_probability_distribution.png", dpi=fig_dpi, bbox_inches='tight')
    plt.close()
    print("✓ Saved: 5_probability_distribution.png")
    
    # SHAP Analysis
    print("\n" + "=" * 80)
    print("SHAP ANALYSIS (CV-AVERAGED)")
    print("=" * 80)
    
    perform_shap_analysis(X, y, cv_models, cat_indices, importance_df, 
                         best_fold_idx, shap_dir, fig_dpi, config)
    
    # Save summary metrics
    prec = precision_score(y, oof_labels)
    rec = recall_score(y, oof_labels)
    f1 = f1_score(y, oof_labels)
    
    summary = {
        "oof_auc": float(roc_auc_val),
        "oof_pr_auc": float(pr_auc_val),
        "baseline_rate": float(baseline),
        "threshold_used": float(threshold),
        "precision_at_threshold": float(prec),
        "recall_at_threshold": float(rec),
        "f1_at_threshold": float(f1),
    }
    
    with open(viz_dir / "viz_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    
    print("\n✓ All visualizations complete!")


def perform_shap_analysis(X, y, cv_models, cat_indices, importance_df, 
                         best_fold_idx, shap_dir, fig_dpi, config):
    """Perform SHAP analysis on CV models."""
    print("\nComputing SHAP values across folds...")
    
    n_folds = config['model']['n_folds']
    random_seed = config['model']['random_seed']
    
    # Reconstruct CV splits
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=random_seed)
    fold_test_indices = [va for tr, va in skf.split(X, y)]
    
    all_shap_values = []
    all_samples = []
    all_base_values = []
    shap_sample_size = 500
    
    for fold_idx in range(n_folds):
        print(f"  Fold {fold_idx + 1}/{n_folds}...")
        
        fold_model = cv_models[fold_idx]
        test_idx = fold_test_indices[fold_idx]
        X_test_fold = X.iloc[test_idx]
        
        if len(X_test_fold) > shap_sample_size:
            X_test_fold = X_test_fold.sample(n=shap_sample_size, random_state=42)
        
        explainer = shap.TreeExplainer(fold_model)
        shap_vals = explainer.shap_values(X_test_fold)
        
        if isinstance(shap_vals, list):
            shap_vals = shap_vals[1]
            base_val = explainer.expected_value[1]
        else:
            base_val = explainer.expected_value
        
        all_shap_values.append(shap_vals)
        all_samples.append(X_test_fold)
        all_base_values.append(base_val)
    
    # Concatenate all OOF SHAP values
    shap_values_oof = np.vstack(all_shap_values)
    X_shap_oof = pd.concat(all_samples, ignore_index=True)
    base_value_avg = np.mean(all_base_values)
    
    print(f"\nTotal SHAP samples: {len(X_shap_oof)}")
    
    # SHAP Summary Plot
    print("\nGenerating SHAP summary plot...")
    top_n_features = config['visualization'].get('top_n_features', 15)
    top_features = importance_df.head(top_n_features)['feature'].tolist()
    top_features = [f for f in top_features if f in X_shap_oof.columns]
    top_idx = [X_shap_oof.columns.get_loc(f) for f in top_features]
    
    plt.figure(figsize=(12, 10))
    shap.summary_plot(
        shap_values_oof[:, top_idx],
        X_shap_oof.iloc[:, top_idx],
        show=False,
        plot_size=(12, 10)
    )
    plt.title('SHAP Feature Importance - CV-Averaged (OOF)',
              fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('SHAP Value (Impact on Escape Prediction)',
               fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(shap_dir / "Fig_SHAP_summary_cv_averaged.png",
                dpi=fig_dpi, bbox_inches='tight')
    plt.savefig(shap_dir / "Fig_SHAP_summary_cv_averaged.pdf",
                bbox_inches='tight')
    plt.close()
    print("✓ Saved: Fig_SHAP_summary_cv_averaged.png/pdf")
    
    # Save SHAP importance rankings
    mean_abs_shap = np.abs(shap_values_oof).mean(axis=0)
    shap_importance_df = pd.DataFrame({
        'feature': X_shap_oof.columns,
        'mean_abs_shap': mean_abs_shap
    }).sort_values('mean_abs_shap', ascending=False)
    
    shap_importance_df.to_csv(shap_dir / "shap_feature_importance_rankings.csv", index=False)
    print("✓ Saved: shap_feature_importance_rankings.csv")


def main():
    """Main execution function."""
    print("=" * 80)
    print("MODEL TRAINING")
    print("=" * 80)
    
    try:
        # Load config and data
        config = load_config()
        X, y = load_data(config)
        
        # Create output directories
        for key in ['models_dir', 'cv_models_dir', 'results_dir', 'figures_dir']:
            Path(config['output'][key]).mkdir(parents=True, exist_ok=True)
        
        # Identify categorical features
        cat_features, cat_indices = identify_categorical_features(X, config)
        
        # Cross-validation training
        cv_models, cv_preds, fold_aucs, feature_importances_folds, best_fold_idx = \
            cross_validation_training(X, y, cat_indices, config)
        
        # Feature importance analysis
        importance_df = analyze_feature_importance(X, feature_importances_folds, 
                                                   cat_features, config)
        
        # Train final model
        final_model = train_final_model(X, y, cat_indices, config)
        
        # Save models
        save_models(final_model, cv_models, cv_preds, y, fold_aucs, best_fold_idx, config)
        
        # Create visualizations
        create_visualizations(X, y, cv_preds, cv_models, cat_indices, fold_aucs,
                            best_fold_idx, importance_df, config)
        
        print("\n" + "=" * 80)
        print("MODEL TRAINING COMPLETE! 🎉")
        print("=" * 80)
        
        cv_auc = roc_auc_score(y, cv_preds)
        mean_auc = np.mean(fold_aucs)
        std_auc = np.std(fold_aucs)
        
        print(f"""
Performance Summary:
  Out-of-Fold AUC: {cv_auc:.4f}
  Mean Fold AUC: {mean_auc:.4f} ± {std_auc:.4f}
  Best Fold: {best_fold_idx+1} (AUC = {fold_aucs[best_fold_idx]:.4f})

Models Saved:
  ✓ Final model (trained on all data)
  ✓ All {config['model']['n_folds']} CV fold models

Outputs Saved:
  ✓ Feature importances (CV-averaged)
  ✓ CV predictions
  ✓ Performance plots
  ✓ SHAP analysis
  ✓ Summary report

All outputs saved to: {config['output']['results_dir']}
        """)
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
