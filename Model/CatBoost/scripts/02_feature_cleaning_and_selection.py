#!/usr/bin/env python3
"""
Feature Cleaning and Selection Script
======================================
This script performs feature engineering and selection:
- Removes identifier and leaky features
- Applies domain-knowledge based rules (AU/GC content, RBP protection, PTC distance)
- Filters highly correlated features
- Handles categorical features
- Imputes missing values

Requirements:
    - Python 3.12+
    - pandas
    - numpy
    - pyyaml

Usage:
    python 02_feature_cleaning_and_selection.py
"""

import pandas as pd
import numpy as np
from pathlib import Path
import yaml
import sys


def load_config(config_path="../config/config.yaml"):
    """Load configuration from YAML file."""
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_data(config):
    """Load merged dataset."""
    print("=" * 80)
    print("FEATURE CLEANING AND SELECTION")
    print("=" * 80)
    
    print("\nLoading data...")
    df = pd.read_csv(config['data']['merged'])
    print(f"✓ Loaded: {df.shape}")
    
    # Separate target and features
    target = config['model']['target']
    y_full = df[target].map({"TRUE": 1, "FALSE": 0, True: 1, False: 0}).astype("float")
    mask = y_full.notna()
    y = y_full.loc[mask].astype(int).reset_index(drop=True)
    X = df.loc[mask].drop(columns=[target]).reset_index(drop=True)
    
    print(f"Valid samples: {len(y)}")
    print(f"  Escapees: {y.sum()} ({y.sum()/len(y)*100:.1f}%)")
    print(f"  NMD: {(~y.astype(bool)).sum()} ({(~y.astype(bool)).sum()/len(y)*100:.1f}%)")
    print(f"Features: {X.shape[1]}")
    
    return X, y


def apply_manual_drops(X, config):
    """Apply custom removal rules based on domain knowledge."""
    print("\n" + "=" * 80)
    print("CUSTOM REMOVAL RULES")
    print("=" * 80)
    
    # Rename typo
    if 'loefu_cat' in X.columns:
        X = X.rename(columns={'loefu_cat': 'loeuf_cat'})
        print("\n✓ Renamed 'loefu_cat' → 'loeuf_cat'")
    
    # Manual drops
    drop_manual = {
        "gene_length": "Identifier",
        "noncoding.length": "Identifier",
        "cds_length.x": "Duplicate",
        "cds_length.y": "Duplicate",
        "GENE_ID": "Identifier",
        "key": "Identifier",
        "txnames": "Identifier",
        "cds.nucloc": "Categorical indicator",
        "penultimate.exon": "Leaky - redundant with last.EJC",
        "last.exon": "Leaky - redundant with last.EJC",
        "downstream": "Redundant - less granular than last.EJC",
        "ALLELE.RAT": "Leaky - used in target calculation",
        "mRNAHalfLifeMin": "Redundant - half_life_PC1 kept",
        "mean_half_life": "Redundant - half_life_PC1 kept",
        "median_half_life": "Redundant - half_life_PC1 kept",
        "PTC.2.start.binning": "Redundant - unbinned version kept",
        "PTC.2.end.binning": "Redundant - unbinned version kept",
        "PTC.2.EJC.binning": "Redundant - unbinned version kept",
        "length.mutated.exon.binning": "Redundant - unbinned version kept",
        "EJC.downstream.cut": "Redundant - binned version not needed",
        "REVEL": "Not applicable - missense score for stopgain variants",
        "FiveUTR_AUcontentlast200": "Redundant - overall fiveUTR_AU_content kept",
        "FiveUTR_AUcontentfirst200": "Redundant - overall fiveUTR_AU_content kept",
        "FiveUTR_UCcontentlast200": "Redundant - overall fiveUTR_UC_content kept",
        "FiveUTR_UCcontentfirst200": "Redundant - overall fiveUTR_UC_content kept",
        "fiveUTR_AUcontentlast100": "Redundant - overall fiveUTR_AU_content kept",
        "fiveUTR_AUcontentfirst100": "Redundant - overall fiveUTR_AU_content kept",
        "fiveUTR_UCcontentlast100": "Redundant - overall fiveUTR_UC_content kept",
        "fiveUTR_UCcontentfirst100": "Redundant - overall fiveUTR_UC_content kept",
        "aug_distance_nt": "Redundant - aug_distance_category captures this information",
        "COHORT_AF": "Data leakage - computed from training data, not generalizable"
    }
    
    # Rule 1: GC content features
    print("\nRule 1: AU/GC content pairs")
    print("  → Keep AU content (more biologically relevant)")
    print("  → Drop GC content (redundant)")
    
    gc_features_to_drop = [col for col in X.columns if '_GC_content' in col or 'GCcontent' in col]
    print(f"  Dropping {len(gc_features_to_drop)} GC content features")
    
    for col in gc_features_to_drop:
        drop_manual[col] = "Redundant - AU content kept (more biologically relevant)"
    
    # Rule 2: RBP features
    print("\nRule 2: RBP features")
    print("  → Keep ALL RBP features (biologically meaningful correlation)")
    print("  → Will be protected from automated correlation removal")
    
    rbp_prefixes = config['features']['rbp_prefixes']
    rbp_features = [col for col in X.columns if any(col.startswith(prefix) for prefix in rbp_prefixes)]
    print(f"  Protecting {len(rbp_features)} RBP features")
    
    # Rule 3: PTC distance
    print("\nRule 3: PTC distance features")
    print("  → Keep PTC.2.EJC (preferred measure)")
    print("  → Drop PTC_dist_exon_end_0b (less informative)")
    
    if 'PTC_dist_exon_end_0b' in X.columns:
        drop_manual['PTC_dist_exon_end_0b'] = "Redundant - PTC.2.EJC kept (preferred)"
        print("  ✓ Will drop PTC_dist_exon_end_0b")
    
    # Apply drops
    features_to_drop = [col for col in drop_manual.keys() if col in X.columns]
    print(f"\nTotal manual drops: {len(features_to_drop)} features")
    
    X_cleaned = X.drop(columns=features_to_drop)
    print(f"After manual drops: {X_cleaned.shape}")
    
    return X_cleaned, rbp_prefixes


def automated_quality_checks(X_cleaned, rbp_prefixes, correlation_threshold):
    """Perform automated quality checks."""
    print("\n" + "=" * 80)
    print("AUTOMATED QUALITY CHECKS")
    print("=" * 80)
    
    automated_drops = []
    
    # Check 1: Duplicate columns
    print("\nCheck 1: Duplicate column names")
    dup_cols = X_cleaned.columns[X_cleaned.columns.duplicated()].tolist()
    if dup_cols:
        print(f"  ⚠️  Found {len(dup_cols)} duplicates:")
        for col in dup_cols:
            print(f"    - {col}")
        automated_drops.extend(dup_cols)
    else:
        print("  ✓ None found")
    
    # Check 2: Zero variance
    print("\nCheck 2: Zero/near-zero variance numeric features")
    zero_var = []
    for col in X_cleaned.select_dtypes(include=[np.number]).columns:
        std = X_cleaned[col].std()
        if std < 1e-10:
            zero_var.append(col)
    
    if len(zero_var) > 0:
        print(f"  Found {len(zero_var)} zero-variance features")
        for col in zero_var[:10]:
            print(f"    - {col}: constant = {X_cleaned[col].iloc[0]}")
        if len(zero_var) > 10:
            print(f"    ... and {len(zero_var)-10} more")
    else:
        print("  ✓ None found")
    
    automated_drops.extend(zero_var)
    
    # Check 4: High correlation (WITH RBP PROTECTION)
    print("\nCheck 4: High correlation (|r| > threshold)")
    print("  Note: RBP features are PROTECTED from removal")
    
    numeric_cols = X_cleaned.select_dtypes(include=[np.number]).columns
    high_corr_drops = []
    rbp_pairs_protected = 0
    
    if len(numeric_cols) > 1:
        print("  Computing correlation matrix...")
        corr_matrix = X_cleaned[numeric_cols].corr().abs()
        
        print("  Checking correlations...")
        for i in range(len(corr_matrix.columns)):
            for j in range(i+1, len(corr_matrix.columns)):
                if corr_matrix.iloc[i, j] > correlation_threshold:
                    feat1 = corr_matrix.columns[i]
                    feat2 = corr_matrix.columns[j]
                    
                    is_rbp1 = any(feat1.startswith(prefix) for prefix in rbp_prefixes)
                    is_rbp2 = any(feat2.startswith(prefix) for prefix in rbp_prefixes)
                    
                    if is_rbp1 and is_rbp2:
                        rbp_pairs_protected += 1
                    elif is_rbp1:
                        if feat2 not in high_corr_drops:
                            high_corr_drops.append(feat2)
                    elif is_rbp2:
                        if feat1 not in high_corr_drops:
                            high_corr_drops.append(feat1)
                    else:
                        if feat2 not in high_corr_drops:
                            high_corr_drops.append(feat2)
        
        print(f"\n  ✓ Protected {rbp_pairs_protected} RBP-RBP pairs")
        print(f"  ✓ Dropping {len(high_corr_drops)} correlated features")
    
    automated_drops.extend(high_corr_drops)
    automated_drops = list(set(automated_drops))
    print(f"\nTotal automated drops: {len(automated_drops)}")
    
    return automated_drops


def prepare_categorical_features(X_cleaned, config):
    """Prepare and handle categorical features."""
    print("\n" + "=" * 80)
    print("CATEGORICAL FEATURE PREPARATION")
    print("=" * 80)
    
    declared_cat = config['features']['categorical']
    print(f"Categorical features from config: {len(declared_cat)}")
    
    # Find binary features
    print("\nDetecting binary features...")
    binary_features = []
    for col in X_cleaned.select_dtypes(include=[np.number]).columns:
        unique_vals = X_cleaned[col].dropna().unique()
        if len(unique_vals) == 2:
            if set(unique_vals).issubset({0, 1, 0.0, 1.0, True, False}):
                binary_features.append(col)
    
    if len(binary_features) > 0:
        print(f"\nFound {len(binary_features)} binary features (will treat as categorical)")
    
    # Build complete categorical list
    cat_features = [c for c in declared_cat if c in X_cleaned.columns]
    cat_features.extend(binary_features)
    
    other_objs = [c for c in X_cleaned.columns 
                  if X_cleaned[c].dtype == "object" and c not in cat_features]
    cat_features.extend(other_objs)
    cat_features = list(set(cat_features))
    
    print(f"\nTotal categorical features: {len(cat_features)}")
    
    # Save missingness info BEFORE imputation
    print("\nRecording missingness before imputation...")
    feature_missingness = {}
    for col in X_cleaned.columns:
        non_null = X_cleaned[col].notna().sum()
        feature_missingness[col] = {
            'non_null_count': non_null,
            'non_null_pct': non_null / len(X_cleaned) * 100
        }
    
    # Convert categoricals to strings
    print("\nConverting categoricals to strings and handling missing values...")
    for col in cat_features:
        X_cleaned[col] = X_cleaned[col].astype(str)
        X_cleaned[col] = X_cleaned[col].replace('nan', 'MISSING')
        if X_cleaned[col].isna().any():
            X_cleaned[col] = X_cleaned[col].fillna('MISSING')
    
    # Impute numeric features
    numeric_cols = X_cleaned.select_dtypes(include=[np.number]).columns
    missing_numeric = X_cleaned[numeric_cols].isna().sum()
    missing_numeric = missing_numeric[missing_numeric > 0]
    
    if len(missing_numeric) > 0:
        print(f"\nImputing {len(missing_numeric)} numeric features with median...")
        for col in missing_numeric.index:
            median_val = X_cleaned[col].median()
            X_cleaned[col] = X_cleaned[col].fillna(median_val)
    
    print(f"\n✓ Remaining NaN: {X_cleaned.isna().sum().sum()}")
    print(f"✓ Final feature breakdown:")
    print(f"  Categorical: {len(cat_features)}")
    print(f"  Numeric: {len(X_cleaned.select_dtypes(include=[np.number]).columns)}")
    
    return X_cleaned, cat_features, feature_missingness


def create_final_dataset(X_cleaned, automated_drops):
    """Apply final drops and create final dataset."""
    print("\n" + "=" * 80)
    print("CREATING FINAL DATASET")
    print("=" * 80)
    
    X_final = X_cleaned.drop(columns=[col for col in automated_drops if col in X_cleaned.columns])
    
    print(f"\nFinal shape: {X_final.shape}")
    print(f"  Samples: {len(X_final)}")
    print(f"  Features: {X_final.shape[1]}")
    
    return X_final


def print_summary(X, X_cleaned, X_final, features_to_drop, automated_drops, 
                  rbp_prefixes, gc_features_to_drop):
    """Print cleaning summary."""
    print("\n" + "=" * 80)
    print("CLEANING SUMMARY")
    print("=" * 80)
    
    print(f"\nFeature counts:")
    print(f"  Original: {X.shape[1]}")
    print(f"  After manual drops: {X_cleaned.shape[1]}")
    print(f"  After automated drops: {X_final.shape[1]}")
    print(f"\n  Total removed: {X.shape[1] - X_final.shape[1]}")
    print(f"    - Manual: {len(features_to_drop)}")
    print(f"    - Automated: {len(automated_drops)}")
    
    # Verify custom rules
    print("\n" + "=" * 80)
    print("CUSTOM RULE VERIFICATION")
    print("=" * 80)
    
    rbp_retained = [col for col in X_final.columns 
                    if any(col.startswith(prefix) for prefix in rbp_prefixes)]
    print(f"\n✓ RBP features retained: {len(rbp_retained)}")
    
    au_retained = [col for col in X_final.columns if '_AU_content' in col or 'AUcontent' in col]
    print(f"\n✓ AU content features retained: {len(au_retained)}")
    
    gc_remaining = [col for col in X_final.columns if '_GC_content' in col or 'GCcontent' in col]
    if len(gc_remaining) == 0:
        print(f"✓ GC content features dropped: {len(gc_features_to_drop)}")
    else:
        print(f"⚠️  GC content features remaining: {len(gc_remaining)}")
    
    if 'PTC.2.EJC' in X_final.columns:
        print(f"\n✓ PTC.2.EJC retained (preferred measure)")
    
    if 'PTC_dist_exon_end_0b' not in X_final.columns:
        print(f"✓ PTC_dist_exon_end_0b dropped")
    
    # Check new features
    print("\n" + "=" * 80)
    print("NEW FEATURE VERIFICATION")
    print("=" * 80)
    
    codon_features = ['CodonOptimalityFraction_CDS', 'CodonOptimalityFraction_PTCpm100nt']
    for feat in codon_features:
        if feat in X_final.columns:
            non_null = X_final[feat].notna().sum()
            print(f"\n✓ {feat}")
            print(f"  Non-null: {non_null}/{len(X_final)} ({non_null/len(X_final)*100:.1f}%)")
            print(f"  Mean: {X_final[feat].mean():.4f}")
    
    if 'AverageCodonRNAUsage' not in X_final.columns:
        print(f"\n✓ AverageCodonRNAUsage dropped (replaced by codon optimality)")


def save_results(X_final, y, cat_features, feature_missingness, config):
    """Save cleaned dataset and feature list."""
    print("\n" + "=" * 80)
    print("SAVING RESULTS")
    print("=" * 80)
    
    target = config['model']['target']
    
    # Combine features and target
    final_df = X_final.copy()
    final_df[target] = y.values
    
    # Save cleaned dataset
    output_path = Path(config['data']['cleaned'])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(output_path, index=False)
    
    print(f"\n✓ Cleaned dataset saved: {output_path}")
    print(f"  Rows: {len(final_df)}")
    print(f"  Columns: {final_df.shape[1]}")
    
    # Create feature list
    feature_list = pd.DataFrame({
        'feature': X_final.columns,
        'dtype': [str(X_final[col].dtype) for col in X_final.columns],
        'is_categorical': [col in cat_features for col in X_final.columns],
        'non_null_count': [feature_missingness.get(col, {}).get('non_null_count', len(X_final)) 
                          for col in X_final.columns],
        'non_null_pct': [feature_missingness.get(col, {}).get('non_null_pct', 100.0) 
                        for col in X_final.columns]
    })
    
    # Add summary stats
    feature_list['mean'] = np.nan
    feature_list['std'] = np.nan
    feature_list['min'] = np.nan
    feature_list['max'] = np.nan
    
    for idx, row in feature_list.iterrows():
        if not row['is_categorical'] and row['dtype'] in ['float64', 'int64']:
            col = row['feature']
            feature_list.at[idx, 'mean'] = X_final[col].mean()
            feature_list.at[idx, 'std'] = X_final[col].std()
            feature_list.at[idx, 'min'] = X_final[col].min()
            feature_list.at[idx, 'max'] = X_final[col].max()
    
    feature_list_path = Path(config['data']['feature_list'])
    feature_list_path.parent.mkdir(parents=True, exist_ok=True)
    feature_list.to_csv(feature_list_path, index=False)
    
    print(f"\n✓ Feature list saved: {feature_list_path}")


def main():
    """Main execution function."""
    print("=" * 80)
    print("FEATURE CLEANING AND SELECTION")
    print("=" * 80)
    
    try:
        # Load config and data
        config = load_config()
        X, y = load_data(config)
        original_shape = X.shape
        
        # Apply manual drops
        X_cleaned, rbp_prefixes = apply_manual_drops(X, config)
        features_to_drop = original_shape[1] - X_cleaned.shape[1]
        
        # Save GC features count for summary
        gc_features_to_drop = [col for col in X.columns if '_GC_content' in col or 'GCcontent' in col]
        
        # Automated quality checks
        correlation_threshold = config['features']['correlation_threshold']
        automated_drops = automated_quality_checks(X_cleaned, rbp_prefixes, correlation_threshold)
        
        # Prepare categorical features
        X_cleaned, cat_features, feature_missingness = prepare_categorical_features(X_cleaned, config)
        
        # Create final dataset
        X_final = create_final_dataset(X_cleaned, automated_drops)
        
        # Print summary
        print_summary(X, X_cleaned, X_final, features_to_drop, automated_drops, 
                     rbp_prefixes, gc_features_to_drop)
        
        # Save results
        save_results(X_final, y, cat_features, feature_missingness, config)
        
        print("\n" + "=" * 80)
        print("COMPLETE! 🎉")
        print("=" * 80)
        print(f"\nDataset ready for model training:")
        print(f"  - Samples: {len(X_final)}")
        print(f"  - Features: {X_final.shape[1]}")
        print(f"\nReady for model training! (next script)")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
