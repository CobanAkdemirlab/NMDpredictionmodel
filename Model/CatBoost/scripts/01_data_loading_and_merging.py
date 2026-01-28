#!/usr/bin/env python3
"""
Data Loading and Merging Script
================================
This script loads and merges multiple data sources for NMD escape prediction:
- Enhanced features
- Codon optimality data
- Half-life (PC1) features
- Readthrough scores

It filters out indel variants and replaces outdated features with updated versions.

Requirements:
    - Python 3.12+
    - pandas
    - numpy
    - pyyaml

Usage:
    python 01_data_loading_and_merging.py
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


def load_datasets(config):
    """Load all required datasets."""
    print("=" * 80)
    print("LOADING DATASETS")
    print("=" * 80)
    
    # Load enhanced features
    df_enhanced = pd.read_csv(config['data']['enhanced_features'])
    print(f"✓ Enhanced features: {df_enhanced.shape}")
    
    # Load codon optimality data
    df_codon = pd.read_csv(config['data']['codon_optimality'], sep='\t')
    print(f"✓ Codon optimality: {df_codon.shape}")
    
    # Load half-life (PC1) data
    df_halflife = pd.read_csv(config['data']['halflife'])
    print(f"✓ Half-life (PC1): {df_halflife.shape}")
    
    # Load readthrough data
    df_readthrough = pd.read_csv(config['data']['readthrough'])
    print(f"✓ Readthrough: {df_readthrough.shape}")
    
    return df_enhanced, df_codon, df_halflife, df_readthrough


def identify_indels(df_readthrough):
    """Identify indel variants from readthrough data."""
    print("\n" + "=" * 80)
    print("IDENTIFYING INDEL KEYS FROM READTHROUGH DATA")
    print("=" * 80)
    
    indel_keys = set()
    
    if 'V7' in df_readthrough.columns and 'V8' in df_readthrough.columns:
        print(f"\nScanning readthrough data for indels...")
        
        # Identify rows with indels
        indel_mask = ((df_readthrough['V7'].astype(str).str.contains('-', na=False)) | 
                      (df_readthrough['V8'].astype(str).str.contains('-', na=False)))
        
        indels_v7 = df_readthrough['V7'].astype(str).str.contains('-', na=False).sum()
        indels_v8 = df_readthrough['V8'].astype(str).str.contains('-', na=False).sum()
        total_indels = indel_mask.sum()
        
        print(f"  Indels detected:")
        print(f"    V7 column: {indels_v7}")
        print(f"    V8 column: {indels_v8}")
        print(f"    Total rows with indels: {total_indels}")
        
        # Extract keys of indel variants
        indel_keys = set(df_readthrough.loc[indel_mask, 'key'].values)
        
        print(f"\n✓ Identified {len(indel_keys)} unique indel variant keys")
        if len(indel_keys) > 0:
            print(f"  Sample indel keys: {list(indel_keys)[:3]}")
    else:
        print("  ⚠️  V7/V8 columns not found - cannot identify indels")
    
    return indel_keys


def filter_indels(df_enhanced, df_codon, df_halflife, df_readthrough, indel_keys):
    """Filter indel variants from all datasets."""
    print("\n" + "=" * 80)
    print("FILTERING INDELS FROM ALL DATASETS")
    print("=" * 80)
    
    if len(indel_keys) == 0:
        print("\n  No indels to filter")
        return df_enhanced, df_codon, df_halflife, df_readthrough
    
    # Filter enhanced dataset
    before = len(df_enhanced)
    df_enhanced = df_enhanced[~df_enhanced['key'].isin(indel_keys)]
    removed = before - len(df_enhanced)
    print(f"\nEnhanced dataset:")
    print(f"  Before: {before} variants")
    print(f"  After: {len(df_enhanced)} variants")
    print(f"  Removed: {removed} indels")
    
    # Filter codon optimality dataset based on PTC_ID format
    if 'PTC_ID' in df_codon.columns:
        before = len(df_codon)
        
        def is_snv(ptc_id):
            parts = ptc_id.split('_')
            if len(parts) == 4:
                ref = parts[2]
                alt = parts[3]
                return len(ref) == 1 and len(alt) == 1
            return False
        
        snv_mask = df_codon['PTC_ID'].apply(is_snv)
        indels_in_codon = (~snv_mask).sum()
        
        print(f"\nCodon optimality - filtering by allele length:")
        print(f"  Total variants: {before}")
        print(f"  SNVs (len=1 for both ref/alt): {snv_mask.sum()}")
        print(f"  Indels (len>1 for ref or alt): {indels_in_codon}")
        
        if indels_in_codon > 0:
            print(f"  Sample indel PTC_IDs: {df_codon.loc[~snv_mask, 'PTC_ID'].head(3).tolist()}")
        
        df_codon = df_codon[snv_mask]
        removed = before - len(df_codon)
        
        print(f"  After filtering: {len(df_codon)} variants")
        print(f"  Removed: {removed} indels")
        
        # Create key column from cleaned SNV data
        def convert_ptc_id_to_key(ptc_id):
            parts = ptc_id.split('_')
            chrom = parts[0]
            pos = parts[1]
            ref = parts[2]
            alt = parts[3]
            return f"{chrom}:{pos}_{ref}>{alt}"
        
        df_codon['key'] = df_codon['PTC_ID'].apply(convert_ptc_id_to_key)
        print(f"  ✓ Created 'key' from {len(df_codon)} SNVs")
        print(f"  Sample keys: {df_codon['key'].head(3).tolist()}")
    
    # Filter half-life dataset
    before = len(df_halflife)
    df_halflife = df_halflife[~df_halflife['key'].isin(indel_keys)]
    removed = before - len(df_halflife)
    print(f"\nHalf-life dataset:")
    print(f"  Before: {before} variants")
    print(f"  After: {len(df_halflife)} variants")
    print(f"  Removed: {removed} indels")
    
    # Filter readthrough dataset
    before = len(df_readthrough)
    df_readthrough = df_readthrough[~df_readthrough['key'].isin(indel_keys)]
    removed = before - len(df_readthrough)
    print(f"\nReadthrough dataset:")
    print(f"  Before: {before} variants")
    print(f"  After: {len(df_readthrough)} variants")
    print(f"  Removed: {removed} indels")
    
    print(f"\n✓ Successfully filtered indel variants from all datasets")
    
    return df_enhanced, df_codon, df_halflife, df_readthrough


def remove_old_features(df_enhanced):
    """Remove outdated features that will be replaced."""
    print("\n" + "=" * 80)
    print("REMOVING OLD FEATURES")
    print("=" * 80)
    
    # Drop AverageCodonRNAUsage (wrong data)
    if 'AverageCodonRNAUsage' in df_enhanced.columns:
        print(f"\n✓ Dropping AverageCodonRNAUsage (incorrect data)")
        print(f"  Old values (first 5): {df_enhanced['AverageCodonRNAUsage'].head().tolist()}")
        df_enhanced = df_enhanced.drop(columns=['AverageCodonRNAUsage'])
    else:
        print(f"\n⚠️  AverageCodonRNAUsage not found (already dropped?)")
    
    # Drop median_half_life if present
    if 'median_half_life' in df_enhanced.columns:
        print(f"✓ Dropping median_half_life (replaced by PC1)")
        df_enhanced = df_enhanced.drop(columns=['median_half_life'])
    else:
        print(f"⚠️  median_half_life not found (already dropped?)")
    
    print(f"\nEnhanced dataset after drops: {df_enhanced.shape}")
    
    return df_enhanced


def merge_datasets(df_enhanced, df_codon, df_halflife, df_readthrough):
    """Merge all datasets."""
    print("\n" + "=" * 80)
    print("MERGING DATASETS")
    print("=" * 80)
    
    # Merge codon optimality
    print("\n1. Merging codon optimality features...")
    codon_cols = ['key', 'CodonOptimalityFraction_CDS', 'CodonOptimalityFraction_PTCpm100nt']
    common_keys = set(df_enhanced['key']) & set(df_codon['key'])
    print(f"  Key overlap: {len(common_keys)}/{len(df_enhanced)} variants")
    
    df_result = df_enhanced.merge(df_codon[codon_cols], on='key', how='left')
    print(f"  Shape after merge: {df_result.shape}")
    
    for col in ['CodonOptimalityFraction_CDS', 'CodonOptimalityFraction_PTCpm100nt']:
        matched = df_result[col].notna().sum()
        pct = matched / len(df_result) * 100
        print(f"  {col}: {matched}/{len(df_result)} ({pct:.1f}%)")
    
    # Merge half-life
    print("\n2. Merging half-life (PC1) features...")
    if 'half_life_PC1' in df_halflife.columns:
        df_result = df_result.merge(df_halflife[['key', 'half_life_PC1']], on='key', how='left')
        matched = df_result['half_life_PC1'].notna().sum()
        print(f"  ✓ Merged half_life_PC1: {matched}/{len(df_result)} matched")
    else:
        print(f"  ⚠️  half_life_PC1 not found in dataset")
    
    # Merge readthrough
    print("\n3. Merging readthrough features...")
    cols_to_merge = ['key', 'readthrough_score_hek293t', 'readthrough_category_hek293t']
    cols_to_merge = [c for c in cols_to_merge if c in df_readthrough.columns]
    
    df_result = df_result.merge(df_readthrough[cols_to_merge], on='key', how='left')
    print(f"  Shape after merge: {df_result.shape}")
    
    for col in cols_to_merge:
        if col != 'key':
            matched = df_result[col].notna().sum()
            pct = matched / len(df_result) * 100
            print(f"  {col}: {matched}/{len(df_result)} ({pct:.1f}%)")
    
    return df_result


def verify_results(df_result, original_shape):
    """Verify final results."""
    print("\n" + "=" * 80)
    print("FINAL VERIFICATION")
    print("=" * 80)
    
    print(f"\nFinal dataset shape: {df_result.shape}")
    print(f"  Original enhanced: {original_shape[1]} columns")
    print(f"  Final: {df_result.shape[1]} columns")
    print(f"  Net change: {df_result.shape[1] - original_shape[1]:+d} columns")
    
    # Verify old features are gone
    print("\n✓ Removed features:")
    if 'AverageCodonRNAUsage' not in df_result.columns:
        print("  - AverageCodonRNAUsage (dropped)")
    if 'median_half_life' not in df_result.columns:
        print("  - median_half_life (dropped)")
    
    # Verify new features are present
    print("\n✓ New features added:")
    for col in ['CodonOptimalityFraction_CDS', 'CodonOptimalityFraction_PTCpm100nt',
                'half_life_PC1', 'readthrough_score_hek293t', 'readthrough_category_hek293t']:
        if col in df_result.columns:
            non_null = df_result[col].notna().sum()
            print(f"  - {col}: {non_null}/{len(df_result)} non-null")
    
    # Check for missing data
    print(f"\nMissing data summary:")
    missing_summary = df_result.isnull().sum()
    missing_summary = missing_summary[missing_summary > 0].sort_values(ascending=False)
    if len(missing_summary) > 0:
        print(f"  Columns with missing data: {len(missing_summary)}")
        print(f"  Top 10 columns with most missing:")
        for col, count in missing_summary.head(10).items():
            pct = count / len(df_result) * 100
            print(f"    {col}: {count} ({pct:.1f}%)")
    else:
        print("  No missing data!")


def save_output(df_result, output_path):
    """Save merged dataset."""
    print("\n" + "=" * 80)
    print("SAVING MERGED DATASET")
    print("=" * 80)
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    df_result.to_csv(output_path, index=False)
    print(f"\n✓ Saved: {output_path}")
    print(f"  Rows: {len(df_result)}")
    print(f"  Columns: {df_result.shape[1]}")


def main():
    """Main execution function."""
    print("=" * 80)
    print("DATA LOADING AND MERGING")
    print("=" * 80)
    
    try:
        # Load configuration
        print("\n✓ Configuration loaded")
        config = load_config()
        
        # Load datasets
        df_enhanced, df_codon, df_halflife, df_readthrough = load_datasets(config)
        original_shape = df_enhanced.shape
        
        # Identify and filter indels
        indel_keys = identify_indels(df_readthrough)
        df_enhanced, df_codon, df_halflife, df_readthrough = filter_indels(
            df_enhanced, df_codon, df_halflife, df_readthrough, indel_keys
        )
        
        # Remove old features
        df_enhanced = remove_old_features(df_enhanced)
        
        # Merge datasets
        df_result = merge_datasets(df_enhanced, df_codon, df_halflife, df_readthrough)
        
        # Verify results
        verify_results(df_result, original_shape)
        
        # Save output
        save_output(df_result, config['data']['merged'])
        
        print("\n" + "=" * 80)
        print("COMPLETE! 🎉")
        print("=" * 80)
        print("\nReady for feature cleaning and selection (next script)")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
