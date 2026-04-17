#!/usr/bin/env python3
"""
01_data_loading_and_merging.py
==============================
Loads and merges all input data sources for NMD escape prediction:
  - Enhanced features (main variant CSV)
  - Codon optimality annotations
  - Readthrough scores
  - EJC occupancy
  - PTC AUG features
  - Conservation scores (PhastCons / PhyloP)

Filters indel variants, removes deprecated features, engineers AUG features,
and merges everything into a single merged CSV for downstream cleaning.

Output: data/TOPMed_merged.csv  (path set in config/config.yaml)

Usage:
    python 01_data_loading_and_merging.py
    python 01_data_loading_and_merging.py --config /path/to/config.yaml
"""

import pandas as pd
import numpy as np
from pathlib import Path
import yaml
import sys
import argparse


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

    # Resolve all data and output paths relative to repo root
    BASE_DIR = config_path.parent.parent
    for key, val in config['data'].items():
        config['data'][key] = str(BASE_DIR / val)
    for key, val in config['output'].items():
        if isinstance(val, str):
            config['output'][key] = str(BASE_DIR / val)

    print(f"✓ Configuration loaded from: {config_path}")
    return config


# ==============================================================================
# LOAD
# ==============================================================================

def load_datasets(config):
    print("=" * 80)
    print("LOADING DATASETS")
    print("=" * 80)

    df_enhanced = pd.read_csv(config['data']['enhanced_features'])
    print(f"✓ Enhanced features: {df_enhanced.shape}")

    df_codon = pd.read_csv(config['data']['codon_optimality'], sep='\t')
    print(f"✓ Codon optimality: {df_codon.shape}")

    df_readthrough = pd.read_csv(config['data']['readthrough'])
    print(f"✓ Readthrough: {df_readthrough.shape}")

    df_ejc = pd.read_csv(config['data']['ejc'], sep='\t')
    print(f"✓ EJC occupancy: {df_ejc.shape}")

    df_ptcaug = pd.read_csv(config['data']['ptc_aug'], sep='\t')
    print(f"✓ PTC AUG: {df_ptcaug.shape}")

    df_conservation = pd.read_csv(config['data']['conservation'])
    print(f"✓ Conservation: {df_conservation.shape}")

    return df_enhanced, df_codon, df_readthrough, df_ejc, df_ptcaug, df_conservation


# ==============================================================================
# INDEL FILTERING
# ==============================================================================

def identify_and_filter_indels(df_enhanced, df_codon, df_readthrough, df_ejc, df_ptcaug):
    print("\n" + "=" * 80)
    print("IDENTIFYING INDEL KEYS FROM READTHROUGH DATA")
    print("=" * 80)

    indel_keys = set()

    if 'V7' in df_readthrough.columns and 'V8' in df_readthrough.columns:
        print("\nScanning readthrough data for indels...")

        indel_mask = (
            df_readthrough['V7'].astype(str).str.contains('-', na=False) |
            df_readthrough['V8'].astype(str).str.contains('-', na=False)
        )

        indels_v7 = df_readthrough['V7'].astype(str).str.contains('-', na=False).sum()
        indels_v8 = df_readthrough['V8'].astype(str).str.contains('-', na=False).sum()
        total_indels = indel_mask.sum()

        print(f"  Indels detected:")
        print(f"    V7 column: {indels_v7}")
        print(f"    V8 column: {indels_v8}")
        print(f"    Total rows with indels: {total_indels}")

        indel_keys = set(df_readthrough.loc[indel_mask, 'variantID'].values)
        print(f"\n✓ Identified {len(indel_keys)} unique indel variantIDs")
        print(f"  Sample indel keys: {list(indel_keys)[:3]}")
    else:
        print("  ⚠️  V7/V8 columns not found - cannot identify indels")

    print("\n" + "=" * 80)
    print("FILTERING INDELS FROM ALL DATASETS")
    print("=" * 80)

    if len(indel_keys) > 0:
        # Enhanced
        before = len(df_enhanced)
        df_enhanced = df_enhanced[~df_enhanced['variantID'].isin(indel_keys)]
        print(f"\nEnhanced dataset: {before} → {len(df_enhanced)} (removed {before - len(df_enhanced)})")

        # Codon optimality — filter by allele length from PTC_ID
        if 'PTC_ID' in df_codon.columns:
            before = len(df_codon)
            def is_snv(ptc_id):
                parts = ptc_id.split('_')
                if len(parts) == 4:
                    return len(parts[2]) == 1 and len(parts[3]) == 1
                return False
            snv_mask = df_codon['PTC_ID'].apply(is_snv)
            indels_in_codon = (~snv_mask).sum()
            print(f"\nCodon optimality - filtering by allele length:")
            print(f"  Total: {before} | SNVs: {snv_mask.sum()} | Indels: {indels_in_codon}")
            if indels_in_codon > 0:
                print(f"  Sample indel PTC_IDs: {df_codon.loc[~snv_mask, 'PTC_ID'].head(3).tolist()}")
            df_codon = df_codon[snv_mask]
            print(f"  After filtering: {len(df_codon)} (removed {before - len(df_codon)})")

            # Standardize key column
            df_codon['variantID'] = df_codon['PTC_ID']
            print(f"  ✓ variantID set from PTC_ID for {len(df_codon)} SNVs")
            print(f"  Sample variantIDs: {df_codon['variantID'].head(3).tolist()}")

        # Readthrough
        before = len(df_readthrough)
        df_readthrough = df_readthrough[~df_readthrough['variantID'].isin(indel_keys)]
        print(f"\nReadthrough dataset: {before} → {len(df_readthrough)} (removed {before - len(df_readthrough)})")

        # EJC
        before = len(df_ejc)
        df_ejc = df_ejc.rename(columns={'PTC_ID': 'variantID'})
        df_ejc = df_ejc[~df_ejc['variantID'].isin(indel_keys)]
        print(f"\nEJC dataset: {before} → {len(df_ejc)} (removed {before - len(df_ejc)})")

        # PTC AUG
        before = len(df_ptcaug)
        df_ptcaug = df_ptcaug.rename(columns={'variant_id': 'variantID'})
        df_ptcaug = df_ptcaug[~df_ptcaug['variantID'].isin(indel_keys)]
        print(f"\nPTC AUG dataset: {before} → {len(df_ptcaug)} (removed {before - len(df_ptcaug)})")

        print("\n✓ Successfully filtered indel variants from all datasets")

    else:
        print("\n  No indels to filter")
        # Still standardize keys
        df_codon  = df_codon.rename(columns={'PTC_ID': 'variantID'})
        df_ejc    = df_ejc.rename(columns={'PTC_ID': 'variantID'})
        df_ptcaug = df_ptcaug.rename(columns={'variant_id': 'variantID'})

    return df_enhanced, df_codon, df_readthrough, df_ejc, df_ptcaug


# ==============================================================================
# REMOVE OLD FEATURES
# ==============================================================================

def remove_old_features(df_enhanced):
    print("\n" + "=" * 80)
    print("REMOVING OLD FEATURES")
    print("=" * 80)

    if 'AverageCodonRNAUsage' in df_enhanced.columns:
        print("\n✓ Dropping AverageCodonRNAUsage (replaced by codon optimality)")
        df_enhanced = df_enhanced.drop(columns=['AverageCodonRNAUsage'])
    else:
        print("\n⚠️  AverageCodonRNAUsage not found (already dropped or not present)")

    if 'median_half_life' in df_enhanced.columns:
        print("✓ Dropping median_half_life (replaced by half_life_PC1 in main CSV)")
        df_enhanced = df_enhanced.drop(columns=['median_half_life'])
    else:
        print("⚠️  median_half_life not found (already dropped or not present)")

    print(f"\nEnhanced dataset after drops: {df_enhanced.shape}")
    return df_enhanced


# ==============================================================================
# MERGES
# ==============================================================================

def merge_codon(df_result, df_codon):
    print("\n" + "=" * 80)
    print("MERGING CODON OPTIMALITY FEATURES")
    print("=" * 80)

    codon_cols = ['variantID', 'CodonOptimalityFraction_CDS', 'CodonOptimalityFraction_PTCpm100nt']
    print(f"\nMerging columns: {codon_cols}")

    common_keys = set(df_result['variantID']) & set(df_codon['variantID'])
    print(f"\nKey overlap: {len(common_keys)}/{len(df_result)} variants")

    df_result = df_result.merge(df_codon[codon_cols], on='variantID', how='left')
    print(f"\n✓ Merged codon optimality features")
    print(f"  Shape after merge: {df_result.shape}")

    for col in ['CodonOptimalityFraction_CDS', 'CodonOptimalityFraction_PTCpm100nt']:
        matched = df_result[col].notna().sum()
        pct = matched / len(df_result) * 100
        print(f"  {col}: {matched}/{len(df_result)} ({pct:.1f}%)")
        if matched > 0:
            print(f"    Mean: {df_result[col].mean():.4f}")
            print(f"    Range: {df_result[col].min():.4f} - {df_result[col].max():.4f}")

    return df_result


def merge_readthrough(df_result, df_readthrough, main_cols):
    print("\n" + "=" * 80)
    print("MERGING READTHROUGH FEATURES")
    print("=" * 80)

    rt_new_cols = ['variantID'] + [c for c in df_readthrough.columns
                   if c not in main_cols and c != 'variantID']

    readthrough_feature_cols = [c for c in rt_new_cols if 'readthrough' in c.lower()]
    print(f"\nNew readthrough feature columns: {readthrough_feature_cols}")
    print(f"Total new columns to merge: {len(rt_new_cols) - 1}")

    cols_to_merge = ['variantID', 'readthrough_score_hek293t', 'readthrough_category_hek293t']
    cols_to_merge = [c for c in cols_to_merge if c in df_readthrough.columns]
    print(f"Merging columns: {cols_to_merge}")

    df_result = df_result.merge(df_readthrough[cols_to_merge], on='variantID', how='left')
    print(f"\n✓ Merged readthrough features")
    print(f"  Shape after merge: {df_result.shape}")

    for col in cols_to_merge:
        if col != 'variantID':
            matched = df_result[col].notna().sum()
            pct = matched / len(df_result) * 100
            print(f"  {col}: {matched}/{len(df_result)} ({pct:.1f}%)")
            if 'category' in col:
                print(f"    Distribution:")
                print(df_result[col].value_counts().to_string(header=False))
            elif matched > 0:
                print(f"    Mean: {df_result[col].mean():.2f}")
                print(f"    Range: {df_result[col].min():.2f} - {df_result[col].max():.2f}")

    return df_result


def merge_ejc(df_result, df_ejc, main_cols):
    print("\n" + "=" * 80)
    print("MERGING EJC OCCUPANCY FEATURES")
    print("=" * 80)

    ejc_new_cols = ['variantID'] + [c for c in df_ejc.columns
                    if c not in main_cols and c != 'variantID']

    print(f"\nEJC columns to merge: {[c for c in ejc_new_cols if c != 'variantID']}")

    common_keys = set(df_result['variantID']) & set(df_ejc['variantID'])
    print(f"Key overlap: {len(common_keys)}/{len(df_result)} variants")

    df_result = df_result.merge(df_ejc[ejc_new_cols], on='variantID', how='left')
    print(f"\n✓ Merged EJC occupancy features")
    print(f"  Shape after merge: {df_result.shape}")

    for col in ['ejc_count_in_window', 'has_ejc_overlap']:
        if col in df_result.columns:
            matched = df_result[col].notna().sum()
            pct = matched / len(df_result) * 100
            print(f"  {col}: {matched}/{len(df_result)} ({pct:.1f}%)")

    return df_result


def engineer_aug_features(df_ptcaug):
    print("\n" + "=" * 80)
    print("STEP: Creating Enhanced AUG Features")
    print("=" * 80)

    aug_distance_na = df_ptcaug['nearest_inframe_aug_distance_nt'].isna().sum()
    aug_kozak_na = df_ptcaug['nearest_inframe_kozak_strength'].isna().sum()
    print(f"\nIn AUG dataframe:")
    print(f"  Distance NA: {aug_distance_na}/{len(df_ptcaug)} ({aug_distance_na/len(df_ptcaug)*100:.1f}%)")
    print(f"  Kozak NA:    {aug_kozak_na}/{len(df_ptcaug)} ({aug_kozak_na/len(df_ptcaug)*100:.1f}%)")
    print(f"  → These are variants with NO in-frame AUG (biologically meaningful!)")

    # Numeric distance
    df_ptcaug['aug_distance_nt'] = df_ptcaug['nearest_inframe_aug_distance_nt']

    # Distance categories
    def categorize_distance(dist):
        if pd.isna(dist):
            return 'no_inframe_AUG'
        elif dist < 50:
            return 'very_close'
        elif dist < 100:
            return 'close'
        elif dist < 200:
            return 'moderate'
        elif dist < 500:
            return 'far'
        else:
            return 'very_far'

    df_ptcaug['aug_distance_category'] = df_ptcaug['aug_distance_nt'].apply(categorize_distance)
    print(f"\nDistance categories created:")
    print(df_ptcaug['aug_distance_category'].value_counts().sort_index())

    # Kozak strength
    df_ptcaug['kozak_strength'] = df_ptcaug['nearest_inframe_kozak_strength'].fillna('no_inframe_AUG').astype(str)
    print(f"\nKozak strength categories:")
    print(df_ptcaug['kozak_strength'].value_counts().sort_index())

    # Frame AUGs
    df_ptcaug['has_plus1_aug'] = df_ptcaug['has_plus1_frame_aug'].fillna(False).astype(str)
    df_ptcaug['has_plus2_aug'] = df_ptcaug['has_plus2_frame_aug'].fillna(False).astype(str)

    # Combined frame status
    def get_aug_frame_status(row):
        plus1 = str(row['has_plus1_aug']) == 'True'
        plus2 = str(row['has_plus2_aug']) == 'True'
        if plus1 and plus2:
            return 'both_frames'
        elif plus1:
            return 'plus1_only'
        elif plus2:
            return 'plus2_only'
        else:
            return 'no_frame_AUG'

    df_ptcaug['aug_frame_status'] = df_ptcaug.apply(get_aug_frame_status, axis=1)
    print(f"\nAUG frame status:")
    print(df_ptcaug['aug_frame_status'].value_counts().sort_index())

    return df_ptcaug


def merge_aug(df_result, df_ptcaug):
    print("\n" + "=" * 80)
    print("MERGING PTC AUG FEATURES")
    print("=" * 80)

    aug_keep = [
        'variantID',
        'aug_distance_category',
        'aug_distance_nt',
        'kozak_strength',
        'has_plus1_aug',
        'has_plus2_aug',
        'aug_frame_status',
    ]
    aug_keep = [c for c in aug_keep if c in df_ptcaug.columns]

    print(f"\nAUG columns to merge: {[c for c in aug_keep if c != 'variantID']}")

    common_keys = set(df_result['variantID']) & set(df_ptcaug['variantID'])
    print(f"Key overlap: {len(common_keys)}/{len(df_result)} variants")

    df_result = df_result.merge(df_ptcaug[aug_keep], on='variantID', how='left')
    print(f"\n✓ Merged PTC AUG features")
    print(f"  Shape after merge: {df_result.shape}")

    for col in ['aug_distance_category', 'kozak_strength', 'aug_frame_status']:
        if col in df_result.columns:
            matched = df_result[col].notna().sum()
            print(f"  {col}: {matched}/{len(df_result)} non-null")

    return df_result


def merge_conservation(df_result, df_conservation):
    print("\n" + "=" * 80)
    print("MERGING CONSERVATION SCORE FEATURES")
    print("=" * 80)

    cons_new_cols = ['variantID'] + [c for c in df_conservation.columns
                     if c.startswith('phastcons_') or c.startswith('phylop_')]

    print(f"\nConservation columns to merge: {len(cons_new_cols) - 1}")
    print(f"  PhastCons: {[c for c in cons_new_cols if c.startswith('phastcons_')]}")
    print(f"  PhyloP:    {[c for c in cons_new_cols if c.startswith('phylop_')]}")

    common_keys = set(df_result['variantID']) & set(df_conservation['variantID'])
    print(f"\nKey overlap: {len(common_keys)}/{len(df_result)} variants")

    df_result = df_result.merge(df_conservation[cons_new_cols], on='variantID', how='left')
    print(f"\n✓ Merged conservation features")
    print(f"  Shape after merge: {df_result.shape}")

    sample_col = 'phastcons_ptc_100bp_median'
    if sample_col in df_result.columns:
        matched = df_result[sample_col].notna().sum()
        print(f"  {sample_col}: {matched}/{len(df_result)} non-null ({matched/len(df_result)*100:.1f}%)")

    return df_result


# ==============================================================================
# VERIFY AND SAVE
# ==============================================================================

def verify_results(df_result, original_shape):
    print("\n" + "=" * 80)
    print("FINAL VERIFICATION")
    print("=" * 80)

    print(f"\nFinal dataset shape: {df_result.shape}")
    print(f"  Original enhanced: {original_shape[1]} columns")
    print(f"  Final: {df_result.shape[1]} columns")
    print(f"  Net change: {df_result.shape[1] - original_shape[1]:+d} columns")

    print("\n✓ New features added:")
    new_features_expected = [
        'CodonOptimalityFraction_CDS',
        'CodonOptimalityFraction_PTCpm100nt',
        'readthrough_score_hek293t',
        'readthrough_category_hek293t',
        'ejc_count_in_window',
        'has_ejc_overlap',
        'aug_distance_category',
        'kozak_strength',
        'aug_frame_status',
        'has_plus1_aug',
        'has_plus2_aug',
        'phastcons_ptc_100bp_median',
        'phylop_ptc_100bp_median',
        'half_life_PC1'
    ]
    for col in new_features_expected:
        if col in df_result.columns:
            non_null = df_result[col].notna().sum()
            print(f"  ✔ {col}: {non_null}/{len(df_result)} non-null")
        else:
            print(f"  ✘ {col}: MISSING")

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
    print("\n" + "=" * 80)
    print("SAVING MERGED DATASET")
    print("=" * 80)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df_result.to_csv(output_path, index=False)
    print(f"\n✓ Saved: {output_path}")
    print(f"  Rows: {len(df_result)}")
    print(f"  Columns: {df_result.shape[1]}")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="Data loading and merging for NMD escape prediction")
    parser.add_argument('--config', default=None, help='Path to config.yaml (default: ../config/config.yaml)')
    args = parser.parse_args()

    print("=" * 80)
    print("DATA LOADING AND MERGING")
    print("=" * 80)

    try:
        config = load_config(args.config)

        print(f"\nInput files:")
        print(f"  Enhanced:      {config['data']['enhanced_features']}")
        print(f"  Codon opt:     {config['data']['codon_optimality']}")
        print(f"  Readthrough:   {config['data']['readthrough']}")
        print(f"  EJC:           {config['data']['ejc']}")
        print(f"  PTC AUG:       {config['data']['ptc_aug']}")
        print(f"  Conservation:  {config['data']['conservation']}")
        print(f"\nOutput:")
        print(f"  {config['data']['merged']}")

        # Load
        df_enhanced, df_codon, df_readthrough, df_ejc, df_ptcaug, df_conservation = load_datasets(config)
        original_shape = df_enhanced.shape

        # Track main columns before merging (used to avoid duplicating cols from annotations)
        main_cols = set(df_enhanced.columns)

        # Filter indels and standardize keys
        df_enhanced, df_codon, df_readthrough, df_ejc, df_ptcaug = identify_and_filter_indels(
            df_enhanced, df_codon, df_readthrough, df_ejc, df_ptcaug
        )

        # Remove deprecated features
        df_enhanced = remove_old_features(df_enhanced)

        # Merge codon optimality (initializes df_result)
        df_result = merge_codon(df_enhanced, df_codon)

        # Merge readthrough
        df_result = merge_readthrough(df_result, df_readthrough, main_cols)

        # Merge EJC
        df_result = merge_ejc(df_result, df_ejc, main_cols)

        # Engineer and merge AUG features
        df_ptcaug = engineer_aug_features(df_ptcaug)
        df_result = merge_aug(df_result, df_ptcaug)

        # Merge conservation
        df_result = merge_conservation(df_result, df_conservation)

        # Verify
        verify_results(df_result, original_shape)

        # Save
        save_output(df_result, config['data']['merged'])

        print("\n" + "=" * 80)
        print("COMPLETE! 🎉")
        print("=" * 80)
        print("\nReady for feature cleaning and selection (02_feature_cleaning_and_selection.py)")

        return 0

    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
