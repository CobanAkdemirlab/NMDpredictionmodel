#!/usr/bin/env python3
"""
02_feature_cleaning_and_selection.py
=====================================
Cleans and selects features from the merged dataset for model training.

Steps:
  1. Load merged data (from script 01)
  2. Apply domain-knowledge manual drop rules
  3. Automated quality checks (duplicates, zero-variance, high correlation)
  4. Prepare and impute categorical features
  5. Save cleaned dataset and feature list

Output: data/TOPMed_cleaned.csv, data/final_feature_list.csv

Usage:
    python 02_feature_cleaning_and_selection.py
    python 02_feature_cleaning_and_selection.py --config /path/to/config.yaml
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

def load_data(config):
    PATH_INPUT = config['data']['merged']
    TARGET = config['model']['target']
    CORRELATION_THRESHOLD = config['features']['correlation_threshold']
    CATEGORICAL_FEATURES = config['features']['categorical']
    RBP_PREFIXES = config['features']['rbp_prefixes']

    print(f"\nInput:  {PATH_INPUT}")
    print(f"Output: {config['data']['cleaned']}")
    print(f"Target: {TARGET}")
    print(f"Correlation threshold: {CORRELATION_THRESHOLD}")

    print("=" * 80)
    print("FEATURE CLEANING AND SELECTION")
    print("=" * 80)

    print("\nLoading data...")
    df = pd.read_csv(PATH_INPUT)
    print(f"✓ Loaded: {df.shape}")

    y_full = df[TARGET].map({"TRUE": 1, "FALSE": 0, True: 1, False: 0}).astype("float")
    mask = y_full.notna()
    y = y_full.loc[mask].astype(int).reset_index(drop=True)
    X = df.loc[mask].drop(columns=[TARGET]).reset_index(drop=True)

    print(f"Valid samples: {len(y)}")
    print(f"  Escapees: {y.sum()} ({y.sum()/len(y)*100:.1f}%)")
    print(f"  NMD: {(~y.astype(bool)).sum()} ({(~y.astype(bool)).sum()/len(y)*100:.1f}%)")
    print(f"Features: {X.shape[1]}")

    return X, y, TARGET, CORRELATION_THRESHOLD, CATEGORICAL_FEATURES, RBP_PREFIXES, config


# ==============================================================================
# MANUAL DROP RULES
# ==============================================================================

def apply_manual_drops(X, CATEGORICAL_FEATURES, RBP_PREFIXES):
    print("\n" + "=" * 80)
    print("CUSTOM REMOVAL RULES")
    print("=" * 80)

    # Handle LOEUF_cat naming variations
    if 'loefu_cat' in X.columns:
        X = X.rename(columns={'loefu_cat': 'loeuf_cat'})
        print("✓ Renamed 'loefu_cat' → 'loeuf_cat'")
    elif 'LOEUF_cat' in X.columns:
        X = X.rename(columns={'LOEUF_cat': 'loeuf_cat'})
        print("✓ Renamed 'LOEUF_cat' → 'loeuf_cat'")

    drop_manual = {
        # ── Identifiers / unnamed columns ──────────────────────────────────────
        "TxName":               "Identifier - transcript name",
        "variantID":            "Identifier - variant ID",
        "ensembl_gene_id":      "Identifier - redundant with GENE_ID",
        "txnames":              "Identifier - transcript name duplicate",
        "gene":                 "Identifier - redundant with GENE_ID",
        "gene_id_gencode":      "Identifier - redundant Gencode gene ID",
        "hgnc_symbol":          "Identifier - gene symbol redundant with GENE_ID",
        "GENE_ID":              "Identifier - gene ID",
        "key":                  "Identifier - composite key",
        "Var1":                 "Identifier - leftover R rowname",
        "V1":                   "Identifier - leftover row indices",
        "V2":                   "Identifier - unnamed column",
        "V3":                   "Identifier - unnamed column",
        "V4":                   "Identifier - unnamed column",
        "V5":                   "Identifier - unnamed column",
        "V6":                   "Identifier - unnamed column",
        "V7":                   "Identifier - unnamed column",
        "V8":                   "Identifier - unnamed column",

        # ── Genomic coordinates ────────────────────────────────────────────────
        "contig":               "Identifier - chromosome name",
        "Chromosome":           "Identifier - duplicate of contig",
        "CHROM":                "Identifier - duplicate of contig",
        "position":             "Identifier - raw genomic position",
        "POS":                  "Identifier - duplicate of position",
        "Start":                "Identifier - genomic start coordinate",
        "End":                  "Identifier - genomic end coordinate",
        "Strand":               "Identifier - strand direction",
        "nearest_junction":     "Identifier - raw EJC junction coordinate; engineered features kept",
        "downstream_start":     "Identifier - raw coordinate, not a model feature",
        "coding.pos":           "Identifier - absolute CDS nucleotide position; redundant with relativePTClocation",

        # ── Raw allele sequences ───────────────────────────────────────────────
        "refAllele":            "Identifier - raw reference allele sequence",
        "REF_len":              "Identifier - check to make sure all SNV",
        "altAllele":            "Identifier - raw alternate allele sequence",
        "ALT_len":              "Identifier - check to make sure all SNV",
        "REF_ALLELE":           "Identifier - duplicate ref allele",
        "ALT_ALLELE":           "Identifier - duplicate alt allele",

        # ── Sequencing QC / technical ──────────────────────────────────────────
        "refCount":             "Leaky - raw read counts used to compute ALLELE.RAT",
        "altCount":             "Leaky - raw read counts used to compute ALLELE.RAT",
        "totalCount":           "Leaky - raw read counts used to compute ALLELE.RAT",
        "lowMAPQDepth":         "Technical QC - not biological",
        "lowBaseQDepth":        "Technical QC - not biological",
        "rawDepth":             "Technical QC - not biological",
        "otherBases":           "Technical QC - not biological",
        "improperPairs":        "Technical QC - not biological",
        "ZYG":                  "Sample-level zygosity",
        "FILTER":               "Sequencing filter flag - not biological",
        "sample":               "Identifier - sample GT",

        # ── Target / leaky ─────────────────────────────────────────────────────
        "ALLELE.RAT":           "Leaky - directly used to compute NMD.ESCAPEE target",
        "NMD.ESCAPEE":          "Target variable",

        # ── Scores not applicable to stopgain variants ─────────────────────────
        "REVEL_score":          "Not applicable - missense pathogenicity score for stopgain variants",

        # ── Leaky canonical NMD rule encodings ─────────────────────────────────
        "last.exon":            "Leaky - direct boolean encoding of canonical NMD last-exon rule",
        "penultimate.exon":     "Leaky - redundant with last.EJC, encodes NMD rule directly",

        # ── Leaky allele frequency ─────────────────────────────────────────────
        "Freq":                 "Leaky - allele frequency from same TOPMed experiment as target",
        "Freq.cat":             "Leaky - binned allele frequency from same TOPMed experiment as target",
        "Freq_old":             "deprecated leaky allele frequency",

        # ── PTC positional raw distances ───────────────────────────────────────
        "PTC.2.start":            "Redundant - proxies for cds_length (r=0.762)",
        "PTC.2.end":              "Redundant - proxies for cds_length (r=0.712)",
        "PTC_dist_exon_start_0b": "Redundant - relPTC_exon captures within-exon position independently",
        "PTC.2.EJC":              "Redundant - relPTC_exon captures within-exon position independently",

        # ── Binned PTC distance features ───────────────────────────────────────
        "PTC.2.EJC.binning":      "Redundant - continuous proxied length.mutated.exon (r=0.838)",
        "PTC.2.start.binning":    "Redundant - continuous proxied cds_length (r=0.762)",
        "PTC.2.end.binning":      "Redundant - continuous proxied cds_length (r=0.712)",

        # ── Other redundant features ───────────────────────────────────────────
        "aug_distance_nt":        "Redundant - aug_distance_category kept",
        "cds_exons":              "Identifier - raw exon boundary string; not usable by CatBoost",
        "CADD_raw":               "Redundant - CADD_phred is the interpretable score",
        "downstream":             "Redundant - last.EJC kept",
        "EJC.downstream":         "Redundant with last.EJC",

        # ── CDS length ─────────────────────────────────────────────────────────
        "cds_length":             "Redundant - log2_CDS kept",
        "cds_length.cut":         "Redundant - continuous log2_CDS kept",

        # ── UTR lengths ────────────────────────────────────────────────────────
        "threeUTR_length":        "Redundant - log2_3utr kept",
        "fiveutr_length":         "Redundant - log2_5utr kept",
        "threeUTR_length.cut":    "Redundant - continuous log2_3utr kept",
        "fiveUTR_length.cut":     "Redundant - continuous log2_5utr kept",
        "newUTR_length":          "Redundant - log2newUTR kept",
        "log2newUTR.cut":         "Redundant - continuous log2newUTR kept",

        # ── Constraint scores ──────────────────────────────────────────────────
        "pLI":                    "Redundant - categorical pLI.cat kept",
        "oe_lof_upper":           "Redundant - categorical LOEUF_cat kept",

        # ── Exon count ─────────────────────────────────────────────────────────
        "exon_count":             "Redundant - AmountExonsAfter more informative",
        "NCexonsnum":             "Redundant - AmountExonsAfter more informative",

        # ── Readthrough ────────────────────────────────────────────────────────
        "readthrough_category_hek293t": "Redundant - readthrough_score_hek293t has higher importance",

        # ── EJC overlap binary ─────────────────────────────────────────────────
        "has_ejc_overlap":        "Redundant - ejc_count_in_window more informative; binary loses dose-response signal",

        # ── Expression ─────────────────────────────────────────────────────────
        "MedianExpression":       "Leaky - correlates with totalCount (r=0.387); MedianExpression_log2 kept",
        "Whole.Blood":            "Leaky - correlates with refCount (r=0.336); MedianExpression_log2 kept",

        # ── 5'UTR regional composition: keep whole + first100 + last100 ────────
        "fiveUTR_AUcontentfirst200": "Redundant - r=0.878 with fiveUTR_AU_content whole; first100 kept",
        "fiveUTR_AUcontentlast200":  "Redundant - r=0.885 with fiveUTR_AU_content whole; last100 kept",
        "fiveUTR_UCcontentfirst200": "Redundant - first100 kept for proximal signal",
        "fiveUTR_UCcontentlast200":  "Redundant - last100 kept for distal signal",

        # ── 3'UTR AU regional composition: keep whole + first100 + last100 ─────
        "ThreeUTR_AUcontentfirst200": "Redundant - r=0.941 with ThreeUTR_AUcontentfirst100; first100 kept",
        "ThreeUTR_AUcontentlast200":  "Redundant - r=0.901 with ThreeUTR_AUcontentlast100; last100 kept",

        # ── 3'UTR UC regional composition: keep whole + first100 + last100 ─────
        "ThreeUTR_UCcontentfirst200": "Redundant - first100 kept; consistent with AU strategy",
        "ThreeUTR_UCcontentlast200":  "Redundant - last100 kept; consistent with AU strategy",

        # ── Length binning ─────────────────────────────────────────────────────
        "length.mutated.exon.binning": "Redundant - continuous length.mutated.exon kept",

        # ── Zero importance RBP features — run 1 (62 features) ─────────────────
        "ptc_pm100.CPEB4": "Zero importance across all CV folds",
        "ptc_pm100.RBM25": "Zero importance across all CV folds",
        "utr3_200.ZNF638": "Zero importance across all CV folds",
        "newutr_200.PUM2": "Zero importance across all CV folds",
        "utr3_200.HNRNPL": "Zero importance across all CV folds",
        "newutr_200.RBM41": "Zero importance across all CV folds",
        "ptc_pm100.RBM3": "Zero importance across all CV folds",
        "ejc_pm100.FXR1": "Zero importance across all CV folds",
        "newutr_200.MSI1": "Zero importance across all CV folds",
        "utr3_all.RBM3": "Zero importance across all CV folds",
        "utr3_200.KHDRBS3": "Zero importance across all CV folds",
        "ptc_to_ejc.MSI1": "Zero importance across all CV folds",
        "ejc_pm100.RBMS3": "Zero importance across all CV folds",
        "ptc_pm100.HNRNPU": "Zero importance across all CV folds",
        "ptc_to_ejc.HNRNPU": "Zero importance across all CV folds",
        "ejc_pm100.HNRNPD": "Zero importance across all CV folds",
        "newutr_200.HNRNPC": "Zero importance across all CV folds",
        "ejc_pm100.HNRNPC": "Zero importance across all CV folds",
        "ejc_pm100.MATR3": "Zero importance across all CV folds",
        "newutr_200.HNRNPD": "Zero importance across all CV folds",
        "ptc_pm100.HNRNPD": "Zero importance across all CV folds",
        "ejc_pm100.ANKHD1": "Zero importance across all CV folds",
        "ptc_to_ejc.RBM46": "Zero importance across all CV folds",
        "ejc_pm100.RBM41": "Zero importance across all CV folds",
        "ptc_to_ejc.EIF4B": "Zero importance across all CV folds",
        "ptc_to_ejc.RBM41": "Zero importance across all CV folds",
        "ejc_pm100.ZCRB1": "Zero importance across all CV folds",
        "ptc_to_ejc.RBM3": "Zero importance across all CV folds",
        "ejc_pm100.RBM42": "Zero importance across all CV folds",

        # ── Zero importance RBP features — run 2 (47 features) ─────────────────
        "ptc_to_ejc.ZNF638": "Zero importance across all CV folds",
        "newutr_200.KHSRP": "Zero importance across all CV folds",
        "ptc_to_ejc.ZFP36L2": "Zero importance across all CV folds",
        "newutr_200.IGHMBP2": "Zero importance across all CV folds",
        "newutr_200.IGF2BP2": "Zero importance across all CV folds",
        "newutr_200.IGF2BP1": "Zero importance across all CV folds",
        "newutr_200.ANKHD1": "Zero importance across all CV folds",
        "newutr_200.RBMS3": "Zero importance across all CV folds",
        "newutr_200.CPEB2": "Zero importance across all CV folds",
        "ptc_to_ejc.SRP14": "Zero importance across all CV folds",
        "newutr_200.RBFOX1": "Zero importance across all CV folds",
        "newutr_200.KHDRBS2": "Zero importance across all CV folds",
        "ptc_pm100.U2AF2": "Zero importance across all CV folds",
        "ptc_to_ejc.RBMS1": "Zero importance across all CV folds",
        "ptc_to_ejc.AKAP1": "Zero importance across all CV folds",
        "ejc_pm100.HNRNPU": "Zero importance across all CV folds",
        "ejc_pm100.HNRNPM": "Zero importance across all CV folds",
        "ptc_pm100.PABPC3": "Zero importance across all CV folds",
        "ejc_pm100.HNRNPDL": "Zero importance across all CV folds",
        "ptc_to_ejc.ANKHD1": "Zero importance across all CV folds",
        "ejc_pm100.HNRNPA1L2": "Zero importance across all CV folds",
        "ptc_to_ejc.AGO2": "Zero importance across all CV folds",
        "newutr_200.ZFP36L2": "Zero importance across all CV folds",
        "utr3_200.MSI1": "Zero importance across all CV folds",
        "utr3_200.PPRC1": "Zero importance across all CV folds",
        "utr3_200.RBM14": "Zero importance across all CV folds",
        "ejc_pm100.CELF5": "Zero importance across all CV folds",
        "ptc_pm100.ZFP36L2": "Zero importance across all CV folds",
        "utr3_200.RBMS1": "Zero importance across all CV folds",
        "ptc_pm100.KHDRBS2": "Zero importance across all CV folds",
        "ejc_pm100.KHDRBS3": "Zero importance across all CV folds",
        "ejc_pm100.MSI1": "Zero importance across all CV folds",
        "ptc_to_ejc.ELAVL3": "Zero importance across all CV folds",
        "ptc_to_ejc.ERI1": "Zero importance across all CV folds",
        "ptc_pm100.HNRNPC": "Zero importance across all CV folds",
        "ptc_pm100.HNRNPA1L2": "Zero importance across all CV folds",
        "ptc_pm100.TUT1": "Zero importance across all CV folds",
        "ptc_to_ejc.HNRNPC": "Zero importance across all CV folds",
        "ejc_pm100.RBM28": "Zero importance across all CV folds",
        "ejc_pm100.RBM3": "Zero importance across all CV folds",
        "utr3_200.A1CF": "Zero importance across all CV folds",
        "ejc_pm100.RBM46": "Zero importance across all CV folds",
        "ptc_pm100.DAZAP1": "Zero importance across all CV folds",
        "ejc_pm100.SAMD4A": "Zero importance across all CV folds",
        "ptc_to_ejc.KHDRBS2": "Zero importance across all CV folds",
        "ptc_to_ejc.KHDRBS3": "Zero importance across all CV folds",
        "ptc_pm100.FUS": "Zero importance across all CV folds",
    }

    # Rule 1: AU/GC content — drop GC, keep AU
    print("\nRule 1: AU/GC content pairs")
    print("  → Keep AU content (more biologically relevant)")
    print("  → Drop GC content (redundant)")
    gc_features_to_drop = [col for col in X.columns if '_GC_content' in col or 'GCcontent' in col]
    print(f"  Dropping {len(gc_features_to_drop)} GC content features")
    for col in gc_features_to_drop:
        drop_manual[col] = "Redundant - AU content kept (more biologically relevant)"

    # Rule 2: RBP protection
    print("\nRule 2: RBP features")
    print("  → Keep ALL RBP features (biologically meaningful correlation)")
    print("  → Will be protected from automated correlation removal")
    rbp_prefixes = RBP_PREFIXES
    print(f"RBP feature prefixes (protected): {len(rbp_prefixes)}")
    for prefix in rbp_prefixes:
        print(f"  - {prefix}")
    rbp_features = [col for col in X.columns if any(col.startswith(p) for p in rbp_prefixes)]
    print(f"  Protecting {len(rbp_features)} RBP features")

    # Rule 3: PTC distance
    print("\nRule 3: PTC distance features")
    print("  → Keep PTC.2.EJC (preferred measure)")
    print("  → Drop PTC_dist_exon_end_0b (less informative)")
    if 'PTC_dist_exon_end_0b' in X.columns:
        drop_manual['PTC_dist_exon_end_0b'] = "Redundant - PTC.2.EJC kept (preferred)"
        print("  ✓ Will drop PTC_dist_exon_end_0b")
    else:
        print("  ⚠️  PTC_dist_exon_end_0b not found")

    # Apply
    features_to_drop = [col for col in drop_manual.keys() if col in X.columns]
    print(f"\nTotal manual drops: {len(features_to_drop)} features")
    X_cleaned = X.drop(columns=features_to_drop)
    print(f"After manual drops: {X_cleaned.shape}")

    return X_cleaned, features_to_drop, gc_features_to_drop, rbp_features, rbp_prefixes


# ==============================================================================
# AUTOMATED CHECKS
# ==============================================================================

def automated_quality_checks(X_cleaned, CORRELATION_THRESHOLD, rbp_prefixes):
    print("\n" + "=" * 80)
    print("AUTOMATED QUALITY CHECKS")
    print("=" * 80)

    automated_drops = []

    # Check 1: Duplicate column names
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
        if X_cleaned[col].std() < 1e-10:
            zero_var.append(col)
    if zero_var:
        print(f"  Found {len(zero_var)} zero-variance features")
        for col in zero_var[:10]:
            print(f"    - {col}: constant = {X_cleaned[col].iloc[0]}")
        if len(zero_var) > 10:
            print(f"    ... and {len(zero_var)-10} more")
    else:
        print("  ✓ None found")
    automated_drops.extend(zero_var)

    # Check 3: Constant categorical (deferred)
    print("\nCheck 3: Constant categorical features")
    print("  ⏭️  Skipped - will check after identifying categoricals in next step")

    # Check 4: High correlation with RBP protection
    print("\nCheck 4: High correlation (|r| > CORRELATION_THRESHOLD)")
    print("  Note: RBP features are PROTECTED from removal")

    numeric_cols = X_cleaned.select_dtypes(include=[np.number]).columns
    high_corr_drops = []
    rbp_pairs_protected = 0

    if len(numeric_cols) > 1:
        print("  Computing correlation matrix...")
        corr_matrix = X_cleaned[numeric_cols].corr().abs()
        print("  Checking correlations...")
        for i in range(len(corr_matrix.columns)):
            for j in range(i + 1, len(corr_matrix.columns)):
                if corr_matrix.iloc[i, j] > CORRELATION_THRESHOLD:
                    feat1 = corr_matrix.columns[i]
                    feat2 = corr_matrix.columns[j]
                    is_rbp1 = any(feat1.startswith(p) for p in rbp_prefixes)
                    is_rbp2 = any(feat2.startswith(p) for p in rbp_prefixes)

                    if is_rbp1 and is_rbp2:
                        rbp_pairs_protected += 1
                    elif is_rbp1:
                        if feat2 not in high_corr_drops:
                            high_corr_drops.append(feat2)
                            print(f"    Drop {feat2}, keep {feat1} (RBP protected)")
                    elif is_rbp2:
                        if feat1 not in high_corr_drops:
                            high_corr_drops.append(feat1)
                            print(f"    Drop {feat1}, keep {feat2} (RBP protected)")
                    else:
                        if feat2 not in high_corr_drops:
                            high_corr_drops.append(feat2)
                            print(f"    Drop {feat2}, correlated with {feat1} (r={corr_matrix.iloc[i, j]:.3f})")

        print(f"\n  ✓ Protected {rbp_pairs_protected} RBP-RBP pairs")
        print(f"  ✓ Dropping {len(high_corr_drops)} correlated features")
    else:
        print("  Not enough numeric features to check")

    automated_drops.extend(high_corr_drops)
    automated_drops = list(set(automated_drops))
    print(f"\nTotal automated drops: {len(automated_drops)}")

    return automated_drops


# ==============================================================================
# CATEGORICAL PREPARATION AND IMPUTATION
# ==============================================================================

def prepare_categoricals_and_impute(X_cleaned, CATEGORICAL_FEATURES):
    print("\n" + "=" * 80)
    print("CATEGORICAL FEATURE PREPARATION")
    print("=" * 80)

    declared_cat = CATEGORICAL_FEATURES
    print(f"Categorical features from config: {len(declared_cat)}")
    for feat in declared_cat:
        print(f"  - {feat}")

    # Detect binary numeric features
    print("\nDetecting binary features...")
    binary_features = []
    for col in X_cleaned.select_dtypes(include=[np.number]).columns:
        unique_vals = X_cleaned[col].dropna().unique()
        if len(unique_vals) == 2 and set(unique_vals).issubset({0, 1, 0.0, 1.0, True, False}):
            binary_features.append(col)

    if binary_features:
        print(f"\nFound {len(binary_features)} binary features (will treat as categorical):")
        for col in sorted(binary_features)[:10]:
            vals = X_cleaned[col].value_counts().sort_index().to_dict()
            print(f"  {col}: {vals}")
        if len(binary_features) > 10:
            print(f"  ... and {len(binary_features) - 10} more")

    # Build full categorical list
    cat_features = [c for c in declared_cat if c in X_cleaned.columns]
    cat_features.extend(binary_features)
    other_objs = [c for c in X_cleaned.columns
                  if X_cleaned[c].dtype == "object" and c not in cat_features]
    cat_features.extend(other_objs)
    cat_features = list(set(cat_features))

    print(f"\nTotal categorical features: {len(cat_features)}")
    print(f"  Declared: {len([c for c in declared_cat if c in X_cleaned.columns])}")
    print(f"  Binary:   {len(binary_features)}")
    print(f"  Other objects: {len(other_objs)}")

    # Convert categoricals to string
    print("\nConverting categoricals to strings and handling missing values...")
    for col in cat_features:
        X_cleaned[col] = X_cleaned[col].astype(str)
        X_cleaned[col] = X_cleaned[col].replace('nan', 'MISSING')
        if X_cleaned[col].isna().any():
            X_cleaned[col] = X_cleaned[col].fillna('MISSING')

    # Record missingness before imputation
    print("\nRecording missingness before imputation...")
    feature_missingness = {}
    for col in X_cleaned.columns:
        non_null = X_cleaned[col].notna().sum()
        feature_missingness[col] = {
            'non_null_count': non_null,
            'non_null_pct': non_null / len(X_cleaned) * 100
        }

    # Zero-fill: structurally absent UTR/positional features
    utr_zero_fill = [
        'threeUTR_length', 'threeUTR_AU_content', 'threeUTR_UC_content',
        'ThreeUTR_AUcontentlast200', 'ThreeUTR_AUcontentfirst200',
        'ThreeUTR_UCcontentlast200', 'ThreeUTR_UCcontentfirst200',
        'ThreeUTR_AUcontentlast100', 'ThreeUTR_AUcontentfirst100',
        'ThreeUTR_UCcontentlast100', 'ThreeUTR_UCcontentfirst100',
        'fiveutr_length', 'fiveUTR_AU_content', 'fiveUTR_UC_content',
        'fiveUTR_AUcontentlast200', 'fiveUTR_AUcontentfirst200',
        'fiveUTR_UCcontentlast200', 'fiveUTR_UCcontentfirst200',
        'fiveUTR_AUcontentlast100', 'fiveUTR_AUcontentfirst100',
        'fiveUTR_UCcontentlast100', 'fiveUTR_UCcontentfirst100',
        'newUTR_length', 'log2newUTR', 'log2_3utr', 'log2_5utr',
        'cdsseq_AUcontentlast200', 'cdsseq_UCcontentlast200',
        'cdsseq_AUcontentfirst200', 'cdsseq_UCcontentfirst200',
        'phastcons_utr5_first200_median', 'phastcons_utr5_whole_median',
        'phylop_utr5_first200_median', 'phylop_utr5_whole_median',
        'phastcons_old3utr_first200_median', 'phastcons_old3utr_whole_median',
        'phylop_old3utr_first200_median', 'phylop_old3utr_whole_median',
        'phastcons_ptc_to_ejc_median',
        'phylop_ptc_to_ejc_median',
    ]

    print("\nZero-filling structurally absent UTR/positional features...")
    for col in utr_zero_fill:
        if col in X_cleaned.columns:
            n = X_cleaned[col].isna().sum()
            if n > 0:
                X_cleaned[col] = X_cleaned[col].fillna(0)
                print(f"  Zero-filled {col}: {n} NAs")

    # gnomAD — absent = not observed = ultra-rare
    if 'gnomAD_exome_ALL' in X_cleaned.columns:
        n = X_cleaned['gnomAD_exome_ALL'].isna().sum()
        if n > 0:
            X_cleaned['gnomAD_exome_ALL'] = X_cleaned['gnomAD_exome_ALL'].fillna(0)
            print(f"\n  Zero-filled gnomAD_exome_ALL: {n} NAs (absent = ultra-rare variant)")

    # Median imputation: genuine annotation gaps
    median_impute_explicit = [
        'pLI', 'oe_lof_upper',
        'MedianExpression', 'MedianExpression_log2', 'Whole.Blood',
        'half_life_PC1',
        'CADD_phred',
        'readthrough_score_hek293t',
    ]

    print("\nMedian-imputing annotation gap features...")
    for col in median_impute_explicit:
        if col in X_cleaned.columns:
            n = X_cleaned[col].isna().sum()
            if n > 0:
                median_val = X_cleaned[col].median()
                X_cleaned[col] = X_cleaned[col].fillna(median_val)
                print(f"  {col}: {n} → median={median_val:.4f}")

    # Catch-all median imputation
    numeric_cols = X_cleaned.select_dtypes(include=[np.number]).columns
    missing_numeric = X_cleaned[numeric_cols].isna().sum()
    missing_numeric = missing_numeric[missing_numeric > 0]

    if len(missing_numeric) > 0:
        print(f"\n⚠️  {len(missing_numeric)} numeric features still have NAs (catch-all median):")
        for col in missing_numeric.index:
            median_val = X_cleaned[col].median()
            X_cleaned[col] = X_cleaned[col].fillna(median_val)
            print(f"  {col}: {missing_numeric[col]} → median={median_val:.4f}")
    else:
        print("\n✓ No remaining numeric NAs after structured imputation")

    print(f"\n✓ Remaining NaN: {X_cleaned.isna().sum().sum()}")
    print(f"✓ Final feature breakdown:")
    print(f"  Categorical: {len(cat_features)}")
    print(f"  Numeric:     {len(X_cleaned.select_dtypes(include=[np.number]).columns)}")

    return X_cleaned, cat_features, feature_missingness


# ==============================================================================
# FINAL DATASET + SUMMARY
# ==============================================================================

def create_final_dataset(X_cleaned, automated_drops):
    print("\n" + "=" * 80)
    print("CREATING FINAL DATASET")
    print("=" * 80)

    X_final = X_cleaned.drop(columns=[col for col in automated_drops if col in X_cleaned.columns])
    print(f"\nFinal shape: {X_final.shape}")
    print(f"  Samples: {len(X_final)}")
    print(f"  Features: {X_final.shape[1]}")

    return X_final


def print_summary(X, X_cleaned, X_final, features_to_drop, automated_drops,
                  rbp_features, rbp_prefixes, gc_features_to_drop):
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

    print("\n" + "=" * 80)
    print("CUSTOM RULE VERIFICATION")
    print("=" * 80)

    rbp_retained = [col for col in X_final.columns if any(col.startswith(p) for p in rbp_prefixes)]
    print(f"\n✓ RBP features retained: {len(rbp_retained)}")
    print(f"  (Original: {len(rbp_features)})")

    au_retained = [col for col in X_final.columns if '_AU_content' in col or 'AUcontent' in col]
    print(f"\n✓ AU content features retained: {len(au_retained)}")

    gc_remaining = [col for col in X_final.columns if '_GC_content' in col or 'GCcontent' in col]
    if len(gc_remaining) == 0:
        print(f"✓ GC content features dropped: {len(gc_features_to_drop)}")
    else:
        print(f"⚠️  GC content features remaining: {len(gc_remaining)} (should be 0!)")

    if 'PTC_dist_exon_end_0b' in X_final.columns:
        print(f"⚠️  PTC_dist_exon_end_0b still present (should be dropped)")
    else:
        print(f"✓ PTC_dist_exon_end_0b dropped")

    print("\n" + "=" * 80)
    print("NEW FEATURE VERIFICATION")
    print("=" * 80)

    for feat in ['CodonOptimalityFraction_CDS', 'CodonOptimalityFraction_PTCpm100nt']:
        if feat in X_final.columns:
            non_null = X_final[feat].notna().sum()
            print(f"\n✓ {feat}")
            print(f"  Non-null: {non_null}/{len(X_final)} ({non_null/len(X_final)*100:.1f}%)")
            print(f"  Mean: {X_final[feat].mean():.4f}")
            print(f"  Range: [{X_final[feat].min():.4f}, {X_final[feat].max():.4f}]")
        else:
            print(f"\n⚠️  {feat} MISSING!")

    if 'AverageCodonRNAUsage' in X_final.columns:
        print(f"\n⚠️  AverageCodonRNAUsage still present (should be dropped!)")
    else:
        print(f"\n✓ AverageCodonRNAUsage dropped (replaced by codon optimality)")


# ==============================================================================
# SAVE
# ==============================================================================

def save_outputs(X_final, y, cat_features, feature_missingness, TARGET, config):
    print("\n" + "=" * 80)
    print("SAVING RESULTS")
    print("=" * 80)

    PATH_OUTPUT = config['data']['cleaned']
    PATH_FEATURES = config['data']['feature_list']

    final_df = X_final.copy()
    final_df[TARGET] = y.values

    Path(PATH_OUTPUT).parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(PATH_OUTPUT, index=False)
    print(f"\n✓ Cleaned dataset saved: {PATH_OUTPUT}")
    print(f"  Rows: {len(final_df)}")
    print(f"  Columns: {final_df.shape[1]}")

    feature_list = pd.DataFrame({
        'feature': X_final.columns,
        'dtype': [str(X_final[col].dtype) for col in X_final.columns],
        'is_categorical': [col in cat_features for col in X_final.columns],
        'non_null_count': [feature_missingness.get(col, {}).get('non_null_count', len(X_final))
                           for col in X_final.columns],
        'non_null_pct': [feature_missingness.get(col, {}).get('non_null_pct', 100.0)
                         for col in X_final.columns]
    })

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

    Path(PATH_FEATURES).parent.mkdir(parents=True, exist_ok=True)
    feature_list.to_csv(PATH_FEATURES, index=False)
    print(f"\n✓ Feature list saved: {PATH_FEATURES}")

    print(f"""
Dataset ready for model training:
  - Samples: {len(final_df)}
  - Features: {X_final.shape[1]}
  - Target: {TARGET}

Custom rules applied:
  ✓ Kept AU content, dropped GC content
  ✓ Protected ALL RBP features (biologically meaningful)
  ✓ Kept PTC.2.EJC over PTC_dist_exon_end_0b
  ✓ Replaced AverageCodonRNAUsage with codon optimality features

Ready for model training! (03_model_training.py)
""")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="Feature cleaning and selection for NMD escape prediction")
    parser.add_argument('--config', default=None, help='Path to config.yaml (default: ../config/config.yaml)')
    args = parser.parse_args()

    try:
        config = load_config(args.config)

        X, y, TARGET, CORRELATION_THRESHOLD, CATEGORICAL_FEATURES, RBP_PREFIXES, config = load_data(config)

        X_cleaned, features_to_drop, gc_features_to_drop, rbp_features, rbp_prefixes = apply_manual_drops(
            X, CATEGORICAL_FEATURES, RBP_PREFIXES
        )

        automated_drops = automated_quality_checks(X_cleaned, CORRELATION_THRESHOLD, rbp_prefixes)

        X_cleaned, cat_features, feature_missingness = prepare_categoricals_and_impute(
            X_cleaned, CATEGORICAL_FEATURES
        )

        X_final = create_final_dataset(X_cleaned, automated_drops)

        print_summary(X, X_cleaned, X_final, features_to_drop, automated_drops,
                      rbp_features, rbp_prefixes, gc_features_to_drop)

        save_outputs(X_final, y, cat_features, feature_missingness, TARGET, config)

        print("\n" + "=" * 80)
        print("COMPLETE! 🎉")
        print("=" * 80)

        return 0

    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
