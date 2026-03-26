#!/usr/bin/env python3
"""
Setup Verification for NMD Escape Prediction Scripts
=====================================================
Checks that config, input data, and dependencies are all in place
before running the three pipeline scripts.

Usage:
    python setup.py
    python setup.py --config /path/to/config.yaml
"""

import sys
import argparse
from pathlib import Path


# ==============================================================================
# CONFIG
# ==============================================================================

def find_config(config_path=None):
    """Locate the config file relative to this script."""
    if config_path:
        p = Path(config_path)
    else:
        p = Path(__file__).resolve().parent.parent / "config" / "config.yaml"
    return p


def check_config(config_path=None):
    """Check config file exists and load it."""
    print("Checking configuration...")

    config_file = find_config(config_path)

    if not config_file.exists():
        print(f"  ✗ Config file not found: {config_file}")
        print(f"  ⚠️  Edit config/config.yaml with your data paths before running.")
        return False, None

    try:
        import yaml
        with open(config_file) as f:
            config = yaml.safe_load(f)
        print(f"  ✓ Found: {config_file}")
        return True, config
    except Exception as e:
        print(f"  ✗ Failed to parse config: {e}")
        return False, None


# ==============================================================================
# INPUT DATA
# ==============================================================================

def check_input_data(config):
    """Check that all required input data files exist."""
    print("\nChecking input data files...")

    BASE_DIR = find_config().parent.parent

    input_keys = {
        'enhanced_features': 'Main variant CSV',
        'codon_optimality':  'Codon optimality TSV',
        'readthrough':       'Readthrough scores CSV',
        'ejc':               'EJC occupancy TSV',
        'ptc_aug':           'PTC AUG features TSV',
        'conservation':      'Conservation scores CSV',
    }

    all_ok = True
    for key, label in input_keys.items():
        raw_path = config.get('data', {}).get(key)
        if raw_path is None:
            print(f"  ✗ {label}: key '{key}' missing from config")
            all_ok = False
            continue

        path = BASE_DIR / raw_path
        if path.exists():
            print(f"  ✓ {label}: {path.name}")
        else:
            print(f"  ✗ {label}: not found")
            print(f"      Expected: {path}")
            all_ok = False

    if not all_ok:
        print("\n  ⚠️  Missing input files — see data/README.md for sourcing instructions.")

    return all_ok


# ==============================================================================
# OUTPUT DIRECTORIES
# ==============================================================================

def check_output_directories(config):
    """Create output directories if they do not exist."""
    print("\nChecking output directories...")

    BASE_DIR = find_config().parent.parent

    output_keys = ['models_dir', 'cv_models_dir', 'results_dir', 'figures_dir']
    for key in output_keys:
        raw_path = config.get('output', {}).get(key)
        if raw_path is None:
            print(f"  ⚠️  output.{key} not found in config")
            continue
        path = BASE_DIR / raw_path
        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)
            print(f"  ✓ Created: {path.relative_to(BASE_DIR)}")
        else:
            print(f"  ✓ Exists:  {path.relative_to(BASE_DIR)}")


# ==============================================================================
# DEPENDENCIES
# ==============================================================================

def check_dependencies():
    """Check all required packages are importable."""
    print("\nChecking dependencies...")

    required = [
        ("pandas",      "pandas",       ">=2.1.0"),
        ("numpy",       "numpy",        ">=1.24.0"),
        ("yaml",        "pyyaml",       ">=6.0"),
        ("catboost",    "catboost",     ">=1.2.8"),
        ("sklearn",     "scikit-learn", ">=1.3.0"),
        ("matplotlib",  "matplotlib",   ">=3.7.0"),
        ("seaborn",     "seaborn",      ">=0.12.0"),
        ("shap",        "shap",         ">=0.43.0"),
        ("tqdm",        "tqdm",         ">=4.66.0"),
        ("joblib",      "joblib",       ">=1.3.0"),
        ("optuna",      "optuna",       ">=3.0.0"),
    ]

    missing = []
    for import_name, package_name, version_req in required:
        try:
            mod = __import__(import_name)
            version = getattr(mod, '__version__', 'unknown')
            print(f"  ✓ {package_name:<20} (installed: {version})")
        except ImportError:
            print(f"  ✗ {package_name:<20} MISSING  (required {version_req})")
            missing.append(package_name)

    if missing:
        print(f"\n  ⚠️  Missing: {', '.join(missing)}")
        print(f"     Install with: pip install -r requirements.txt")
        return False

    print("\n  ✓ All dependencies installed!")
    return True


# ==============================================================================
# SUMMARY
# ==============================================================================

def print_summary(config_ok, data_ok, deps_ok):
    print("\n" + "=" * 80)

    all_ok = config_ok and data_ok and deps_ok

    if all_ok:
        print("SETUP VERIFIED — READY TO RUN!")
        print("=" * 80)
        print("""
Run the pipeline scripts in order:

    python 01_data_loading_and_merging.py
    python 02_feature_cleaning_and_selection.py
    python 03_model_training.py

Each script accepts an optional --config flag:

    python 01_data_loading_and_merging.py --config /path/to/config.yaml
""")
    else:
        print("SETUP INCOMPLETE — fix the issues above before running.")
        print("=" * 80)
        if not config_ok:
            print("\n  1. Ensure config/config.yaml exists and is valid YAML.")
        if not data_ok:
            print("\n  2. Place input data files in the locations specified in config.yaml.")
            print("     See data/README.md for sourcing instructions.")
        if not deps_ok:
            print("\n  3. Install missing dependencies:")
            print("     pip install -r requirements.txt")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="NMD escape prediction — setup verification")
    parser.add_argument('--config', default=None, help='Path to config.yaml')
    args = parser.parse_args()

    print("=" * 80)
    print("NMD Escape Prediction Model — Setup Verification")
    print("=" * 80)
    print()

    config_ok, config = check_config(args.config)

    data_ok = False
    if config_ok:
        data_ok = check_input_data(config)
        check_output_directories(config)

    deps_ok = check_dependencies()
    print_summary(config_ok, data_ok, deps_ok)


if __name__ == "__main__":
    main()
