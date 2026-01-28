#!/usr/bin/env python3
"""
Setup Script for NMD Escape Prediction Model Scripts
====================================================
This script helps verify the setup when running from the scripts/ folder.
It assumes config/config.yaml already exists from the notebooks setup.

Usage:
    python setup.py
"""

import os
from pathlib import Path


def check_config():
    """Check if config file exists."""
    print("Checking configuration...")
    
    config_file = Path("../config/config.yaml")
    
    if config_file.exists():
        print(f"  ✓ Found config file: {config_file}")
        print(f"  ℹ️  Using shared configuration with notebooks")
    else:
        print(f"  ✗ Config file not found at: {config_file}")
        print(f"  ⚠️  Please set up config/config.yaml first")
        print(f"      (See notebooks setup or use config_example.yaml as template)")
        return False
    
    return True


def check_output_directories():
    """Check/create output directories from scripts/ folder."""
    print("\nChecking output directories...")
    
    # These paths are relative to the CatBoost/ root, not scripts/
    base_dirs = [
        "../data/processed",
        "../output/models",
        "../output/models/cv_folds",
        "../output/results",
        "../output/figures"
    ]
    
    for directory in base_dirs:
        dir_path = Path(directory)
        if not dir_path.exists():
            dir_path.mkdir(parents=True, exist_ok=True)
            print(f"  ✓ Created: {directory}")
        else:
            print(f"  ✓ Exists: {directory}")


def check_dependencies():
    """Check if required packages are installed."""
    print("\nChecking dependencies...")
    
    required_packages = [
        ("pandas", "pandas"),
        ("numpy", "numpy"),
        ("yaml", "pyyaml"),
        ("catboost", "catboost"),
        ("sklearn", "scikit-learn"),
        ("matplotlib", "matplotlib"),
        ("seaborn", "seaborn"),
        ("shap", "shap"),
        ("tqdm", "tqdm"),
        ("joblib", "joblib")
    ]
    
    missing = []
    for import_name, package_name in required_packages:
        try:
            __import__(import_name)
            print(f"  ✓ {package_name}")
        except ImportError:
            print(f"  ✗ {package_name} (missing)")
            missing.append(package_name)
    
    if missing:
        print(f"\n⚠️  Missing packages: {', '.join(missing)}")
        print(f"   Install with: pip install -r requirements.txt")
        return False
    else:
        print("\n✓ All dependencies installed!")
        return True


def print_next_steps(config_ok, deps_ok):
    """Print next steps for user."""
    print("\n" + "=" * 80)
    
    if config_ok and deps_ok:
        print("SETUP VERIFIED - READY TO RUN!")
        print("=" * 80)
        print("""
You can now run the scripts:

    cd scripts/  # (if not already here)
    python 01_data_loading_and_merging.py
    python 02_feature_cleaning_and_selection.py
    python 03_model_training.py

The scripts will use the configuration from ../config/config.yaml
        """)
    else:
        print("SETUP INCOMPLETE")
        print("=" * 80)
        print("\nPlease fix the issues above before running the scripts.")
        
        if not config_ok:
            print("\n1. Set up ../config/config.yaml with your data paths")
        if not deps_ok:
            print("\n2. Install missing dependencies: pip install -r requirements.txt")


def main():
    """Main setup function."""
    print("=" * 80)
    print("NMD Escape Prediction Model - Scripts Setup Verification")
    print("=" * 80)
    print()
    
    config_ok = check_config()
    check_output_directories()
    deps_ok = check_dependencies()
    print_next_steps(config_ok, deps_ok)


if __name__ == "__main__":
    main()
