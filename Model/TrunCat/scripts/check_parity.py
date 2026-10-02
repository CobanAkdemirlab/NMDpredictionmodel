#!/usr/bin/env python3
"""
check_parity.py
===============
Checks that the standalone scripts (01-03) reproduce the committed TrunCat results.

The scripts run in a scratch directory, so nothing in the repo is overwritten. Their outputs
are then compared with the committed files:

  stage 01  merged table          same columns, same column order, same values
  stage 02  cleaned table, final_feature_list.csv, training_medians.json
            same values, same feature list and order, same categorical flags
  stage 03  out-of-fold AUC, mean fold AUC, per-fold AUCs, importance ranking, CV predictions
            equal within a tolerance (CatBoost results vary slightly with library versions
            and thread counts, so exact equality is not expected)

Usage (from anywhere):
    python Model/TrunCat/scripts/check_parity.py                 # stages 01-02, fast and exact
    python Model/TrunCat/scripts/check_parity.py --through 03    # also trains the model (slow; includes SHAP)
    python Model/TrunCat/scripts/check_parity.py --workdir /tmp/parity --no-run --through 03
                                                                 # compare outputs from an earlier run

Stage 01 needs the annotation files in data/annotations/ (git-ignored). If they are missing, stage 01 is
skipped and stage 02 starts from the committed data/TOPMed_merged.csv.

Exit code 0 = every check passed, 1 = at least one check failed.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

SCRIPTS = Path(__file__).resolve().parent
TRUNCAT = SCRIPTS.parent
INPUT_KEYS = ["enhanced_features", "codon_optimality", "readthrough", "ejc", "ptc_aug", "conservation"]
ANNOTATION_KEYS = INPUT_KEYS[1:]
STAGES = {"01": "01_data_loading_and_merging.py",
          "02": "02_feature_cleaning_and_selection.py",
          "03": "03_model_training.py"}

RESULTS = []


def record(ok, name, detail=""):
    RESULTS.append((ok, name, detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))


def versions():
    out = {"python": sys.version.split()[0]}
    for mod in ("pandas", "numpy", "sklearn", "catboost", "shap"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:
            out[mod] = "not installed"
    return out


# ------------------------------------------------------------------------------ scratch config
def committed(rel):
    return (TRUNCAT / rel).resolve()


def make_scratch_config(workdir, copy_merged):
    cfg = yaml.safe_load(open(TRUNCAT / "config" / "config.yaml"))
    for key in INPUT_KEYS:                       # inputs stay in the repo (absolute paths)
        cfg["data"][key] = str(committed(cfg["data"][key]))
    (workdir / "config").mkdir(parents=True, exist_ok=True)
    shutil.copy(TRUNCAT / "config" / "feature_column_order.txt", workdir / "config")
    path = workdir / "config" / "config.yaml"    # outputs are relative -> land in the scratch dir
    with open(path, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    if copy_merged:                              # stage 01 skipped: start stage 02 from the committed merged table
        dst = workdir / cfg["data"]["merged"]
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(committed("data/TOPMed_merged.csv"), dst)
    return path, cfg


def run_stage(stage, cfg_path, workdir):
    log = workdir / "logs" / f"{stage}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(log, "w") as f:
        proc = subprocess.run([sys.executable, str(SCRIPTS / STAGES[stage]), "--config", str(cfg_path)],
                              stdout=f, stderr=subprocess.STDOUT, env={**os.environ, "MPLBACKEND": "Agg"})
    ok = proc.returncode == 0
    record(ok, f"script {stage} ran to completion", f"{time.time() - t0:.0f} s, log: {log}")
    if not ok:
        print("\n".join("      " + l for l in log.read_text().splitlines()[-12:]))
    return ok


# ------------------------------------------------------------------------------ comparisons
def same_table(a, b, label):
    if list(a.columns) != list(b.columns):
        only_a = [c for c in a.columns if c not in set(b.columns)][:3]
        only_b = [c for c in b.columns if c not in set(a.columns)][:3]
        pos = next((i for i, (x, y) in enumerate(zip(a.columns, b.columns)) if x != y), None)
        record(False, f"{label}: columns and column order", f"first order difference at {pos}; only committed {only_a}; only script {only_b}")
        return
    record(True, f"{label}: columns and column order", f"{a.shape[1]} columns")
    if a.shape != b.shape:
        record(False, f"{label}: shape", f"{a.shape} vs {b.shape}")
        return
    bad = []
    for c in a.columns:
        x, y = a[c], b[c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            ok = np.isclose(x.astype(float), y.astype(float), rtol=1e-9, atol=1e-12, equal_nan=True).all()
        else:
            ok = (x.where(x.notna(), "<NA>").astype(str).values == y.where(y.notna(), "<NA>").astype(str).values).all()
        if not ok:
            bad.append(c)
    record(not bad, f"{label}: values", "all columns match" if not bad else f"{len(bad)} columns differ, e.g. {bad[:5]}")


def read_csv(path):
    return pd.read_csv(path, low_memory=False)


def compare_stage01(workdir, cfg):
    print("\nStage 01 (merged table)")
    same_table(read_csv(committed("data/TOPMed_merged.csv")), read_csv(workdir / cfg["data"]["merged"]), "merged table")


def compare_stage02(workdir, cfg):
    print("\nStage 02 (cleaning)")
    same_table(read_csv(committed("data/TOPMed_cleaned.csv")), read_csv(workdir / cfg["data"]["cleaned"]), "cleaned table")
    a, b = read_csv(committed("data/final_feature_list.csv")), read_csv(workdir / cfg["data"]["feature_list"])
    record(a.feature.tolist() == b.feature.tolist(), "final_feature_list: feature names and order", f"{len(a)} vs {len(b)} features")
    if a.feature.tolist() == b.feature.tolist():
        for col in ("dtype", "is_categorical", "non_null_count"):
            record((a[col].astype(str) == b[col].astype(str)).all(), f"final_feature_list: {col}")
    ma = json.load(open(committed("predict/training_medians.json")))
    mb = json.load(open(workdir / "predict" / "training_medians.json"))
    record(ma.keys() == mb.keys() and all(abs(ma[k] - mb[k]) < 1e-9 for k in ma), "training_medians.json", str(mb))


def parse_summary(path):
    txt = Path(path).read_text()
    oof = float(re.search(r"Out-of-Fold AUC:\s*([0-9.]+)", txt).group(1))
    mean = float(re.search(r"Mean Fold AUC:\s*([0-9.]+)", txt).group(1))
    folds = [float(x) for x in re.findall(r"([0-9]\.[0-9]+)", re.search(r"Fold AUCs:\s*(\[.*?\])", txt).group(1))]
    return oof, mean, folds


def compare_stage03(workdir, tol):
    print("\nStage 03 (model training; tolerances, not exact equality)")
    oa, ma, fa = parse_summary(committed("results/model_performance_summary.txt"))
    ob, mb, fb = parse_summary(workdir / "results" / "model_performance_summary.txt")
    record(abs(oa - ob) <= tol, "out-of-fold AUC", f"committed {oa:.4f}, script {ob:.4f}, tolerance {tol}")
    record(abs(ma - mb) <= tol, "mean fold AUC", f"committed {ma:.4f}, script {mb:.4f}")
    worst = max(abs(x - y) for x, y in zip(fa, fb))
    record(len(fa) == len(fb) and worst <= max(0.005, 2 * tol), "per-fold AUCs", f"largest difference {worst:.4f}")
    ca, cb = read_csv(committed("results/cv_predictions.csv")), read_csv(workdir / "results" / "cv_predictions.csv")
    record(ca["key"].tolist() == cb["key"].tolist() and (ca["true_label"] == cb["true_label"]).all(), "cv_predictions: same variants, same order, same labels")
    r = float(np.corrcoef(ca["predicted_prob"], cb["predicted_prob"])[0, 1])
    record(r >= 0.95, "cv_predictions: correlation of probabilities", f"r = {r:.4f}")
    ia = read_csv(committed("results/feature_importances_cv_averaged.csv")).sort_values("importance_mean", ascending=False)
    ib = read_csv(workdir / "results" / "feature_importances_cv_averaged.csv").sort_values("importance_mean", ascending=False)
    overlap = len(set(ia.feature.head(10)) & set(ib.feature.head(10)))
    record(overlap >= 8, "feature importance: top-10 overlap", f"{overlap} of 10")
    record((workdir / "model" / "truncat.cbm").exists(), "final model file written")


# ------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description="Check that scripts 01-03 reproduce the committed TrunCat results.")
    ap.add_argument("--through", choices=["01", "02", "03"], default="02", help="last stage to run and check (default 02)")
    ap.add_argument("--workdir", help="scratch directory (default: a temporary one, deleted afterwards)")
    ap.add_argument("--keep", action="store_true", help="keep the scratch directory")
    ap.add_argument("--no-run", action="store_true", help="do not run the scripts; compare outputs already in --workdir")
    ap.add_argument("--skip-01", action="store_true", help="start stage 02 from the committed merged table")
    ap.add_argument("--tol-auc", type=float, default=0.003, help="AUC tolerance for stage 03 (default 0.003)")
    args = ap.parse_args()
    if args.no_run and not args.workdir:
        ap.error("--no-run needs --workdir")

    print("Environment:", versions())
    try:
        major, minor = (int(x) for x in __import__("sklearn").__version__.split(".")[:2])
        if (major, minor) >= (1, 8):
            print("WARNING: scikit-learn >= 1.8 builds different gene-grouped CV folds; stage 03 will not match. Use scikit-learn<1.8.")
    except Exception:
        pass

    stages = [s for s in STAGES if s <= args.through]
    annotations_ok = all(committed(yaml.safe_load(open(TRUNCAT / "config" / "config.yaml"))["data"][k]).exists() for k in INPUT_KEYS)
    run_01 = "01" in stages and annotations_ok and not args.skip_01
    if "01" in stages and not run_01:
        print("Stage 01 skipped (annotation files missing or --skip-01): stage 02 starts from the committed merged table.")
    workdir = Path(args.workdir).resolve() if args.workdir else Path(tempfile.mkdtemp(prefix="truncat_parity_"))
    workdir.mkdir(parents=True, exist_ok=True)
    print("Scratch directory:", workdir)

    if args.no_run:
        cfg = yaml.safe_load(open(workdir / "config" / "config.yaml"))
        cfg_path = workdir / "config" / "config.yaml"
    else:
        cfg_path, cfg = make_scratch_config(workdir, copy_merged=not run_01)
        for stage in stages:
            if stage == "01" and not run_01:
                continue
            print(f"\nRunning script {stage} ...")
            if not run_stage(stage, cfg_path, workdir):
                break

    if run_01 or (args.no_run and "01" in stages and (workdir / cfg["data"]["merged"]).exists() and not args.skip_01 and annotations_ok):
        compare_stage01(workdir, cfg)
    if "02" in stages and (workdir / cfg["data"]["cleaned"]).exists():
        compare_stage02(workdir, cfg)
    if "03" in stages and (workdir / "results" / "model_performance_summary.txt").exists():
        compare_stage03(workdir, args.tol_auc)

    n_fail = sum(1 for ok, _, _ in RESULTS if not ok)
    print(f"\n{'=' * 70}\n{len(RESULTS) - n_fail} passed, {n_fail} failed")
    if n_fail:
        print("Logs and outputs kept in", workdir)
    elif not args.workdir and not args.keep:
        shutil.rmtree(workdir, ignore_errors=True)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
