"""
Reusable utilities for benchmarking TrunCat / TrunKitten against NMDetective-B.

Kept as a module rather than inlined in the notebook so:
- The notebook stays readable (markdown + a few cells per step)
- The metrics functions are unit-testable and reusable
- Future benchmarks can import the same helpers

All functions assume:
- y_true: 1 = escape, 0 = NMD-sensitive
- All score columns: higher = more escape

Bootstrap CIs can be gene-clustered: pass `groups` (one gene ID per row) and
whole genes are resampled instead of individual variants. This matches the
gene-grouped cross-validation used to train the models.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import (
    roc_auc_score,
    roc_curve,
    precision_recall_curve,
    average_precision_score,
    confusion_matrix,
    brier_score_loss,
    f1_score,
    balanced_accuracy_score,
)


# ----------------------------------------------------------------------
# Bootstrap CIs
# ----------------------------------------------------------------------

def bootstrap_metric(
    y_true: np.ndarray,
    y_score: np.ndarray,
    metric_fn,
    n_boot: int = 1000,
    seed: int = 42,
    ci: float = 0.95,
    groups: Optional[np.ndarray] = None,
) -> Tuple[float, float, float]:
    """
    Bootstrap a metric. Returns (point_estimate, lower, upper).

    groups=None  -> resample individual variants.
    groups given -> resample whole groups (genes) with replacement, keeping all
                    variants of each drawn gene. Use this when variants within
                    a gene are correlated (they are here).

    The RNG depends only on `seed` and the group structure, so calling this
    repeatedly with the same seed and the same rows gives the same resamples
    for every model.
    """
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    point = metric_fn(y_true, y_score)

    rng = np.random.default_rng(seed)

    if groups is None:
        n = len(y_true)

        def draw():
            return rng.integers(0, n, size=n)
    else:
        groups = np.asarray(groups)
        codes, uniq = pd.factorize(groups)
        if (codes < 0).any():
            raise ValueError("groups contains NaN; every row needs a gene ID")
        order = np.argsort(codes, kind="stable")
        bounds = np.searchsorted(codes[order], np.arange(len(uniq) + 1))
        by_group = [order[bounds[g]:bounds[g + 1]] for g in range(len(uniq))]
        n_groups = len(uniq)

        def draw():
            picked = rng.integers(0, n_groups, size=n_groups)
            return np.concatenate([by_group[g] for g in picked])

    boot_vals = np.empty(n_boot, dtype=np.float64)
    for i in range(n_boot):
        idx = draw()
        try:
            boot_vals[i] = metric_fn(y_true[idx], y_score[idx])
        except ValueError:
            # Resample lacks both classes; mark NaN
            boot_vals[i] = np.nan

    valid = boot_vals[~np.isnan(boot_vals)]
    alpha = (1 - ci) / 2
    lo, hi = np.quantile(valid, [alpha, 1 - alpha])
    return float(point), float(lo), float(hi)


def auc_with_ci(y_true, y_score, n_boot=1000, seed=42, groups=None):
    return bootstrap_metric(y_true, y_score, roc_auc_score, n_boot, seed,
                            groups=groups)


def pr_auc_with_ci(y_true, y_score, n_boot=1000, seed=42, groups=None):
    return bootstrap_metric(y_true, y_score, average_precision_score, n_boot,
                            seed, groups=groups)


# ----------------------------------------------------------------------
# Threshold helpers
# ----------------------------------------------------------------------

def youden_threshold(y_true, y_score) -> float:
    """Threshold maximizing Youden's J = sensitivity + specificity − 1."""
    fpr, tpr, thr = roc_curve(y_true, y_score)
    j = tpr - fpr
    return float(thr[np.argmax(j)])


def threshold_at_specificity(y_true, y_score, target_specificity: float) -> float:
    """Threshold with the highest sensitivity among those achieving >= target specificity."""
    fpr, tpr, thr = roc_curve(y_true, y_score)
    spec = 1 - fpr
    # Find thresholds where spec >= target; pick the one with highest sensitivity
    mask = spec >= target_specificity
    if not mask.any():
        return float(thr[0])  # fallback: most stringent
    valid_thr = thr[mask]
    valid_tpr = tpr[mask]
    return float(valid_thr[np.argmax(valid_tpr)])


def threshold_at_sensitivity(y_true, y_score, target_sensitivity: float) -> float:
    fpr, tpr, thr = roc_curve(y_true, y_score)
    mask = tpr >= target_sensitivity
    if not mask.any():
        return float(thr[-1])
    valid_thr = thr[mask]
    valid_fpr = fpr[mask]
    return float(valid_thr[np.argmin(valid_fpr)])


# ----------------------------------------------------------------------
# Binary classification metric block
# ----------------------------------------------------------------------

def binary_metrics(
    y_true: np.ndarray,
    y_score: np.ndarray,
    threshold: float,
) -> Dict[str, float]:
    """
    Compute all binary metrics at a given threshold.
    y_pred = 1 (escape) if y_score >= threshold else 0.
    """
    y_true = np.asarray(y_true).astype(int)
    y_pred = (np.asarray(y_score) >= threshold).astype(int)

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    sens = tp / (tp + fn) if (tp + fn) else np.nan       # recall for escape
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    ppv  = tp / (tp + fp) if (tp + fp) else np.nan
    npv  = tn / (tn + fn) if (tn + fn) else np.nan
    acc  = (tp + tn) / (tp + tn + fp + fn)
    f1   = f1_score(y_true, y_pred, zero_division=0)
    bal  = balanced_accuracy_score(y_true, y_pred)

    return {
        "threshold":          float(threshold),
        "TP": int(tp), "FP": int(fp), "TN": int(tn), "FN": int(fn),
        "sensitivity":        float(sens),
        "specificity":        float(spec),
        "ppv":                float(ppv),
        "npv":                float(npv),
        "accuracy":           float(acc),
        "f1":                 float(f1),
        "balanced_accuracy":  float(bal),
    }


# ----------------------------------------------------------------------
# Continuous metric block
# ----------------------------------------------------------------------

def continuous_metrics(
    y_true: np.ndarray,
    y_score: np.ndarray,
    n_boot: int = 1000,
    seed: int = 42,
    groups: Optional[np.ndarray] = None,
) -> Dict[str, float]:
    """ROC-AUC, PR-AUC (with bootstrap CIs), and Brier score."""
    auc_p, auc_lo, auc_hi = auc_with_ci(y_true, y_score, n_boot, seed, groups)
    pr_p,  pr_lo,  pr_hi  = pr_auc_with_ci(y_true, y_score, n_boot, seed, groups)

    # Brier score requires probabilities in [0,1]. NMDetective-B's escape score
    # (1 − raw) IS in [0, 1] because raw is in [0, 0.65]. So Brier is defined,
    # but interpret with care: NMDetective-B is not probability-calibrated.
    if np.all((y_score >= 0) & (y_score <= 1)):
        brier = float(brier_score_loss(y_true, y_score))
    else:
        brier = np.nan

    return {
        "roc_auc":          auc_p,
        "roc_auc_lo":       auc_lo,
        "roc_auc_hi":       auc_hi,
        "pr_auc":           pr_p,
        "pr_auc_lo":        pr_lo,
        "pr_auc_hi":        pr_hi,
        "brier_score":      brier,
    }


# ----------------------------------------------------------------------
# Within-bin stratified analysis
# ----------------------------------------------------------------------

def within_bin_auc(
    df: pd.DataFrame,
    bin_col: str,
    score_col: str,
    label_col: str,
    min_per_class: int = 10,
    group_col: Optional[str] = None,
    n_boot: int = 1000,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Compute ROC-AUC and PR-AUC of `score_col` against `label_col` within each
    level of `bin_col`. Skips bins with fewer than min_per_class of either
    class — AUC is undefined or noisy with few examples.

    If group_col is given (e.g. "GENE_ID"), bootstrap CIs resample whole genes
    within each bin.
    """
    rows = []
    for bin_val, sub in df.groupby(bin_col, dropna=False):
        n = len(sub)
        n_pos = int(sub[label_col].sum())
        n_neg = n - n_pos
        rec = {
            "bin":       bin_val,
            "n":         n,
            "n_escape":  n_pos,
            "n_sensitive": n_neg,
        }
        if n_pos < min_per_class or n_neg < min_per_class:
            rec.update({
                "roc_auc": np.nan, "roc_auc_lo": np.nan, "roc_auc_hi": np.nan,
                "pr_auc":  np.nan, "pr_auc_lo":  np.nan, "pr_auc_hi":  np.nan,
                "skipped_reason": f"need >={min_per_class} per class",
            })
        else:
            g = sub[group_col].values if group_col else None
            auc_p, auc_lo, auc_hi = auc_with_ci(
                sub[label_col].values, sub[score_col].values,
                n_boot=n_boot, seed=seed, groups=g,
            )
            pr_p, pr_lo, pr_hi = pr_auc_with_ci(
                sub[label_col].values, sub[score_col].values,
                n_boot=n_boot, seed=seed, groups=g,
            )
            rec.update({
                "roc_auc": auc_p, "roc_auc_lo": auc_lo, "roc_auc_hi": auc_hi,
                "pr_auc":  pr_p,  "pr_auc_lo":  pr_lo,  "pr_auc_hi":  pr_hi,
                "skipped_reason": "",
            })
        rows.append(rec)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------

# Categorical color scheme; matches your manuscript visualization standards
MODEL_COLORS = {
    "TrunCat":       "#2E5C8A",   # deep blue
    "TrunKitten":    "#5BA8D4",   # lighter blue (refined-from-TrunCat feel)
    "NMDetective-B": "#D17A22",   # orange — distinct family
}


def plot_roc_overlay(
    results: Dict[str, dict],   # name -> {y_true, y_score, legend}
    title: str = "",
    figsize=(7, 7),
    save_path: Optional[Path] = None,
):
    """
    Overlay ROC curves.

    NMDetective-B has only 5 distinct score values, so its ROC has only a few
    operating points. It is drawn with straight segments (matching how the
    tie-aware trapezoidal AUC is computed) plus markers at the actual points.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    for name, d in results.items():
        fpr, tpr, _ = roc_curve(d["y_true"], d["y_score"])
        color = MODEL_COLORS.get(name, "#444")
        if name == "NMDetective-B":
            ax.plot(fpr, tpr, label=d["legend"], color=color, linewidth=2.2,
                    marker="o", markersize=6)
        else:
            ax.plot(fpr, tpr, label=d["legend"], color=color, linewidth=2.2)
    ax.plot([0, 1], [0, 1], color="gray", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xlabel("False positive rate (1 − specificity)", fontsize=12)
    ax.set_ylabel("True positive rate (sensitivity)", fontsize=12)
    if title:
        ax.set_title(title, fontsize=13)
    ax.legend(loc="lower right", fontsize=11, frameon=False)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    if save_path:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        fig.savefig(save_path.with_suffix(".pdf"), bbox_inches="tight")
    return fig, ax

def plot_pr_overlay(
    results: Dict[str, dict],
    title: str = "",
    figsize=(7, 7),
    save_path: Optional[Path] = None,
):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    for name, d in results.items():
        prec, rec, _ = precision_recall_curve(d["y_true"], d["y_score"])
        color = MODEL_COLORS.get(name, "#444")
        if name == "NMDetective-B":
            # Step interpolation matches how average precision is defined
            ax.plot(rec, prec, label=d["legend"], color=color, linewidth=2.2,
                    drawstyle="steps-post")
        else:
            # Drop the recall=0 points so a false positive at the top of the
            # ranking doesn't draw a vertical line down the y-axis
            keep = rec > 0
            ax.plot(rec[keep], prec[keep], label=d["legend"], color=color,
                    linewidth=2.2)
    # Baseline = prevalence of positive class
    if results:
        any_y = next(iter(results.values()))["y_true"]
        prevalence = float(np.mean(any_y))
        ax.axhline(prevalence, color="gray", linestyle="--", linewidth=1,
                   alpha=0.6, label=f"Prevalence = {prevalence:.2f}")
    ax.set_xlabel("Recall (sensitivity)", fontsize=12)
    ax.set_ylabel("Precision (PPV)", fontsize=12)
    if title:
        ax.set_title(title, fontsize=13)
    ax.legend(loc="lower left", fontsize=11, frameon=False)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    if save_path:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        fig.savefig(save_path.with_suffix(".pdf"), bbox_inches="tight")
    return fig, ax

def plot_within_bin_strip(
    df: pd.DataFrame,
    bin_col: str,
    score_col: str,
    label_col: str,
    bin_order: Optional[List[str]] = None,
    title: str = "",
    figsize=(10, 6),
    save_path: Optional[Path] = None,
):
    """
    Box+strip plot: TrunCat/TrunKitten escape probability stratified by
    NMDetective-B rule bin, with empirical escape rate overlay.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    if bin_order is None:
        bin_order = list(df[bin_col].unique())

    fig, ax = plt.subplots(figsize=figsize, dpi=150)

    # Box (no whisker outliers; we'll add strip below)
    sns.boxplot(
        data=df, x=bin_col, y=score_col, order=bin_order,
        ax=ax, color="#dde6ed", fliersize=0, width=0.55,
    )
    # Strip colored by escape label
    sns.stripplot(
        data=df, x=bin_col, y=score_col, hue=label_col, order=bin_order,
        palette={1: "#2E5C8A", 0: "#D17A22"},
        ax=ax, alpha=0.35, jitter=0.25, size=2.5, dodge=False,
    )
    # Empirical escape rate per bin overlaid as black diamonds
    rates = df.groupby(bin_col)[label_col].mean().reindex(bin_order)
    ax.scatter(
        range(len(bin_order)), rates.values,
        marker="D", s=70, color="black", zorder=5, label="Empirical escape rate",
    )
    ax.set_xlabel("NMDetective-B rule bin", fontsize=12)
    ax.set_ylabel(f"Predicted escape probability ({score_col})", fontsize=12)
    if title:
        ax.set_title(title, fontsize=13)
    handles, labels = ax.get_legend_handles_labels()
    # Relabel legend ints -> human names
    new_labels = ["NMD-sensitive" if l == "0" else "Escape" if l == "1" else l
                  for l in labels]
    ax.legend(handles, new_labels, loc="best", fontsize=10, frameon=False)
    ax.grid(True, axis="y", alpha=0.25)
    plt.xticks(rotation=20, ha="right")
    fig.tight_layout()
    if save_path:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        fig.savefig(save_path.with_suffix(".pdf"), bbox_inches="tight")
    return fig, ax