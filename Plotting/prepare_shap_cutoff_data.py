"""
Build the three CSVs used by Panel G of NMD_Figure4_main_v7.R from the
OOF SHAP pickle. Uses the same method and settings as shap_cutoff_analysis.ipynb
(LOESS frac = 0.10 on a 500-point grid; rate-of-change peaks with prominence
>= 0.30 x max |derivative|).

Usage:  python prepare_shap_cutoff_data.py [path/to/shap_values_for_sharing.pkl]
Writes: shap_cutoff_curve_data.csv, shap_cutoff_loess.csv, shap_cutoff_markers.csv
"""
import sys
import pickle
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.signal import find_peaks

PKL = sys.argv[1] if len(sys.argv) > 1 else "shap_values_for_sharing.pkl"
FEATURE = "relativePTClocation"
LOESS_FRAC = 0.10
GRID_POINTS = 500
PEAK_PROMINENCE_FRAC = 0.30

with open(PKL, "rb") as f:
    d = pickle.load(f)

X = d["feature_data"]
y = d["shap_values"][:, X.columns.get_loc(FEATURE)]
x = pd.to_numeric(X[FEATURE], errors="coerce").values
ok = ~np.isnan(x)
x, y = x[ok], y[ok]

# LOESS smooth on a regular grid
sm_xy = sm.nonparametric.lowess(y, x, frac=LOESS_FRAC, return_sorted=True)
grid = np.linspace(x.min(), x.max(), GRID_POINTS)
smooth = np.interp(grid, sm_xy[:, 0], sm_xy[:, 1])
dy = np.gradient(smooth, grid)

# Zero crossings of the smooth
markers = []
for i in range(len(smooth) - 1):
    y1, y2 = smooth[i], smooth[i + 1]
    if y1 * y2 < 0:
        zx = grid[i] - y1 * (grid[i + 1] - grid[i]) / (y2 - y1)
        markers.append({"marker_type": "zero_crossing", "x": zx,
                        "label": "pos_to_neg" if y1 > 0 else "neg_to_pos"})

# Rate-of-change peaks (steepest transitions)
peaks, _ = find_peaks(np.abs(dy), prominence=PEAK_PROMINENCE_FRAC * np.abs(dy).max())
for p in peaks:
    markers.append({"marker_type": "rate_peak", "x": grid[p],
                    "label": "increasing" if dy[p] > 0 else "decreasing"})

pd.DataFrame({FEATURE: x, "shap_value": y}).to_csv("shap_cutoff_curve_data.csv", index=False)
pd.DataFrame({"x_grid": grid, "y_smooth": smooth, "dy_dx": dy}).to_csv("shap_cutoff_loess.csv", index=False)
pd.DataFrame(markers).to_csv("shap_cutoff_markers.csv", index=False)

for m in markers:
    print(f"{m['marker_type']:14s} {m['label']:11s} {m['x']:.4f}")
print(f"n variants = {len(x)}")
