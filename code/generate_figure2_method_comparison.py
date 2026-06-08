#!/usr/bin/env python3
"""Generate Figure 2E: five-method AUROC comparison on GSE26050.

The script first looks for the optional archived input
analysis_p1_R/E1_real_data.csv. When that archive is not distributed, it uses
the AUROC values reported in the manuscript to regenerate the submitted panel.
"""

import numpy as np
import pandas as pd
import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, ".."))
FIG_DIR = os.path.join(PROJECT_ROOT, "figures")
PNG_DIR = os.path.join(PROJECT_ROOT, "manuscript", "docx_embedded_figures_png")

data_file = os.path.join(PROJECT_ROOT, "analysis_p1_R", "E1_real_data.csv")
if os.path.exists(data_file):
    df = pd.read_csv(data_file)
else:
    # Keep this figure reproducible from the submission package even when
    # the archived analysis_p1_R directory is not distributed.
    df = pd.DataFrame({
        "Method": ["CADE", "Unadjusted", "SVA-like", "RUVg-like", "PCA top-4"],
        "AUROC": [0.886, 1.000, 1.000, 0.705, 0.698],
    })

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "figure.dpi": 300,
})

fig, ax = plt.subplots(figsize=(5.5, 4))

colors = {
    "CADE": "#B2182B",
    "Unadjusted": "#7F7F7F",
    "PCA top-4": "#F4A582",
    "SVA-like": "#2166AC",
    "RUVg-like": "#92C5DE",
}

methods = df["Method"].tolist()
aurocs = df["AUROC"].tolist()
bar_colors = [colors.get(m, "#7F7F7F") for m in methods]

bars = ax.barh(methods, aurocs, color=bar_colors, edgecolor="white", linewidth=0.5)

for bar, val in zip(bars, aurocs):
    if val >= 0.96:
        ax.text(val - 0.025, bar.get_y() + bar.get_height()/2,
                f"{val:.3f}", va="center", ha="right", fontsize=9, color="white",
                fontweight="bold")
    else:
        ax.text(val + 0.02, bar.get_y() + bar.get_height()/2,
                f"{val:.3f}", va="center", ha="left", fontsize=9)

ax.set_xlim(0, 1.15)
ax.set_xlabel("AUROC (cell-intrinsic DE detection)")
ax.set_title("E. Five-method comparison on real GSE26050", fontweight="bold")
ax.axvline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.3)
plt.tight_layout()
out = os.path.join(FIG_DIR, "Figure2E_5method_AUROC.tif")
plt.savefig(out, dpi=300, bbox_inches="tight", format="tif")
os.makedirs(PNG_DIR, exist_ok=True)
png_out = os.path.join(PNG_DIR, "Figure2E_5method_AUROC.png")
plt.savefig(png_out, dpi=300, bbox_inches="tight", format="png")
print(f"Saved Figure 2E to: {out}")
print(f"Saved Figure 2E PNG to: {png_out}")
plt.close()
