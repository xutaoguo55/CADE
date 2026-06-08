#!/usr/bin/env python3
"""Regenerate Figure 2 with separated panel titles and Word-safe spacing."""

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PNG_OUT = ROOT / "manuscript" / "docx_embedded_figures_png" / "Figure2_CADE_Benchmark.png"
TIF_OUT = ROOT / "figures" / "Figure2_CADE_Benchmark.tif"


def read_captioned_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, skiprows=1)


def bias_to_numeric(series: pd.Series) -> np.ndarray:
    return series.astype(str).str.extract(r"([0-9]+)").astype(float)[0].to_numpy()


def panel_label(ax, label: str, title: str) -> None:
    ax.text(-0.08, 1.08, label, transform=ax.transAxes, fontsize=13,
            weight="bold", va="bottom", ha="left")
    ax.set_title("\n".join(textwrap.wrap(title, width=46)), fontsize=12,
                 weight="bold", loc="left", pad=10)


def main() -> None:
    table_a = read_captioned_csv(ROOT / "tables" / "Table_01A_Benchmark_Gradient.csv")
    table_b = read_captioned_csv(ROOT / "tables" / "Table_01B_scRNA_Pseudobulk_Benchmark.csv")
    table_c = read_captioned_csv(ROOT / "tables" / "Table_01C_MAS_Pseudobulk_Validation.csv")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.28,
        "grid.linewidth": 0.8,
    })

    green = "#0b7f35"
    red = "#d35f77"
    purple = "#3d2c8d"
    yellow = "#dfc766"
    gray = "#8c8c8c"

    fig, axes = plt.subplots(2, 2, figsize=(15, 10.5), dpi=300)
    fig.suptitle(
        "Figure 2 | CADE benchmark validation across composition-bias gradients and ground-truth datasets",
        fontsize=16, weight="bold", y=0.985,
    )

    # Panel A
    ax = axes[0, 0]
    x = bias_to_numeric(table_a["Bias Level"])
    ax.plot(x, table_a["AUROC CADE"], "-o", color=green, lw=2.4, ms=7, label="CADE")
    ax.plot(x, table_a["AUROC Unadj"], "-s", color=red, lw=2.4, ms=7, label="Unadjusted limma")
    for xi, yi in zip(x, table_a["AUROC CADE"]):
        ax.text(xi, yi + 0.010, f"{yi:.3f}", color=green, fontsize=8.5, ha="center")
    for xi, yi in zip(x, table_a["AUROC Unadj"]):
        ax.text(xi, yi - 0.022, f"{yi:.3f}", color=red, fontsize=8.5, ha="center")
    ax.set_ylim(0.84, 1.025)
    ax.set_xlabel("Composition bias level (%)")
    ax.set_ylabel("AUROC for cell-intrinsic DE detection")
    ax.legend(frameon=True, loc="lower left", fontsize=9)
    panel_label(ax, "A", "Dose-response synthetic benchmark (3,000 genes; 44 samples; 5 cell types)")

    # Panel B
    ax = axes[0, 1]
    xb = bias_to_numeric(table_b["Bias Level"])
    series = [
        ("CADE", "CADE (Intrinsic)", green, "o"),
        ("limma", "limma (Intrinsic)", red, "s"),
        ("MarkerScore", "MarkerScore", purple, "^"),
        ("SVA", "SVA", yellow, "v"),
        ("Matched nnls*", "Matched nnls*", gray, "D"),
    ]
    for label, col, color, marker in series:
        ax.plot(xb, table_b[col], marker=marker, lw=2.0, ms=6.5, color=color, label=label)
    ax.set_ylim(0.50, 1.05)
    ax.set_xlim(-8, 158)
    ax.set_xlabel("Composition bias level (%)")
    ax.set_ylabel("AUROC for cell-intrinsic DE")
    ax.legend(frameon=True, loc="lower left", ncol=2, fontsize=8.5)
    ax.text(150, 0.52, "*matched reference,\nupper-bound estimate",
            ha="right", va="bottom", color="#777777", fontsize=8.5, style="italic")
    panel_label(ax, "B", "Sparse-marker scRNA-seq pseudobulk benchmark (449 genes; 18 markers; 60 samples)")

    # Panel C
    ax = axes[1, 0]
    true_de = table_a["CCI True DE"].to_numpy(dtype=float)
    bg = table_a["CCI Background"].to_numpy(dtype=float)
    ax.plot(x, true_de, "-o", color=green, lw=2.4, ms=7, label="True cell-intrinsic DE genes")
    ax.plot(x, bg, "-s", color=red, lw=2.4, ms=7, label="Background genes")
    ax.fill_between(x, true_de, bg, color=yellow, alpha=0.28, label="Separation")
    for xi, yi, sep in zip(x, (true_de + bg) / 2, table_a["CCI Separation"]):
        ax.text(xi, yi, f"Delta={sep:.2f}", color="#805700", fontsize=8.5,
                ha="center", va="center")
    ax.set_ylim(0.0, 0.84)
    ax.set_xlabel("Composition bias level (%)")
    ax.set_ylabel("Mean CCI")
    ax.legend(frameon=True, loc="upper left", fontsize=8.8)
    panel_label(ax, "C", "CCI noise-floor calibration (true DE vs background)")

    # Panel D
    ax = axes[1, 1]
    categories = ["B cells", "CD4 T cells", "CD8 T cells", "NK cells", "Monocytes",
                  "Macrophages", "Neutrophils", "Erythrocytes"]
    pbmc = np.array([1.00, 1.00, 0.99, 0.99, 0.99, 0.99, 0.98, 0.98])
    mas_lookup = dict(zip(table_c["CADE Cell Type"], table_c["Pearson r"]))
    mas = np.array([
        float(mas_lookup["B cells"]),
        float(mas_lookup["CD4 T cells"]),
        float(mas_lookup["CD8 T cells"]),
        float(mas_lookup["NK cells"]),
        float(mas_lookup["Monocytes"]),
        float(mas_lookup["Macrophages"]),
        np.nan,
        float(mas_lookup["Erythrocytes"]),
    ])
    idx = np.arange(len(categories))
    width = 0.36
    ax.bar(idx - width / 2, pbmc, width, color=green, label="PBMC synthetic ground truth")
    ax.bar(idx + width / 2, mas, width, color=purple, label="MAS pseudobulk vs scRNA-seq GT")
    for xi, val in zip(idx - width / 2, pbmc):
        ax.text(xi, val + 0.015, f"{val:.2f}", ha="center", va="bottom", fontsize=8.5)
    for xi, val in zip(idx + width / 2, mas):
        if np.isfinite(val):
            ax.text(xi, val + 0.015, f"{val:.2f}", ha="center", va="bottom", fontsize=8.5)
    ax.text(idx[6] + width / 2, 0.12, "no\nGT", ha="center", va="center",
            fontsize=8, color="#555555")
    ax.axhline(0.80, ls=":", lw=1.1, color="#999999")
    ax.text(len(categories) - 0.1, 0.82, "r = 0.80", ha="right", va="bottom",
            fontsize=9, color="#666666")
    ax.set_ylim(0.0, 1.08)
    ax.set_ylabel("Pearson r (estimated weights vs ground truth)")
    ax.set_xticks(idx)
    ax.set_xticklabels(categories, rotation=28, ha="right")
    ax.legend(frameon=True, loc="lower center", fontsize=8.5)
    panel_label(ax, "D", "Cell-type weight-estimation accuracy (two ground-truth datasets)")

    fig.subplots_adjust(left=0.07, right=0.985, top=0.90, bottom=0.09,
                        wspace=0.30, hspace=0.42)
    PNG_OUT.parent.mkdir(parents=True, exist_ok=True)
    TIF_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_OUT, dpi=300)
    fig.savefig(TIF_OUT, dpi=300)
    plt.close(fig)
    print(f"Wrote {PNG_OUT}")
    print(f"Wrote {TIF_OUT}")


if __name__ == "__main__":
    main()
