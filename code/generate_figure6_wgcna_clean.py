#!/usr/bin/env python3
"""Regenerate Figure 6 with readable WGCNA panels."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
PNG_OUT = ROOT / "manuscript" / "docx_embedded_figures_png" / "Figure6_WGCNA_Modules.png"
TIF_OUT = ROOT / "figures" / "Figure6_WGCNA_Modules.tif"


GENE_PANEL = {
    "SLC7A11", "TFRC", "FTH1", "FTL", "SLC25A37", "SLC40A1", "GPX4",
    "GCLM", "HMOX1", "NCOA4", "NFE2L2", "JAK2", "STAT1", "STAT3",
    "IFNG", "IL1B", "IL6", "TNF", "CXCL8",
}
IMMUNE_OR_DEATH = {
    "PYCARD", "WAS", "BID", "TICAM1", "FERMT3", "RHOG", "SBNO2",
    "MAP1S", "ARHGDIA", "PDLIM2",
}


def panel_label(ax, label: str, title: str) -> None:
    ax.text(-0.11, 1.06, label, transform=ax.transAxes, fontsize=13,
            weight="bold", va="bottom", ha="left")
    ax.set_title(title, fontsize=12, weight="bold", loc="left", pad=10)


def p_text(p: float) -> str:
    if p < 1e-4:
        return f"P={p:.1e}"
    if p < 0.001:
        return f"P={p:.3g}"
    return f"P={p:.3f}"


def gene_color(gene: str) -> str:
    if gene in GENE_PANEL:
        return "#c95d73"
    if gene in IMMUNE_OR_DEATH:
        return "#0b7f35"
    return "#80c8e8"


def radial_network(ax, hub: pd.DataFrame, module: int, corr: float) -> None:
    ax.set_axis_off()
    sub = hub[hub["Module"] == module].sort_values("Rank").head(15).copy()
    center_gene = sub.iloc[0]["Gene"]
    others = sub.iloc[1:]
    ax.scatter([0], [0], s=560, color="#3d2c8d", edgecolor="#111111", zorder=5)
    ax.text(0, 0, center_gene if len(center_gene) <= 7 else f"M{module}",
            ha="center", va="center", fontsize=8.5, color="white", weight="bold", zorder=6)
    angles = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, len(others), endpoint=False)
    radius = 1.0
    for angle, (_, row) in zip(angles, others.iterrows()):
        x, y = radius * np.cos(angle), radius * np.sin(angle)
        ax.plot([0, x], [0, y], color="#b8b8b8", lw=1.1, zorder=1)
        ax.scatter([x], [y], s=320, color=gene_color(row["Gene"]),
                   edgecolor="#111111", zorder=3)
        ha = "left" if x > 0.10 else "right" if x < -0.10 else "center"
        va = "bottom" if y > 0.65 else "top" if y < -0.65 else "center"
        ax.text(x * 1.17, y * 1.17, row["Gene"], fontsize=8.8, ha=ha, va=va)
    ax.set_xlim(-1.60, 1.60)
    ax.set_ylim(-1.45, 1.45)
    ax.text(0.02, 0.98, f"Module {module}: r={corr:+.2f}; top 15 by kME",
            transform=ax.transAxes, ha="left", va="top", fontsize=9.2)


def main() -> None:
    modules = pd.read_csv(ROOT / "supplementary" / "Table_S09_WGCNA_Module_FHL.csv")
    hub = pd.read_csv(ROOT / "supplementary" / "Table_S11_WGCNA_Hub_Genes.csv")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
    })
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 11), dpi=300)
    fig.suptitle("Figure 6 | WGCNA co-expression modules associated with FHL PBMC",
                 fontsize=16, weight="bold", y=0.985)

    # Panel A: all modules, sorted by correlation.
    ax = axes[0, 0]
    mod_sorted = modules.sort_values("Cor_FHL")
    y = np.arange(len(mod_sorted))
    colors = np.where(mod_sorted["Cor_FHL"] >= 0, "#c95d73", "#80c8e8")
    ax.barh(y, mod_sorted["Cor_FHL"], color=colors, edgecolor="#333333", height=0.70)
    ax.axvline(0, color="#333333", lw=0.8)
    for yi, (_, row) in zip(y, mod_sorted.iterrows()):
        x = row["Cor_FHL"]
        ha = "left" if x >= 0 else "right"
        offset = 0.025 if x >= 0 else -0.025
        stars = "***" if row["FDR"] < 0.001 else "**" if row["FDR"] < 0.01 else "*" if row["FDR"] < 0.05 else ""
        ax.text(x + offset, yi, f"{x:+.2f}{stars}",
                va="center", ha=ha, fontsize=8.2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"Module {int(m)}" for m in mod_sorted["Module"]], fontsize=8.8)
    ax.set_xlim(-1.05, 0.93)
    ax.set_xlabel("Pearson correlation with FHL status")
    panel_label(ax, "A", "Module-FHL status correlation (all 15 modules)")

    # Panel B: module size vs absolute correlation.
    ax = axes[0, 1]
    fdr = modules["FDR"].clip(lower=1e-300)
    neglog = -np.log10(fdr)
    size_scale = 45 + 520 * (neglog / neglog.max())
    sc = ax.scatter(modules["Size"], modules["Cor_FHL"].abs(), s=size_scale,
                    c=modules["Cor_FHL"], cmap="RdBu_r", vmin=-0.9, vmax=0.9,
                    edgecolor="#333333")
    label_offsets = {
        1: (18, 0.010), 2: (230, 0.015), 3: (22, 0.018), 4: (-45, 0.012),
        7: (240, 0.025), 8: (235, -0.015), 9: (38, 0.026),
    }
    for _, row in modules.iterrows():
        module = int(row["Module"])
        if module not in label_offsets:
            continue
        dx, dy = label_offsets[module]
        ha = "right" if dx < 0 else "left"
        ax.text(row["Size"] + dx, abs(row["Cor_FHL"]) + dy,
                f"M{module}", fontsize=8.7, weight="bold", ha=ha)
    ax.axhline(0.3, ls=":", color="#999999", lw=1)
    ax.axhline(0.5, ls="-", color="#bbbbbb", lw=0.9)
    ax.text(modules["Size"].max(), 0.515, "|r| = 0.50", ha="right",
            va="bottom", color="#777777", fontsize=9)
    ax.set_xlabel("Module size (number of genes)")
    ax.set_ylabel("|Pearson r| with FHL status")
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("Pearson r")
    ax.text(0.98, 0.06,
            "Bubble size reflects -log10(FDR)\n"
            "WGCNA context layer: n=44,\n"
            "scale-free fit R2=0.856",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.7, bbox=dict(facecolor="white", edgecolor="#bbbbbb", alpha=0.92))
    panel_label(ax, "B", "Module size, effect magnitude, and significance")

    # Panels C-D: radial hub-gene summaries.
    corr_lookup = dict(zip(modules["Module"], modules["Cor_FHL"]))
    radial_network(axes[1, 0], hub, 2, corr_lookup.get(2, np.nan))
    panel_label(axes[1, 0], "C", "Module 2 hub-gene network")
    radial_network(axes[1, 1], hub, 6, corr_lookup.get(6, np.nan))
    panel_label(axes[1, 1], "D", "Module 6 hub-gene network")

    legend_handles = [
        Patch(facecolor="#c95d73", edgecolor="#111111", label="Ferroptosis / iron / immune panel gene"),
        Patch(facecolor="#0b7f35", edgecolor="#111111", label="Immune or cell-death related hub"),
        Patch(facecolor="#80c8e8", edgecolor="#111111", label="Other hub gene"),
        Patch(facecolor="#3d2c8d", edgecolor="#111111", label="Top hub"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 0.025), fontsize=9.2)

    fig.subplots_adjust(left=0.08, right=0.94, top=0.90, bottom=0.15,
                        wspace=0.30, hspace=0.38)
    PNG_OUT.parent.mkdir(parents=True, exist_ok=True)
    TIF_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_OUT, dpi=300)
    fig.savefig(TIF_OUT, dpi=300)
    plt.close(fig)
    print(f"Wrote {PNG_OUT}")
    print(f"Wrote {TIF_OUT}")


if __name__ == "__main__":
    main()
