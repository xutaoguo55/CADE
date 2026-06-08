#!/usr/bin/env python3
"""Regenerate Figure 5 with separated titles and non-overlapping labels."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PNG_OUT = ROOT / "manuscript" / "docx_embedded_figures_png" / "Figure5_Pathway_Analysis.png"
TIF_OUT = ROOT / "figures" / "Figure5_Pathway_Analysis.tif"


def panel_label(ax, label: str, title: str) -> None:
    ax.text(-0.14, 1.06, label, transform=ax.transAxes, fontsize=13,
            weight="bold", va="bottom", ha="left")
    ax.set_title(title, fontsize=12, weight="bold", loc="left", pad=10)


def signif(fdr: float) -> str:
    if fdr < 0.001:
        return "***"
    if fdr < 0.01:
        return "**"
    if fdr < 0.05:
        return "*"
    return "ns"


def pretty(name: str) -> str:
    return name.replace("_", " ")


def main() -> None:
    ferro = pd.read_csv(ROOT / "supplementary" / "Table_S06_Ferroptosis_GSVA_Scores.csv")
    death = pd.read_csv(ROOT / "supplementary" / "Table_S07_CellDeath_Pathways.csv")

    order_ferro = ["Iron_Homeostasis", "Ferroptosis_Core", "Ferroptosis_Defense",
                   "GSH_Metabolism", "Ferroptosis_Drivers"]
    ferro = ferro.set_index("GeneSet").loc[order_ferro].reset_index()
    order_death = ["Autophagy", "Ferroptosis", "Necroptosis", "Apoptosis_Intrinsic",
                   "Programmed_Cell_Death", "Apoptosis_Extrinsic", "Pyroptosis"]
    death = death.set_index("Pathway").loc[order_death].reset_index()

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
    })
    rose = "#c95d73"
    green = "#0b7f35"
    gold = "#dccb70"
    gray = "#bfbfbf"
    blue = "#7fc4e4"
    purple = "#372b8c"

    fig, axes = plt.subplots(1, 3, figsize=(18, 6.3), dpi=300,
                             gridspec_kw={"width_ratios": [1.28, 1.28, 1.55]})
    fig.suptitle("Figure 5 | Pathway-level patterns from per-sample GSVA",
                 fontsize=16, weight="bold", y=0.98)

    # Panel A
    ax = axes[0]
    y = np.arange(len(ferro))
    ax.barh(y, ferro["Delta"], color=rose, edgecolor="#333333", height=0.72)
    ax.axvline(0, color="#333333", lw=0.8)
    for yi, delta, fdr in zip(y, ferro["Delta"], ferro["FDR"]):
        ax.text(delta + 0.018, yi, f"+{delta:.2f}  {signif(fdr)}",
                va="center", ha="left", fontsize=9.5)
    ax.set_yticks(y)
    ax.set_yticklabels([pretty(v) for v in ferro["GeneSet"]])
    ax.invert_yaxis()
    ax.set_xlim(-0.2, 0.72)
    ax.set_xlabel("Delta GSVA score (FHL - HC)")
    panel_label(ax, "A", "Ferroptosis and iron-related pathways")

    # Panel B
    ax = axes[1]
    y = np.arange(len(death))
    colors = [green if n in ["Autophagy", "Ferroptosis"] else
              gold if n in ["Necroptosis", "Apoptosis_Intrinsic"] else
              gray if v >= 0 else blue
              for n, v in zip(death["Pathway"], death["Delta"])]
    ax.barh(y, death["Delta"], color=colors, edgecolor="#333333", height=0.72)
    ax.axvline(0, color="#333333", lw=0.8)
    for yi, delta, fdr in zip(y, death["Delta"], death["FDR"]):
        ha = "left" if delta >= 0 else "right"
        x = delta + 0.018 if delta >= 0 else delta - 0.018
        ax.text(x, yi, f"{delta:+.2f}  {signif(fdr)}",
                va="center", ha=ha, fontsize=9.5)
    ax.set_yticks(y)
    ax.set_yticklabels([pretty(v) for v in death["Pathway"]])
    ax.invert_yaxis()
    ax.set_xlim(-0.22, 0.72)
    ax.set_xlabel("Delta GSVA score (FHL - HC)")
    panel_label(ax, "B", "Multi-cell-death pathway comparison")

    # Panel C
    ax = axes[2]
    scatter = pd.concat([
        ferro.rename(columns={"GeneSet": "Name"}).assign(Group="Ferroptosis / iron sets"),
        death.rename(columns={"Pathway": "Name"}).assign(Group="Cell-death sets"),
    ], ignore_index=True)
    scatter["neglog10FDR"] = -np.log10(scatter["FDR"].clip(lower=1e-300))
    for group, color in [("Ferroptosis / iron sets", rose), ("Cell-death sets", purple)]:
        sub = scatter[scatter["Group"] == group]
        ax.scatter(sub["Delta"], sub["neglog10FDR"], s=95, color=color,
                   edgecolor="#222222", alpha=0.9, label=group, zorder=3)
    ax.axhline(-np.log10(0.05), color="#999999", lw=1.0, ls="--")
    ax.text(0.69, -np.log10(0.05) + 0.06, "FDR=0.05", color="#777777",
            fontsize=9, ha="right")
    offsets = {
        "Ferroptosis_Core": (0.018, 0.20), "Ferroptosis_Defense": (0.018, 0.28),
        "Ferroptosis_Drivers": (0.018, -0.18), "Iron_Homeostasis": (0.018, -0.04),
        "GSH_Metabolism": (0.018, 0.08), "Ferroptosis": (0.018, 0.12),
        "Autophagy": (0.015, -0.05), "Apoptosis_Intrinsic": (0.018, -0.12),
        "Programmed_Cell_Death": (0.018, 0.03), "Necroptosis": (0.018, 0.12),
        "Apoptosis_Extrinsic": (0.018, 0.02), "Pyroptosis": (0.018, -0.02),
    }
    labeled_low_points = {"Pyroptosis", "Apoptosis_Extrinsic"}
    for _, row in scatter.iterrows():
        if row["FDR"] >= 0.05 and row["Name"] not in labeled_low_points:
            continue
        dx, dy = offsets.get(row["Name"], (0.015, 0.0))
        label = "Cell-death ferroptosis" if row["Name"] == "Ferroptosis" else pretty(row["Name"])
        ax.text(row["Delta"] + dx, row["neglog10FDR"] + dy, label,
                fontsize=8.4, ha="left", va="center")
    ax.set_xlim(-0.24, 0.92)
    ax.set_ylim(0, max(scatter["neglog10FDR"]) + 0.9)
    ax.set_xlabel("Delta GSVA score (FHL - HC)")
    ax.set_ylabel("-log10(FDR)")
    ax.legend(frameon=True, loc="lower right", fontsize=8.8)
    panel_label(ax, "C", "Pathway-level Delta vs significance")

    fig.subplots_adjust(left=0.12, right=0.985, top=0.86, bottom=0.14,
                        wspace=0.48)
    PNG_OUT.parent.mkdir(parents=True, exist_ok=True)
    TIF_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_OUT, dpi=300)
    fig.savefig(TIF_OUT, dpi=300)
    plt.close(fig)
    print(f"Wrote {PNG_OUT}")
    print(f"Wrote {TIF_OUT}")


if __name__ == "__main__":
    main()
