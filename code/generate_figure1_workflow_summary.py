#!/usr/bin/env python3
"""Regenerate Figure 1 with non-overlapping labels for the manuscript DOCX."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
PNG_OUT = ROOT / "manuscript" / "docx_embedded_figures_png" / "Figure1_Workflow_Summary.png"
TIF_OUT = ROOT / "figures" / "Figure1_Workflow_Summary.tif"

PURPLE = "#3b2a8f"
BLUE = "#6ec5e9"
ORANGE = "#f6a045"
MAGENTA = "#9a145d"
GRAY = "#a0a0a0"
BOX = "#f7f9fc"
LINE = "#333333"


def add_panel_label(ax, label, title):
    ax.text(0.0, 1.03, label, transform=ax.transAxes, fontsize=15, weight="bold",
            va="bottom", ha="left")
    ax.text(0.065, 1.03, title, transform=ax.transAxes, fontsize=15, weight="bold",
            va="bottom", ha="left")


def rounded_box(ax, xy, width, height, text="", fc="white", ec=PURPLE,
                lw=1.6, fontsize=10, weight="normal", color="black"):
    patch = FancyBboxPatch(
        xy, width, height,
        boxstyle="round,pad=0.018,rounding_size=0.018",
        linewidth=lw, edgecolor=ec, facecolor=fc,
        transform=ax.transAxes, clip_on=False,
    )
    ax.add_patch(patch)
    if text:
        ax.text(
            xy[0] + width / 2, xy[1] + height / 2, text,
            transform=ax.transAxes, ha="center", va="center",
            fontsize=fontsize, weight=weight, color=color, linespacing=1.14,
        )
    return patch


def panel_a(ax):
    ax.set_axis_off()
    add_panel_label(ax, "A", "The cell-composition confounding problem")

    categories = [
        ("CD8/CD4/B/NK", BLUE),
        ("Myeloid", ORANGE),
        ("Erythroid", MAGENTA),
        ("Other", GRAY),
    ]
    healthy = [0.52, 0.18, 0.12, 0.18]
    fhl = [0.12, 0.39, 0.31, 0.18]
    x_positions = [0.12, 0.30]
    bar_w = 0.11
    bottom_y = 0.34
    bar_h = 0.44

    for x, vals, label in zip(x_positions, [healthy, fhl], ["Healthy control", "FHL"]):
        y = bottom_y
        for value, (_, color) in zip(vals, categories):
            h = bar_h * value
            ax.add_patch(Rectangle((x, y), bar_w, h, transform=ax.transAxes,
                                   facecolor=color, edgecolor=LINE, lw=0.8))
            y += h
        ax.add_patch(Rectangle((x, bottom_y), bar_w, bar_h, transform=ax.transAxes,
                               fill=False, edgecolor=LINE, lw=1.0))
        ax.text(x + bar_w / 2, bottom_y - 0.055, label, transform=ax.transAxes,
                ha="center", va="top", fontsize=11, weight="bold")

    legend_x = 0.45
    for idx, (name, color) in enumerate(categories):
        y = 0.68 - idx * 0.07
        ax.add_patch(Rectangle((legend_x, y), 0.025, 0.04, transform=ax.transAxes,
                               facecolor=color, edgecolor=LINE, lw=0.7))
        ax.text(legend_x + 0.035, y + 0.02, name, transform=ax.transAxes,
                va="center", ha="left", fontsize=10)

    ax.add_patch(FancyArrowPatch((0.38, 0.55), (0.60, 0.64), transform=ax.transAxes,
                                 arrowstyle="-|>", mutation_scale=12, lw=1.4,
                                 color="#555555"))
    rounded_box(
        ax, (0.60, 0.53), 0.34, 0.18,
        text=r"Bulk RNA:  $Y_g = \sum_{k=1}^{K} p_k e_{gk} + \epsilon$"
             "\n(weighted mixture of cell types)",
        fc="#f2f2fb", ec=PURPLE, fontsize=11
    )
    rounded_box(
        ax, (0.06, 0.09), 0.88, 0.14,
        text="Standard DE on bulk expression can conflate two signals:\n"
             "cell-intrinsic regulation and disease-associated composition shift",
        fc="#fff3ed", ec="#e68a9a", fontsize=11
    )


def panel_b(ax):
    ax.set_axis_off()
    add_panel_label(ax, "B", "CADE three-stage framework (reference-free, marker-based)")

    box_y = 0.30
    box_w = 0.25
    box_h = 0.42
    xs = [0.05, 0.375, 0.70]
    texts = [
        "Stage 1\nMarker-derived weights\n\nCurated marker sets\n7-9 markers per lineage\nsoftmax or QP refinement\nrelative weights $w_{ki}$",
        "Stage 2\nComposition-adjusted DE\n\nlimma model:\nexpression ~ group + weights\nempirical-Bayes moderation\nadjusted logFC",
        "Stage 3\nCCI sensitivity index\n\ncoefficient-change ratio\nclipped to [0,1]\nbootstrap and permutation\nrank-stability summaries",
    ]
    for x, text in zip(xs, texts):
        rounded_box(ax, (x, box_y), box_w, box_h, text=text, fc="white",
                    ec=PURPLE, fontsize=10, weight="normal", color=PURPLE)

    for x1, x2 in [(xs[0] + box_w, xs[1]), (xs[1] + box_w, xs[2])]:
        ax.add_patch(FancyArrowPatch((x1 + 0.015, box_y + box_h / 2),
                                     (x2 - 0.015, box_y + box_h / 2),
                                     transform=ax.transAxes, arrowstyle="-|>",
                                     mutation_scale=13, lw=1.6, color=PURPLE))

    rounded_box(
        ax, (0.05, 0.07), 0.90, 0.13,
        text="Input: bulk expression matrix  +  group labels  +  curated marker gene sets\n"
             "No matched single-cell reference required",
        fc="#eef4f8", ec="#555555", fontsize=10.5
    )


def panel_c(ax):
    ax.set_axis_off()
    add_panel_label(ax, "C", "CCI interpretation")

    grad = np.linspace(0, 1, 256).reshape(1, -1)
    ax.imshow(grad, extent=(0.15, 0.90, 0.62, 0.70), transform=ax.transAxes,
              cmap="RdYlGn_r", aspect="auto")
    ax.add_patch(Rectangle((0.15, 0.62), 0.75, 0.08, transform=ax.transAxes,
                           fill=False, ec=LINE, lw=0.9))
    for xpos, label in [(0.15, "0.0"), (0.34, "0.2"), (0.53, "0.5"),
                        (0.72, "0.8"), (0.90, "1.0")]:
        ax.text(xpos, 0.58, label, transform=ax.transAxes, ha="center",
                va="top", fontsize=9)
    ax.text(0.525, 0.49, "CCI", transform=ax.transAxes, ha="center",
            va="center", fontsize=12, weight="bold")

    tier_y = 0.26
    tiers = [
        (0.22, "Lowest\nCCI < 0.2\ncoefficient stable\nafter adjustment", "#008a3b"),
        (0.52, "Low-moderate\n0.2-0.5\nmixed sensitivity", "#a26a00"),
        (0.81, "High\nCCI > 0.5\nstrongly\nmodel-sensitive", "#d7191c"),
    ]
    for x, text, color in tiers:
        ax.text(x, tier_y, text, transform=ax.transAxes, ha="center",
                va="center", fontsize=10, color=color, linespacing=1.08)

    ax.text(0.525, 0.08, "CCI reports coefficient sensitivity,\nnot mediation fraction.",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=10, style="italic", color="#555555")


def panel_d(ax):
    ax.set_axis_off()
    add_panel_label(ax, "D", "CADE outputs")

    items = [
        ("Per-gene table", "logFC_unadj, logFC_adj,\nDelta logFC, CCI, direction"),
        ("Bootstrap intervals", "marker resampling\nand rank stability"),
        ("Permutation null", "two-sided empirical\nP-values and FDR"),
        ("Cell-type weights", "marker-derived\nrelative composition"),
    ]
    y0 = 0.69
    for idx, (title, detail) in enumerate(items):
        y = y0 - idx * 0.20
        rounded_box(ax, (0.08, y), 0.84, 0.14, fc="#fbf8ef", ec="#e5c85a", lw=1.1)
        ax.text(0.13, y + 0.09, title, transform=ax.transAxes, ha="left",
                va="center", fontsize=11, weight="bold")
        ax.text(0.88, y + 0.05, detail, transform=ax.transAxes, ha="right",
                va="center", fontsize=9.5, color="#555555", linespacing=1.05)


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "axes.unicode_minus": False,
    })
    fig = plt.figure(figsize=(16, 10), dpi=300, facecolor="white")
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1.9, 1.0],
        height_ratios=[1.0, 1.0],
        left=0.035, right=0.985, bottom=0.05, top=0.90,
        wspace=0.12, hspace=0.23,
    )
    fig.suptitle("Figure 1 | CADE: concept and workflow",
                 fontsize=17, weight="bold", y=0.975)

    panel_a(fig.add_subplot(gs[0, 0]))
    panel_c(fig.add_subplot(gs[0, 1]))
    panel_b(fig.add_subplot(gs[1, 0]))
    panel_d(fig.add_subplot(gs[1, 1]))

    PNG_OUT.parent.mkdir(parents=True, exist_ok=True)
    TIF_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_OUT, dpi=300)
    fig.savefig(TIF_OUT, dpi=300)
    plt.close(fig)
    print(f"Wrote {PNG_OUT}")
    print(f"Wrote {TIF_OUT}")


if __name__ == "__main__":
    main()
