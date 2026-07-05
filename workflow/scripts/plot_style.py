#!/usr/bin/env python3
"""Shared publication styling and helpers for the traitar3-smk figures.

All figures are written as a matched PNG (raster, high DPI) + SVG pair. SVG text
is kept as editable ``<text>`` elements (``svg.fonttype = 'none'``) so the files
open cleanly in Illustrator/Inkscape for manuscript tweaking.

Design rule for 342k genomes: never draw one mark per genome. Every figure here
consumes an aggregate (per-trait column reductions, histograms/ECDFs of
per-genome statistics, 2D density, or a trait x trait summary), so cost and
readability are independent of the number of genomes.
"""
from __future__ import annotations

import os
import tarfile

import matplotlib as mpl

mpl.use("Agg")  # headless / cluster
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# journal column widths (inches)
WIDTH_1COL = 3.5
WIDTH_1_5COL = 5.5
WIDTH_FULL = 7.2

# combined-call code colours (0/1/2/3) and per-model colours
CODE_COLORS = {
    0: "#d9d9d9",  # negative in both
    1: "#4C72B0",  # phypat only
    2: "#DD8452",  # phypat+PGL only
    3: "#55A868",  # both
}
CODE_LABELS = {
    0: "negative (both)",
    1: "phypat only",
    2: "phypat+PGL only",
    3: "both models",
}
MODEL_COLORS = {"phypat": "#4C72B0", "phypat+PGL": "#DD8452", "reliable": "#55A868"}


def apply_style(dpi: int = 600) -> None:
    """Classic-academic matplotlib defaults; keeps SVG text editable."""
    try:
        import seaborn as sns
        sns.set_context("paper", font_scale=1.15)
        sns.set_style("ticks")
    except Exception:
        pass
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "axes.linewidth": 1.0,
        "figure.dpi": 150,
        "savefig.dpi": dpi,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "svg.fonttype": "none",   # <-- editable text in SVG
        "pdf.fonttype": 42,
        "figure.autolayout": False,
    })


def save_fig(fig, name: str, outdir: str, formats=("png", "svg"), dpi: int = 600) -> list:
    """Write fig to <outdir>/<name>.<fmt> for each requested format."""
    os.makedirs(outdir, exist_ok=True)
    written = []
    for fmt in formats:
        path = os.path.join(outdir, f"{name}.{fmt}")
        fig.savefig(path, format=fmt, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
        written.append(path)
    plt.close(fig)
    return written


def panel_label(ax, label: str, x: float = -0.12, y: float = 1.06) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")


def annotate_n(ax, n: int, loc: str = "upper right") -> None:
    xy = {"upper right": (0.98, 0.98, "right", "top"),
          "upper left": (0.02, 0.98, "left", "top")}[loc]
    ax.text(xy[0], xy[1], f"n = {n:,} genomes", transform=ax.transAxes,
            ha=xy[2], va=xy[3], fontsize=8, color="0.35")


# ---------------------------------------------------------------------------
# scale-safe loaders
# ---------------------------------------------------------------------------
def load_calls(path: str, dtype="int16") -> pd.DataFrame:
    """Load a genome x trait matrix compactly (calls are small integers)."""
    head = pd.read_csv(path, sep="\t", index_col=0, nrows=0)
    cols = {c: dtype for c in head.columns}
    df = pd.read_csv(path, sep="\t", index_col=0, dtype=cols)
    df.index = df.index.astype(str)
    return df


def load_trait_categories(model_tar: str) -> pd.Series:
    """trait_name -> category, from a phenotype model archive's pt2acc.txt."""
    with tarfile.open(model_tar, "r:gz") as tf:
        pt2acc = pd.read_csv(tf.extractfile("pt2acc.txt"), sep="\t", index_col=0)
    # columns: accession (trait name), category
    s = pt2acc.set_index("accession")["category"]
    s.index = s.index.astype(str)
    return s


def category_palette(categories) -> dict:
    """Stable colour per trait category (sorted for determinism)."""
    cats = sorted(set(categories))
    try:
        import seaborn as sns
        colors = sns.color_palette("tab20", max(len(cats), 3))
    except Exception:
        cmap = plt.get_cmap("tab20")
        colors = [cmap(i / max(len(cats) - 1, 1)) for i in range(len(cats))]
    return {c: colors[i] for i, c in enumerate(cats)}


def ecdf(values: np.ndarray):
    """Return (sorted_values, cumulative_fraction) for an ECDF step plot."""
    v = np.sort(np.asarray(values, dtype=float))
    if v.size == 0:
        return v, v
    y = np.arange(1, v.size + 1) / v.size
    return v, y
