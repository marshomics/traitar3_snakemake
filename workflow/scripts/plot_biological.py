#!/usr/bin/env python3
"""Biological interpretation figures (aggregated over all genomes).

  bio_prevalence_by_category - reliable (phypat+PGL) trait prevalence grouped and
                               coloured by biological category
  bio_trait_cooccurrence     - 67x67 clustered heatmap of trait co-occurrence
                               (phi / Pearson correlation of reliable calls),
                               revealing structure such as mutually exclusive
                               Gram stain or co-utilised sugars
  bio_key_traits_composition - dataset-level prevalence of curated key traits
                               (oxygen relationship, Gram stain, morphology, ...)
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

import plot_style as S

KEY_TRAITS = [
    ("Oxygen", ["Aerobe", "Anaerobe", "Facultative", "Capnophilic"]),
    ("Gram stain", ["Gram positive", "Gram negative"]),
    ("Morphology", ["Coccus", "Bacillus or coccobacillus", "Motile", "Spore formation"]),
    ("Key enzymes", ["Catalase", "Oxidase"]),
]


def _reliable(c: pd.DataFrame) -> np.ndarray:
    return np.isin(c.to_numpy(), [2, 3]).astype(np.float32)


def _write_captions(outdir, n):
    txt = f"""# Biological figure captions  (n = {n:,} genomes)

**bio_prevalence_by_category** — Fraction of genomes positive for each trait
(reliable phypat+PGL calls), grouped by biological category and sorted within
each category. Provides a biological overview of which capabilities are common
or rare across the dataset.

**bio_trait_cooccurrence** — Hierarchically clustered heatmap of pairwise trait
co-occurrence across genomes (Pearson/phi correlation of reliable binary calls).
Positive values indicate traits that tend to co-occur; negative values indicate
mutually exclusive traits (e.g. Gram positive vs Gram negative, aerobe vs
anaerobe). Side colour bars denote trait category.

**bio_key_traits_composition** — Prevalence of curated key traits describing
oxygen relationship, Gram stain, cell morphology and marker enzymes, giving a
dataset-level biological composition.
"""
    with open(os.path.join(outdir, "FIGURES.biological.md"), "w") as fh:
        fh.write(txt)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--combined", required=True)
    ap.add_argument("--model-tar", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--formats", default="png,svg")
    ap.add_argument("--dpi", type=int, default=600)
    args = ap.parse_args(argv)
    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    S.apply_style(args.dpi)
    import matplotlib.pyplot as plt
    import seaborn as sns

    c = S.load_calls(args.combined)
    n = c.shape[0]
    traits = list(c.columns)
    cats = S.load_trait_categories(args.model_tar).reindex(traits).fillna("Other")
    rel = _reliable(c)
    prev = pd.Series(rel.mean(axis=0), index=traits)
    catcol = S.category_palette(cats)

    # ---- Figure 1: prevalence grouped by category ----------------------------
    ordered, colors, ylabels, sep = [], [], [], []
    pos = 0
    for cat in sorted(cats.unique()):
        members = prev[cats == cat].sort_values()
        for t in members.index:
            ordered.append(prev[t]); colors.append(catcol[cat]); ylabels.append(t); pos += 1
        sep.append(pos - 0.5)
    fig, ax = plt.subplots(figsize=(S.WIDTH_1_5COL, 11))
    yy = np.arange(len(ordered))
    ax.barh(yy, ordered, color=colors, edgecolor="white", linewidth=0.2)
    ax.set_yticks(yy)
    ax.set_yticklabels(ylabels, fontsize=7)
    for s in sep[:-1]:
        ax.axhline(s, color="0.85", lw=0.6)
    ax.set_xlim(0, 1)
    ax.set_ylim(-1, len(ordered))
    ax.set_xlabel("fraction of genomes positive (reliable calls)")
    ax.set_title("Trait prevalence by category")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=catcol[k], label=k) for k in sorted(catcol)],
              loc="upper center", bbox_to_anchor=(0.5, -0.055), ncol=3, frameon=False,
              fontsize=7, title="trait category")
    S.annotate_n(ax, n, loc="upper left")
    S.save_fig(fig, "bio_prevalence_by_category", args.outdir, formats, args.dpi)

    # ---- Figure 2: trait co-occurrence clustered heatmap ---------------------
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = np.corrcoef(rel.T)
    corr = np.nan_to_num(corr, nan=0.0)
    np.fill_diagonal(corr, 1.0)
    corr_df = pd.DataFrame(corr, index=traits, columns=traits)
    col_colors = pd.Series({t: catcol[cats[t]] for t in traits}, name="category")
    g = sns.clustermap(corr_df, cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                       row_colors=col_colors, col_colors=col_colors,
                       xticklabels=True, yticklabels=True,
                       figsize=(11, 11), cbar_kws={"label": "co-occurrence (r)"},
                       dendrogram_ratio=0.12, linewidths=0)
    g.ax_heatmap.tick_params(labelsize=6)
    g.fig.suptitle(f"Trait co-occurrence across genomes (n = {n:,})", y=1.01, fontweight="bold")
    os.makedirs(args.outdir, exist_ok=True)
    for fmt in formats:
        g.savefig(os.path.join(args.outdir, f"bio_trait_cooccurrence.{fmt}"),
                  dpi=args.dpi, bbox_inches="tight")
    plt.close(g.fig)

    # ---- Figure 3: key traits composition ------------------------------------
    fig, ax = plt.subplots(figsize=(S.WIDTH_1_5COL, 4.6))
    labels, vals, colors, group_sep, group_mid, group_name = [], [], [], [], [], []
    pos = 0
    grp_palette = S.category_palette([g for g, _ in KEY_TRAITS])
    for gname, members in KEY_TRAITS:
        present = [t for t in members if t in prev.index]
        start = pos
        for t in present:
            labels.append(t); vals.append(prev[t]); colors.append(grp_palette[gname]); pos += 1
        if present:
            group_mid.append((start + pos - 1) / 2)
            group_name.append(gname)
            group_sep.append(pos - 0.5)
    yy = np.arange(len(labels))
    ax.barh(yy, vals, color=colors, edgecolor="white", linewidth=0.3)
    ax.set_yticks(yy); ax.set_yticklabels(labels, fontsize=8)
    for s in group_sep[:-1]:
        ax.axhline(s, color="0.8", lw=0.7)
    ax.set_xlim(0, 1)
    ax.set_ylim(-1, len(labels))
    ax.set_xlabel("fraction of genomes positive (reliable calls)")
    ax.set_title("Key trait composition of the dataset")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=grp_palette[g], label=g) for g, _ in KEY_TRAITS],
              loc="lower right", frameon=False, fontsize=8)
    S.annotate_n(ax, n, loc="upper right")
    S.save_fig(fig, "bio_key_traits_composition", args.outdir, formats, args.dpi)

    _write_captions(args.outdir, n)
    print(f"[plot_biological] wrote biological figures for {n:,} genomes -> {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
