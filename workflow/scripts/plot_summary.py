#!/usr/bin/env python3
"""Prediction summary figures (aggregated over all genomes).

  summary_trait_prevalence  - per-trait positive prevalence, phypat vs phypat+PGL
                              (dumbbell), trait labels coloured by category
  summary_model_agreement   - per-trait stacked fractions of the combined 0/1/2/3
                              code (negative / phypat-only / PGL-only / both)
  summary_confidence        - committee vote support: overall per-cell vote
                              distribution and support among positive calls
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

import plot_style as S


def _prevalence(c: pd.DataFrame):
    v = c.to_numpy()
    phypat = np.isin(v, [1, 3]).mean(axis=0)
    pgl = np.isin(v, [2, 3]).mean(axis=0)
    return (pd.Series(phypat, index=c.columns), pd.Series(pgl, index=c.columns))


def _write_captions(outdir, n, thr):
    txt = f"""# Summary figure captions  (n = {n:,} genomes)

**summary_trait_prevalence** — Fraction of genomes called positive for each of
the 67 traits by the phypat and phypat+PGL classifiers (dumbbell endpoints),
sorted by phypat+PGL prevalence. Trait labels are coloured by biological
category.

**summary_model_agreement** — For each trait, the fraction of genomes in each
combined call class: negative in both, phypat only, phypat+PGL only, or both
models. Sorted by the fraction supported by both models. Quantifies concordance
between the two classifiers.

**summary_confidence** — Committee support. (A) Distribution of per-cell voter
agreement (0-5 of the 5 SVM voters) for each model across all genome x trait
cells. (B) Among positive calls (>= {thr} of 5 voters), the fraction at each
support level, showing that positive calls are predominantly strongly supported
rather than marginal.
"""
    with open(os.path.join(outdir, "FIGURES.summary.md"), "w") as fh:
        fh.write(txt)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--combined", required=True)
    ap.add_argument("--phypat-single", required=True)
    ap.add_argument("--phypatpgl-single", required=True)
    ap.add_argument("--model-tar", required=True, help="phypat.tar.gz (for categories)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--formats", default="png,svg")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--majority-threshold", type=float, default=3)
    args = ap.parse_args(argv)
    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    S.apply_style(args.dpi)
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    c = S.load_calls(args.combined)
    n = c.shape[0]
    cats = S.load_trait_categories(args.model_tar)
    prev_phypat, prev_pgl = _prevalence(c)
    order = prev_pgl.sort_values().index.tolist()
    catcol = S.category_palette(cats.reindex(order).fillna("Other"))

    # ---- Figure 1: prevalence dumbbell ---------------------------------------
    fig, ax = plt.subplots(figsize=(S.WIDTH_1_5COL, 11))
    y = np.arange(len(order))
    for yi, t in zip(y, order):
        ax.plot([prev_phypat[t], prev_pgl[t]], [yi, yi], color="0.7", lw=1.2, zorder=1)
    ax.scatter(prev_phypat[order], y, s=22, color=S.MODEL_COLORS["phypat"], label="phypat", zorder=2)
    ax.scatter(prev_pgl[order], y, s=22, color=S.MODEL_COLORS["phypat+PGL"], label="phypat+PGL", zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    for tick, t in zip(ax.get_yticklabels(), order):
        tick.set_color(catcol.get(cats.get(t, "Other"), "0.2"))
        tick.set_fontsize(7)
    ax.set_xlabel("fraction of genomes called positive")
    ax.set_ylim(-1, len(order))
    ax.set_title("Trait prevalence by classifier")
    model_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=S.MODEL_COLORS["phypat"], label="phypat", markersize=7),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=S.MODEL_COLORS["phypat+PGL"], label="phypat+PGL", markersize=7),
    ]
    leg_model = ax.legend(handles=model_handles, loc="lower right", frameon=False)
    ax.add_artist(leg_model)
    cat_handles = [Patch(color=catcol[k], label=k) for k in sorted(catcol)]
    ax.legend(handles=cat_handles, loc="upper center", bbox_to_anchor=(0.5, -0.055),
              ncol=3, frameon=False, fontsize=7, title="trait category")
    S.annotate_n(ax, n, loc="upper left")
    S.save_fig(fig, "summary_trait_prevalence", args.outdir, formats, args.dpi)

    # ---- Figure 2: model agreement stacked -----------------------------------
    v = c.to_numpy()
    frac = {code: np.isin(v, [code]).mean(axis=0) for code in (0, 1, 2, 3)}
    fdf = pd.DataFrame(frac, index=c.columns)
    order2 = fdf[3].sort_values().index.tolist()
    fig, ax = plt.subplots(figsize=(S.WIDTH_1_5COL, 11))
    y = np.arange(len(order2))
    left = np.zeros(len(order2))
    for code in (0, 1, 2, 3):
        vals = fdf.loc[order2, code].to_numpy()
        ax.barh(y, vals, left=left, color=S.CODE_COLORS[code], label=S.CODE_LABELS[code],
                edgecolor="white", linewidth=0.2)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(order2, fontsize=7)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of genomes")
    ax.set_title("Classifier agreement per trait")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.055), ncol=4,
              frameon=False, fontsize=8)
    S.annotate_n(ax, n, loc="upper left")
    S.save_fig(fig, "summary_model_agreement", args.outdir, formats, args.dpi)

    # ---- Figure 3: confidence -------------------------------------------------
    sp = S.load_calls(args.phypat_single, dtype="int8").to_numpy().ravel()
    sg = S.load_calls(args.phypatpgl_single, dtype="int8").to_numpy().ravel()
    fig = plt.figure(figsize=(S.WIDTH_FULL, 3.6), layout="constrained")
    gs = fig.add_gridspec(1, 2)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    levels = np.arange(0, 6)
    w = 0.4
    for off, data, name in [(-w / 2, sp, "phypat"), (w / 2, sg, "phypat+PGL")]:
        frac = np.array([(data == L).mean() for L in levels])
        axA.bar(levels + off, frac, width=w, color=S.MODEL_COLORS[name], label=name)
    axA.set_xlabel("voters in agreement (of 5)")
    axA.set_ylabel("fraction of all cells")
    axA.set_title("Overall vote support")
    axA.legend(frameon=False)
    S.panel_label(axA, "A")
    thr = args.majority_threshold
    pos_levels = [L for L in levels if L >= thr]
    for off, data, name in [(-w / 2, sp, "phypat"), (w / 2, sg, "phypat+PGL")]:
        pos = data[data >= thr]
        frac = np.array([(pos == L).mean() if pos.size else 0 for L in pos_levels])
        axB.bar(np.array(pos_levels) + off, frac, width=w, color=S.MODEL_COLORS[name], label=name)
    axB.set_xlabel(f"voters in agreement (positive calls, ≥ {thr:g})")
    axB.set_ylabel("fraction of positive calls")
    axB.set_title("Support among positive calls")
    axB.set_xticks(pos_levels)
    S.panel_label(axB, "B")
    fig.suptitle("Committee vote support", fontweight="bold")
    S.save_fig(fig, "summary_confidence", args.outdir, formats, args.dpi)

    _write_captions(args.outdir, n, args.majority_threshold)
    print(f"[plot_summary] wrote summary figures for {n:,} genomes -> {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
