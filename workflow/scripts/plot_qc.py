#!/usr/bin/env python3
"""QC figures: annotation completeness and its effect on phenotype calls.

Inputs are aggregates over all genomes, so these render the same whether there
are 2 or 342,000 proteomes. Figures:

  qc_annotation_completeness  - distributions of Pfam families / annotation depth
                                / model-feature coverage per genome (with a
                                low-annotation flag line)
  qc_completeness_vs_calls    - 2D density of model-feature coverage vs number of
                                positive traits (are sparsely annotated genomes
                                driving the calls?)
  qc_calls_per_genome         - number of positive traits per genome, any-model
                                vs the reliable (phypat+PGL) set
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

import plot_style as S


def _write_captions(outdir, n, low):
    txt = f"""# QC figure captions  (n = {n:,} genomes)

**qc_annotation_completeness** — Per-genome annotation summary. (A) Number of
distinct Pfam families detected per genome; (B) total annotation depth (summed
per-family protein hits); (C) empirical cumulative distribution and (D) histogram
of model-feature coverage, the fraction of the phenotype models' ~Pfam feature
set present in each genome. Dashed lines mark the median and quartiles; the
shaded region in (A) flags genomes with fewer than {low} annotated families as
potentially incomplete.

**qc_completeness_vs_calls** — Joint density (log-scaled counts) of model-feature
coverage (x) versus the number of positively predicted traits per genome (y).
The white line is the median number of calls in each coverage bin. Used to check
whether calls are concentrated among well-annotated genomes.

**qc_calls_per_genome** — Distribution of the number of positively predicted
traits per genome, for the any-model set (combined code >= 1) and the reliable
phypat+PGL set (codes 2-3). Genomes with zero calls are annotated separately as
a QC flag.
"""
    with open(os.path.join(outdir, "FIGURES.qc.md"), "w") as fh:
        fh.write(txt)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation-stats", required=True)
    ap.add_argument("--combined", required=True, help="predictions_majority-vote_combined.tsv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--formats", default="png,svg")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--low-families", type=int, default=50)
    args = ap.parse_args(argv)
    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    S.apply_style(args.dpi)
    import matplotlib.pyplot as plt

    qc = pd.read_csv(args.annotation_stats, sep="\t", index_col=0)
    n = qc.shape[0]
    fam = qc["n_pfam_families"].to_numpy()
    depth = qc["total_pfam_hits"].to_numpy()
    cov = qc["frac_model_coverage"].to_numpy()

    c = S.load_calls(args.combined)
    n_any = (c.to_numpy() >= 1).sum(axis=1)
    n_rel = np.isin(c.to_numpy(), [2, 3]).sum(axis=1)
    # align calls-per-genome to qc order where possible (else just use as-is)
    cov_for_join = pd.Series(cov, index=qc.index)
    calls_series = pd.Series(n_any, index=c.index)
    joined = pd.concat([cov_for_join.rename("cov"), calls_series.rename("calls")],
                       axis=1, join="inner").dropna()

    # ---- Figure 1: completeness distributions (2x2) --------------------------
    fig = plt.figure(figsize=(S.WIDTH_FULL, 6), layout="constrained")
    gs = fig.add_gridspec(2, 2)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axC, axD = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    axA.hist(fam, bins=60, color="#4C72B0", edgecolor="white", linewidth=0.2)
    for q, ls in [(50, "-"), (25, "--"), (75, "--")]:
        axA.axvline(np.percentile(fam, q), color="0.25", ls=ls, lw=1)
    axA.axvspan(0, args.low_families, color="#C44E52", alpha=0.10)
    axA.set_xlabel("Pfam families per genome")
    axA.set_ylabel("genomes")
    S.panel_label(axA, "A")
    S.annotate_n(axA, n)

    axB.hist(np.log10(np.clip(depth, 1, None)), bins=60, color="#8172B3",
             edgecolor="white", linewidth=0.2)
    axB.set_xlabel("annotation depth (log$_{10}$ total Pfam hits)")
    axB.set_ylabel("genomes")
    S.panel_label(axB, "B")

    xv, yv = S.ecdf(cov)
    axC.plot(xv, yv, color="#55A868", lw=1.8)
    for q in (0.1, 0.5, 0.9):
        xq = np.quantile(cov, q)
        axC.plot([xq, xq], [0, q], color="0.5", ls=":", lw=1)
    axC.set_xlabel("model-feature coverage (fraction)")
    axC.set_ylabel("cumulative fraction of genomes")
    axC.set_ylim(0, 1)
    S.panel_label(axC, "C")

    axD.hist(cov, bins=60, color="#55A868", edgecolor="white", linewidth=0.2)
    axD.axvline(np.median(cov), color="0.25", lw=1)
    axD.set_xlabel("model-feature coverage (fraction)")
    axD.set_ylabel("genomes")
    S.panel_label(axD, "D")

    fig.suptitle("Annotation completeness across genomes", fontweight="bold")
    S.save_fig(fig, "qc_annotation_completeness", args.outdir, formats, args.dpi)

    # ---- Figure 2: coverage vs number of calls (hexbin density) --------------
    fig, ax = plt.subplots(figsize=(S.WIDTH_1_5COL, 4.2))
    hb = ax.hexbin(joined["cov"].to_numpy(), joined["calls"].to_numpy(),
                   gridsize=55, cmap="viridis", bins="log", mincnt=1)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label("genomes (log scale)")
    # median calls per coverage bin
    nb = 25
    edges = np.linspace(joined["cov"].min(), joined["cov"].max(), nb + 1)
    idx = np.clip(np.digitize(joined["cov"].to_numpy(), edges) - 1, 0, nb - 1)
    med = [np.median(joined["calls"].to_numpy()[idx == b]) if np.any(idx == b) else np.nan
           for b in range(nb)]
    ax.plot(0.5 * (edges[:-1] + edges[1:]), med, color="white", lw=2)
    ax.set_xlabel("model-feature coverage (fraction)")
    ax.set_ylabel("positive traits per genome")
    ax.set_title("Annotation coverage vs. phenotype calls")
    S.annotate_n(ax, joined.shape[0], loc="upper left")
    S.save_fig(fig, "qc_completeness_vs_calls", args.outdir, formats, args.dpi)

    # ---- Figure 3: calls per genome ------------------------------------------
    fig, ax = plt.subplots(figsize=(S.WIDTH_1_5COL, 4.0))
    bins = np.arange(0, 68) - 0.5
    ax.hist(n_any, bins=bins, color=S.MODEL_COLORS["phypat"], alpha=0.55,
            label="any model (code ≥ 1)")
    ax.hist(n_rel, bins=bins, color=S.MODEL_COLORS["reliable"], alpha=0.7,
            label="reliable (phypat+PGL)")
    zero = int((n_any == 0).sum())
    ax.set_xlabel("number of positive traits per genome")
    ax.set_ylabel("genomes")
    ax.set_title("Positive phenotype calls per genome")
    ax.legend(frameon=False)
    ax.text(0.98, 0.80, f"{zero:,} genomes with 0 calls",
            transform=ax.transAxes, ha="right", fontsize=8, color="#C44E52")
    S.annotate_n(ax, n)
    S.save_fig(fig, "qc_calls_per_genome", args.outdir, formats, args.dpi)

    _write_captions(args.outdir, n, args.low_families)
    print(f"[plot_qc] wrote QC figures for {n:,} genomes -> {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
