#!/usr/bin/env python3
"""Write one majority-vote output file per genome for a single batch.

Consumes the per-batch, per-model majority-vote tables that predict_batch
already produced (phypat and phypat+PGL, genomes x 67 traits, 0/1), and writes
one file per genome into the shared output directory:

    <outdir>/<sample>.majority_vote.tsv
      phenotype   category   phypat   phypat+PGL   combined

`combined` uses traitar's 0/1/2/3 encoding (0 negative in both, 1 phypat only,
2 phypat+PGL only, 3 both) -- identical to the corresponding row of
predictions_majority-vote_combined.tsv. A per-batch sentinel is written last so
the workflow can track completion without declaring 350k individual outputs.

This reads existing prediction tables only; it does not recompute anything and
does not touch the combined matrix that the plots depend on.
"""
from __future__ import annotations

import argparse
import csv
import os
import tarfile

import pandas as pd


def load_categories(model_tar: str) -> dict:
    with tarfile.open(model_tar, "r:gz") as tf:
        pt2acc = pd.read_csv(tf.extractfile("pt2acc.txt"), sep="\t", index_col=0)
    s = pt2acc.set_index("accession")["category"]
    return {str(k): str(v) for k, v in s.items()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--primary", required=True, help="phypat per-batch majority_vote.tsv")
    ap.add_argument("--secondary", required=True, help="phypat+PGL per-batch majority_vote.tsv")
    ap.add_argument("--model-tar", required=True, help="phypat.tar.gz (for trait categories)")
    ap.add_argument("--primary-name", default="phypat")
    ap.add_argument("--secondary-name", default="phypat+PGL")
    ap.add_argument("--outdir", required=True, help="directory for the per-genome files")
    ap.add_argument("--marker", required=True, help="sentinel file written on completion")
    ap.add_argument("--suffix", default=".majority_vote.tsv")
    args = ap.parse_args(argv)

    p = pd.read_csv(args.primary, sep="\t", index_col=0)
    s = pd.read_csv(args.secondary, sep="\t", index_col=0)
    p.index = p.index.astype(str)
    s.index = s.index.astype(str)
    traits = list(p.columns)
    s = s.reindex(index=p.index, columns=traits, fill_value=0)

    cats = load_categories(args.model_tar)
    cat_col = [cats.get(t, "") for t in traits]

    P = p.to_numpy()
    S = s.to_numpy()
    combined = P * 1 + S * 2  # 0/1/2/3

    os.makedirs(args.outdir, exist_ok=True)
    header = ["phenotype", "category", args.primary_name, args.secondary_name, "combined"]
    for i, genome in enumerate(p.index):
        path = os.path.join(args.outdir, f"{genome}{args.suffix}")
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header)
            for j, t in enumerate(traits):
                w.writerow([t, cat_col[j], int(P[i, j]), int(S[i, j]), int(combined[i, j])])

    os.makedirs(os.path.dirname(os.path.abspath(args.marker)) or ".", exist_ok=True)
    with open(args.marker, "w") as fh:
        fh.write(f"{len(p.index)} genomes written to {args.outdir}\n")
    print(f"[split_per_genome] wrote {len(p.index)} files -> {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
