#!/usr/bin/env python3
"""Derive a binary call set from traitar's combined 0/1/2/3 majority-vote matrix.

The combined matrix encodes, per genome x trait:
    0 negative in both models
    1 positive in phypat only
    2 positive in phypat+PGL only
    3 positive in both
A "call set" selects a subset of those codes and emits 1 where the cell is in the
set, else 0. For example the phypat+PGL "reliable" set is codes {2, 3}; the
strict both-agree set is {3}; the phypat "sensitive" set is {1, 3}.

Outputs a genome x trait 0/1 matrix and, optionally, a long/tidy list of the
positive calls.
"""
from __future__ import annotations

import argparse
import os

import pandas as pd


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--combined", required=True,
                    help="predictions_majority-vote_combined.tsv (0/1/2/3 matrix)")
    ap.add_argument("--codes", required=True,
                    help="comma-separated codes that count as positive, e.g. 2,3")
    ap.add_argument("--out", required=True, help="output genome x trait 0/1 matrix")
    ap.add_argument("--out-long", default=None,
                    help="optional tidy list of positive calls (sample, phenotype, call_set)")
    ap.add_argument("--name", default="call_set", help="label used in the long output")
    args = ap.parse_args(argv)

    codes = {int(c) for c in args.codes.split(",") if c.strip() != ""}

    m = pd.read_csv(args.combined, sep="\t", index_col=0)
    m.index = m.index.astype(str)
    codes_int = m.fillna(0).round().astype(int)
    binary = codes_int.isin(codes).astype(int)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    binary.to_csv(args.out, sep="\t")

    if args.out_long:
        rows = []
        for genome in binary.index:
            for trait in binary.columns:
                if binary.loc[genome, trait] == 1:
                    rows.append((genome, trait, args.name))
        pd.DataFrame(rows, columns=["sample", "phenotype", "call_set"]).to_csv(
            args.out_long, sep="\t", index=False
        )

    n_pos = int(binary.to_numpy().sum())
    print(f"[{args.name}] codes={sorted(codes)}: {n_pos} positive calls across "
          f"{binary.shape[0]} genomes x {binary.shape[1]} traits -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
