#!/usr/bin/env python3
"""Concatenate per-batch prediction tables into one table (genomes x phenotypes).

All batch tables for a given model + vote type share the same phenotype columns,
so this is a row-wise concatenation. Reads paths from a manifest file to stay
well under command-line length limits at large genome counts.
"""
from __future__ import annotations

import argparse

import pandas as pd


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs-list", required=True,
                    help="file listing per-batch table paths (one per line)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    with open(args.inputs_list) as fh:
        paths = [ln.strip() for ln in fh if ln.strip()]

    frames = [pd.read_csv(p, sep="\t", index_col=0) for p in paths]
    full = pd.concat(frames, axis=0)
    # stable column order across batches (union; missing -> 0)
    cols = sorted(set().union(*[set(f.columns) for f in frames]))
    full = full.reindex(columns=cols, fill_value=0)
    full.to_csv(args.out, sep="\t")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
