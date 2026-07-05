#!/usr/bin/env python3
"""Combine phypat and phypat+PGL predictions into traitar's joint calls.

Faithful reimplementation of ``traitar3.merge_preds.comb_preds``.

Combined majority-vote encoding (per phenotype shared by both models):
    0  negative in both
    1  positive in the primary model only          (phypat)
    2  positive in the secondary model only         (phypat+PGL)
    3  positive in both
For a phenotype present in only one model, traitar emits 2*primary or
1*secondary respectively (kept here for fidelity; with the shipped models both
collections cover the same 67 phenotypes, so every phenotype is "shared").

Combined single-votes is the elementwise sum of the two models' vote counts.
Flat (long) versions list every non-zero call as one row.
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd


def _read(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df.fillna(0)


def _align(a: pd.DataFrame, b: pd.DataFrame):
    idx = list(a.index) + [i for i in b.index if i not in set(a.index)]
    cols = sorted(set(a.columns) | set(b.columns))
    A = a.reindex(index=idx, columns=cols, fill_value=0)
    B = b.reindex(index=idx, columns=cols, fill_value=0)
    return A, B, idx, cols


def combine_majority(m1: pd.DataFrame, m2: pd.DataFrame) -> pd.DataFrame:
    M1, M2, idx, cols = _align(m1, m2)
    only1 = [c for c in cols if c in m1.columns and c not in m2.columns]
    only2 = [c for c in cols if c in m2.columns and c not in m1.columns]
    shared = [c for c in cols if c in m1.columns and c in m2.columns]
    comb = pd.DataFrame(0.0, index=idx, columns=cols)
    if shared:
        # (1,1)->3, (1,0)->1, (0,1)->2, (0,0)->0
        comb[shared] = M1[shared] * 1 + M2[shared] * 2
    if only1:
        comb[only1] = M1[only1] * 2
    if only2:
        comb[only2] = M2[only2] * 1
    return comb


def combine_single(s1: pd.DataFrame, s2: pd.DataFrame) -> pd.DataFrame:
    S1, S2, _idx, _cols = _align(s1, s2)
    return S1 + S2


def _long_nonzero(df: pd.DataFrame, name: str) -> pd.DataFrame:
    """Every non-zero cell of df as rows (sample, phenotype, score, model).

    Vectorised with numpy.nonzero instead of a per-cell .loc loop, so it stays
    fast at hundreds of thousands of genomes (the loop version is ~O(genomes x
    phenotypes) label lookups and does not scale).
    """
    vals = df.to_numpy()
    r, c = np.nonzero(vals)
    return pd.DataFrame({
        "sample": np.asarray(df.index)[r],
        "phenotype": np.asarray(df.columns)[c],
        "score": vals[r, c],
        "phenotype_model": name,
    })


_FLAT_COLS = ["sample", "phenotype", "score", "phenotype_model"]


def flatten(df1: pd.DataFrame, df2: pd.DataFrame, name1: str, name2: str,
            out: str, full: bool = True) -> None:
    if not full:
        # header-only placeholder so the declared output still exists
        pd.DataFrame(columns=_FLAT_COLS).to_csv(out, sep="\t", index=False)
        return
    long = pd.concat([_long_nonzero(df1, name1), _long_nonzero(df2, name2)],
                     ignore_index=True)
    long.to_csv(out, sep="\t", index=False)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--primary-majority", required=True)
    ap.add_argument("--primary-single", required=True)
    ap.add_argument("--secondary-majority", required=True)
    ap.add_argument("--secondary-single", required=True)
    ap.add_argument("--primary-name", default="phypat")
    ap.add_argument("--secondary-name", default="phypat+PGL")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--no-flat", action="store_true",
                    help="write header-only placeholders for the long-format flat "
                         "files instead of the full tables (the single-votes flat is "
                         "tens of millions of rows / several GB at many genomes)")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    m1_maj, m2_maj = _read(args.primary_majority), _read(args.secondary_majority)
    m1_sv, m2_sv = _read(args.primary_single), _read(args.secondary_single)

    out = os.path.join(args.out_dir, "predictions_majority-vote_combined.tsv")
    combine_majority(m1_maj, m2_maj).to_csv(out, sep="\t")
    out = os.path.join(args.out_dir, "predictions_single-votes_combined.tsv")
    combine_single(m1_sv, m2_sv).to_csv(out, sep="\t")
    full = not args.no_flat
    flatten(m1_maj, m2_maj, args.primary_name, args.secondary_name,
            os.path.join(args.out_dir, "predictions_flat_majority-votes_combined.tsv"), full=full)
    flatten(m1_sv, m2_sv, args.primary_name, args.secondary_name,
            os.path.join(args.out_dir, "predictions_flat_single-votes_combined.tsv"), full=full)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
