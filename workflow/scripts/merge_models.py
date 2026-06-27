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


def flatten(df1: pd.DataFrame, df2: pd.DataFrame, name1: str, name2: str, out: str) -> None:
    rows = []
    genomes = list(df1.index) + [i for i in df2.index if i not in set(df1.index)]
    phenos = sorted(set(df1.columns) | set(df2.columns))
    for g in genomes:
        for p in phenos:
            if g in df1.index and p in df1.columns and df1.loc[g, p] != 0:
                rows.append((g, p, df1.loc[g, p], name1))
            if g in df2.index and p in df2.columns and df2.loc[g, p] != 0:
                rows.append((g, p, df2.loc[g, p], name2))
    pd.DataFrame(rows, columns=["sample", "phenotype", "score", "phenotype_model"]).to_csv(
        out, sep="\t", index=False
    )


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
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    m1_maj, m2_maj = _read(args.primary_majority), _read(args.secondary_majority)
    m1_sv, m2_sv = _read(args.primary_single), _read(args.secondary_single)

    out = os.path.join(args.out_dir, "predictions_majority-vote_combined.tsv")
    combine_majority(m1_maj, m2_maj).to_csv(out, sep="\t")
    out = os.path.join(args.out_dir, "predictions_single-votes_combined.tsv")
    combine_single(m1_sv, m2_sv).to_csv(out, sep="\t")
    flatten(m1_maj, m2_maj, args.primary_name, args.secondary_name,
            os.path.join(args.out_dir, "predictions_flat_majority-votes_combined.tsv"))
    flatten(m1_sv, m2_sv, args.primary_name, args.secondary_name,
            os.path.join(args.out_dir, "predictions_flat_single-votes_combined.tsv"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
