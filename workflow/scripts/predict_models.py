#!/usr/bin/env python3
"""Predict traitar phenotypes for a batch of genomes from per-genome Pfam counts.

Vectorised, dependency-light reimplementation of traitar3's
``predict.majority_predict`` + ``predict.aggregate`` for a single phenotype-model
archive. It produces the same numbers as traitar3 but:

  * uses numpy instead of the removed ``pandas.np`` API, so it runs on current
    pandas/numpy;
  * scores a whole batch of genomes in one matrix multiply, and loads the model
    once per batch (instead of once per genome), which is what makes 300k
    genomes tractable.

Decision value of voter ``i`` for phenotype ``p`` on genome ``g``::

    d = bias[p, i] + sum_f  W[p, f, i] * 1[count(g, f) > 0]

(annotation counts are binarised to presence/absence, exactly as traitar does
with ``test_data_n = (test_data > 0)``). ``single-votes`` is the number of voters
with ``d > 0`` (0..k). ``majority-vote`` is ``single-votes >= threshold``.

On the majority threshold
-------------------------
traitar3's literal Python-3 code computes ``votes >= k/2 + 1`` which, with
``k = 5``, is ``votes >= 3.5`` (i.e. 4 of 5 voters) -- a regression from the
Python-2 integer division it was ported from. The original traitar and the
reference outputs shipped in the repo use a true majority of ``>= 3`` of 5.
The default here is the true majority ``floor(k/2) + 1``; pass
``--majority-threshold 3.5`` (or ``--literal-traitar3``) to reproduce the
literal traitar3 numbers.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import tarfile

import numpy as np
import pandas as pd


class PhenotypeModel:
    """Reader for a traitar phenotype-collection ``.tar.gz`` archive."""

    def __init__(self, tar_path: str, k: int):
        self.path = tar_path
        self.k = k
        self.tar = tarfile.open(tar_path, "r:gz")
        cfg = pd.read_csv(self.tar.extractfile("config.txt"), sep="\t", index_col=0)
        self.name = str(cfg.loc["archive_name", "value"])
        pt2acc = pd.read_csv(self.tar.extractfile("pt2acc.txt"), sep="\t", index_col=0)
        pt2acc.index = pt2acc.index.astype(str)
        self.pt2acc = pt2acc                      # index = pt id, col 0 = phenotype name
        self.pt_ids = list(pt2acc.index)          # preserves traitar's column order

    def _bias(self, pt: str) -> np.ndarray:
        b = pd.read_csv(self.tar.extractfile(f"{pt}_bias.txt"),
                        sep="\t", index_col=0, header=None)
        return b.iloc[: self.k, 0].to_numpy(dtype=float)

    def _feats(self, pt: str) -> pd.DataFrame:
        f = pd.read_csv(self.tar.extractfile(f"{pt}_feats.txt"), sep="\t", index_col=0)
        f.index = f.index.astype(str)
        return f.iloc[:, : self.k]

    def load(self):
        """Return (pt_ids, pt_names, feature_index, W, B).

        W: (n_pt, n_feat, k) weights;  B: (n_pt, k) biases.
        Phenotypes whose model files are missing from the archive are skipped
        (matching traitar's KeyError handling).
        """
        kept_ids, biases, feats = [], [], []
        for pt in self.pt_ids:
            try:
                b = self._bias(pt)
                f = self._feats(pt)
            except KeyError:
                sys.stderr.write(f"  [{self.name}] no model for phenotype {pt}; skipping\n")
                continue
            kept_ids.append(pt)
            biases.append(b)
            feats.append(f)

        # shared feature universe across all phenotypes (sorted for determinism)
        feature_index = sorted(set().union(*[set(f.index) for f in feats]))
        pos = {acc: i for i, acc in enumerate(feature_index)}
        n_pt, n_feat = len(kept_ids), len(feature_index)

        W = np.zeros((n_pt, n_feat, self.k), dtype=float)
        B = np.zeros((n_pt, self.k), dtype=float)
        for p, (f, b) in enumerate(zip(feats, biases)):
            rows = [pos[acc] for acc in f.index]
            W[p, rows, :] = f.to_numpy(dtype=float)
            B[p, : len(b)] = b

        pt_names = [str(self.pt2acc.loc[pt].iloc[0]) for pt in kept_ids]
        return kept_ids, pt_names, feature_index, W, B


def read_counts_matrix(counts_files: list[str], names: list[str],
                       feature_index: list[str]) -> np.ndarray:
    """Binary presence/absence matrix (n_genomes, n_feat) from per-genome counts."""
    pos = {acc: i for i, acc in enumerate(feature_index)}
    X = np.zeros((len(counts_files), len(feature_index)), dtype=np.float32)
    for g, path in enumerate(counts_files):
        with open(path) as fh:
            header = fh.readline()  # pfam_acc<TAB>count
            for line in fh:
                acc, _, _count = line.partition("\t")
                j = pos.get(acc.strip())
                if j is not None:
                    X[g, j] = 1.0  # presence/absence
    return X


def matrix_from_summary(summary_path: str, feature_index: list[str]
                        ) -> tuple[list[str], np.ndarray]:
    """Build the binary matrix from a traitar-style summary.dat (sample x Pfam)."""
    m = pd.read_csv(summary_path, sep="\t", index_col=0)
    m = m.reindex(columns=feature_index, fill_value=0)
    X = (m.to_numpy() > 0).astype(np.float32)
    return [str(i) for i in m.index], X


def score(X: np.ndarray, W: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Voter decision values -> (n_genomes, n_pt, k)."""
    n_pt, n_feat, k = W.shape
    w_flat = W.transpose(0, 2, 1).reshape(n_pt * k, n_feat)   # (P*k, F)
    d = X.astype(float) @ w_flat.T                            # (G, P*k)
    d += B.reshape(1, n_pt * k)
    return d.reshape(X.shape[0], n_pt, k)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, help="phenotype model .tar.gz")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--counts-list", help="file listing per-genome Pfam count TSVs (one path per line)")
    src.add_argument("--matrix", help="traitar-style summary.dat (sample x Pfam counts)")
    ap.add_argument("--out-prefix", required=True,
                    help="output prefix; writes <prefix>.single_votes.tsv and <prefix>.majority_vote.tsv")
    ap.add_argument("--voters", type=int, default=5, help="committee size k (default: %(default)s)")
    ap.add_argument("--majority-threshold", type=float, default=None,
                    help="min positive votes for a majority call (default: floor(k/2)+1)")
    ap.add_argument("--literal-traitar3", action="store_true",
                    help="use traitar3's literal threshold k/2+1 (== 3.5 for k=5)")
    ap.add_argument("--write-raw", action="store_true", help="also write <prefix>.raw.tsv decision values")
    args = ap.parse_args(argv)

    k = args.voters
    if args.literal_traitar3:
        threshold = k / 2 + 1
    elif args.majority_threshold is not None:
        threshold = args.majority_threshold
    else:
        threshold = math.floor(k / 2) + 1

    model = PhenotypeModel(args.model, k)
    pt_ids, pt_names, feature_index, W, B = model.load()

    if args.counts_list:
        with open(args.counts_list) as fh:
            counts_files = [ln.strip() for ln in fh if ln.strip()]
        genomes = [_sample_name(p) for p in counts_files]
        X = read_counts_matrix(counts_files, genomes, feature_index)
    else:
        genomes, X = matrix_from_summary(args.matrix, feature_index)

    decisions = score(X, W, B)                       # (G, P, k)
    votes = (decisions > 0).sum(axis=2).astype(int)  # (G, P)

    sv = pd.DataFrame(votes, index=genomes, columns=pt_names)
    mv = (sv >= threshold).astype(int)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)) or ".", exist_ok=True)
    sv.to_csv(args.out_prefix + ".single_votes.tsv", sep="\t")
    mv.to_csv(args.out_prefix + ".majority_vote.tsv", sep="\t")
    if args.write_raw:
        cols = [f"{pt}_{i}" for pt in pt_ids for i in range(k)]
        raw = pd.DataFrame(decisions.reshape(len(genomes), -1), index=genomes, columns=cols)
        raw.to_csv(args.out_prefix + ".raw.tsv", sep="\t", float_format="%.3f")

    sys.stderr.write(
        f"[{model.name}] {len(genomes)} genomes x {len(pt_names)} phenotypes; "
        f"k={k} threshold>={threshold}\n"
    )
    return 0


def _sample_name(path: str) -> str:
    """Recover the genome/sample name from a counts file path.

    Counts files are named ``<sample>.pfam_counts.tsv`` by the workflow; fall
    back to stripping a single extension for other inputs.
    """
    base = os.path.basename(path)
    if base.endswith(".pfam_counts.tsv"):
        return base[: -len(".pfam_counts.tsv")]
    stem, _, _ = base.partition(".pfam_counts")
    if stem != base:
        return stem
    return os.path.splitext(base)[0]


if __name__ == "__main__":
    raise SystemExit(main())
