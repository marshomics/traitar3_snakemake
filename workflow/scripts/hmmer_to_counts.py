#!/usr/bin/env python3
"""Convert one hmmsearch --domtblout file into a per-genome Pfam count vector.

Faithful reproduction of the two traitar3 steps that turn raw HMMER output into
the per-genome row of the annotation summary matrix:

  * ``hmmer2filtered_best.filter_pfam`` keeps domain hits whose independent
    domain E-value (``i-Evalue``) is <= 0.01 AND whose domain bit score is
    >= 25.
  * ``hmmer2filtered_best.aggregate_domain_hits`` keeps, for each
    (protein, Pfam-family) pair, only the single best-scoring domain hit.
  * ``domtblout2gene_generic.gene2hmm`` then counts, per Pfam accession, how
    many proteins retained a hit (``summary.dat`` cell value).

The output is the sparse equivalent of one ``summary.dat`` row: a TSV with
columns ``pfam_acc<TAB>count`` for every Pfam accession with count > 0.

HMMER3 ``--domtblout`` column indices (0-based) used here:
    0  target name   (protein id)
    3  query name    (Pfam family name)
    4  query acc     (Pfam accession, e.g. PF00389.25)
   12  i-Evalue      (independent, this-domain E-value)
   13  score         (this-domain bit score)
"""
from __future__ import annotations

import argparse
import csv
import sys

COL_TARGET = 0
COL_QUERY_ACC = 4
COL_I_EVALUE = 12
COL_DOM_SCORE = 13
_MIN_FIELDS = COL_DOM_SCORE + 1  # need indices 0..13


def count_pfams(domtbl_path: str, evalue_max: float, score_min: float) -> dict[str, int]:
    """Return {pfam_accession: n_proteins_with_a_best_hit}."""
    # best score seen for each (protein, pfam) that passes the thresholds
    best: dict[tuple[str, str], float] = {}
    with open(domtbl_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            # columns 0..13 are single whitespace-delimited tokens; a plain
            # split is safe (only the trailing description column has spaces).
            f = line.split()
            if len(f) < _MIN_FIELDS:
                continue
            try:
                ievalue = float(f[COL_I_EVALUE])
                score = float(f[COL_DOM_SCORE])
            except ValueError:
                continue
            if ievalue > evalue_max or score < score_min:
                continue
            gene = f[COL_TARGET]
            pfam = f[COL_QUERY_ACC].split(".")[0]
            key = (gene, pfam)
            prev = best.get(key)
            if prev is None or score > prev:
                best[key] = score

    counts: dict[str, int] = {}
    for _gene, pfam in best:
        counts[pfam] = counts.get(pfam, 0) + 1
    return counts


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("domtbl", help="hmmsearch --domtblout file for one genome")
    ap.add_argument("out_tsv", help="output Pfam count TSV (pfam_acc<TAB>count)")
    ap.add_argument("--evalue-max", type=float, default=0.01,
                    help="max i-Evalue to keep a domain hit (default: %(default)s)")
    ap.add_argument("--score-min", type=float, default=25.0,
                    help="min domain bit score to keep a hit (default: %(default)s)")
    args = ap.parse_args(argv)

    counts = count_pfams(args.domtbl, args.evalue_max, args.score_min)
    with open(args.out_tsv, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["pfam_acc", "count"])
        for pfam in sorted(counts):
            w.writerow([pfam, counts[pfam]])
    sys.stderr.write(
        f"{args.domtbl}: {len(counts)} Pfam families with hits -> {args.out_tsv}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
