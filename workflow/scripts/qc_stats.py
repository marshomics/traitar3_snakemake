#!/usr/bin/env python3
"""Per-genome annotation QC statistics for one batch of genomes.

Reads the small per-genome Pfam count files (already produced by the annotate
step; this does not touch hmmsearch) and, for each genome, records:

  * n_pfam_families    - distinct Pfam families with >=1 hit
  * total_pfam_hits    - sum of per-family protein counts (annotation depth)
  * n_model_features   - families that are in the phenotype models' feature set
  * frac_model_coverage- n_model_features / (size of the model feature set)

`frac_model_coverage` is the QC metric that matters most for this pipeline: a
genome covering few of the model's Pfam families has little for the classifiers
to act on, so its phenotype calls are less reliable (the Traitar paper shows
accuracy drops with genome incompleteness).
"""
from __future__ import annotations

import argparse
import csv
import os


def sample_name(path: str) -> str:
    base = os.path.basename(path)
    for suf in (".pfam_counts.tsv", ".tsv"):
        if base.endswith(suf):
            return base[: -len(suf)]
    return os.path.splitext(base)[0]


def stats_for_file(path: str, model_set: set[str]) -> tuple[int, int, int]:
    n_fam = 0
    total = 0
    n_model = 0
    with open(path) as fh:
        fh.readline()  # header: pfam_acc<TAB>count
        for line in fh:
            if not line.strip():
                continue
            acc, _, cnt = line.partition("\t")
            acc = acc.strip()
            n_fam += 1
            try:
                total += int(cnt)
            except ValueError:
                pass
            if acc in model_set:
                n_model += 1
    return n_fam, total, n_model


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts-list", required=True, help="file of per-genome count TSV paths")
    ap.add_argument("--model-accessions", required=True,
                    help="model_pfam_accessions.txt (one Pfam accession per line)")
    ap.add_argument("--out", required=True, help="output per-batch stats TSV")
    args = ap.parse_args(argv)

    with open(args.model_accessions) as fh:
        model_set = {ln.strip().split(".")[0] for ln in fh if ln.strip()}
    n_model_total = max(len(model_set), 1)

    with open(args.counts_list) as fh:
        paths = [ln.strip() for ln in fh if ln.strip()]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample", "n_pfam_families", "total_pfam_hits",
                    "n_model_features", "frac_model_coverage"])
        for p in paths:
            n_fam, total, n_model = stats_for_file(p, model_set)
            w.writerow([sample_name(p), n_fam, total, n_model,
                        round(n_model / n_model_total, 6)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
