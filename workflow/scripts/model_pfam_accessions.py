#!/usr/bin/env python3
"""Write the Pfam accessions that actually affect traitar's phenotype calls.

The phenotype predictors are linear models: a decision value is
``bias + sum_f weight[f] * presence[f]``. A Pfam family whose weight is zero in
every voter of every phenotype cannot change any decision, so it never needs to
be searched. Restricting ``hmmsearch`` to the families with a non-zero weight
therefore gives byte-for-byte identical predictions while scanning far fewer
profiles.

(The models' ``pf2acc_desc.txt`` lists the full Pfam 27.0 feature space -- all
~14,831 families -- so it is *not* a useful subset. The non-zero-weight union of
the ``*_feats.txt`` predictors is much smaller, typically a few thousand.)

Accessions are written unversioned (PF00389, not PF00389.25), one per line.
"""
from __future__ import annotations

import argparse
import sys
import tarfile

import pandas as pd


def nonzero_accessions(tar_path: str) -> set[str]:
    """Pfam accessions with a non-zero weight in any voter of any phenotype."""
    accs: set[str] = set()
    with tarfile.open(tar_path, "r:gz") as tf:
        pt2acc = pd.read_csv(tf.extractfile("pt2acc.txt"), sep="\t", index_col=0)
        pt2acc.index = pt2acc.index.astype(str)
        for pt in pt2acc.index:
            try:
                feats = pd.read_csv(tf.extractfile(f"{pt}_feats.txt"), sep="\t", index_col=0)
            except KeyError:
                continue
            feats.index = feats.index.astype(str)
            nz = feats.index[(feats != 0).any(axis=1)]
            accs |= {str(x).split(".")[0] for x in nz}
    return accs


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model_tar", nargs="+", help="one or more model .tar.gz archives")
    ap.add_argument("-o", "--out", required=True, help="output accession list (one per line)")
    args = ap.parse_args(argv)

    accs: set[str] = set()
    for tar_path in args.model_tar:
        accs |= nonzero_accessions(tar_path)

    with open(args.out, "w") as fh:
        for acc in sorted(accs):
            fh.write(acc + "\n")
    sys.stderr.write(
        f"{len(accs)} Pfam families with non-zero predictor weight -> {args.out}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
