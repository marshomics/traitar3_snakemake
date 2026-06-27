#!/usr/bin/env python3
"""Subset a Pfam-A.hmm file to a given set of accessions.

Streams the (large) HMMER3 ASCII database one profile record at a time and keeps
only records whose ``ACC`` (compared without the version suffix) is in the
requested set. Because each profile carries its own gathering threshold (used by
``hmmsearch --cut_ga``), the hits reported for a kept family are byte-for-byte the
same whether the search runs against the full or the subset database.

A profile record runs from a ``HMMER3/`` header line to the next ``//`` line.
Memory use is one record at a time, so this handles the ~1.2 GB Pfam-A.hmm fine.
"""
from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("acc_list", help="file of Pfam accessions to keep (one per line, unversioned)")
    ap.add_argument("hmm_in", help="input Pfam-A.hmm")
    ap.add_argument("hmm_out", help="output subset HMM")
    args = ap.parse_args(argv)

    with open(args.acc_list) as fh:
        keep = {line.strip().split(".")[0] for line in fh if line.strip()}

    n_kept = 0
    n_total = 0
    buf: list[str] = []
    acc: str | None = None
    with open(args.hmm_in) as fin, open(args.hmm_out, "w") as fout:
        for line in fin:
            buf.append(line)
            if line.startswith("ACC "):
                acc = line.split()[1].split(".")[0]
            elif line.startswith("//"):
                n_total += 1
                if acc is not None and acc in keep:
                    fout.writelines(buf)
                    n_kept += 1
                buf = []
                acc = None

    sys.stderr.write(
        f"kept {n_kept}/{n_total} profiles ({len(keep)} requested) -> {args.hmm_out}\n"
    )
    if n_kept == 0:
        sys.stderr.write("ERROR: no profiles matched; check the accession list / HMM file\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
