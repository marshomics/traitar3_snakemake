#!/usr/bin/env bash
# Offline correctness check (no HMMER / Pfam needed).
#
# Drives the prediction + merge scripts from traitar3's own committed
# annotation summary and compares the result to the committed reference
# predictions. Exits non-zero on any mismatch.
#
# Run from the repository root:  bash test/run_offline_check.sh
set -euo pipefail
cd "$(dirname "$0")/.."

SC=workflow/scripts
EXP=test/expected
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "[1/3] predict from committed summary.dat (both models)"
python3 "$SC/predict_models.py" --model resources/models/phypat.tar.gz \
    --matrix "$EXP/summary.dat" --out-prefix "$TMP/phypat"
python3 "$SC/predict_models.py" --model "resources/models/phypat+PGL.tar.gz" \
    --matrix "$EXP/summary.dat" --out-prefix "$TMP/phypatPGL"

echo "[2/3] merge -> combined calls"
python3 "$SC/merge_models.py" \
    --primary-majority "$TMP/phypat.majority_vote.tsv" --primary-single "$TMP/phypat.single_votes.tsv" \
    --secondary-majority "$TMP/phypatPGL.majority_vote.tsv" --secondary-single "$TMP/phypatPGL.single_votes.tsv" \
    --out-dir "$TMP"

echo "[3/3] compare to committed reference outputs"
python3 - "$TMP" "$EXP" <<'PY'
import sys, pandas as pd
tmp, exp = sys.argv[1], sys.argv[2]

def load(p):
    df = pd.read_csv(p, sep='\t', index_col=0); df.index = df.index.astype(str); return df

# combined majority-vote must match exactly
mine = load(f"{tmp}/predictions_majority-vote_combined.tsv")
ref  = load(f"{exp}/predictions_majority-vote_combined.txt")
cols = sorted(set(mine.columns) & set(ref.columns))
assert set(mine.columns) == set(ref.columns), "phenotype set differs"
a = mine.reindex(index=ref.index, columns=cols).astype(float)
b = ref.reindex(columns=cols).astype(float)
ndiff = int((a.values != b.values).sum())
print(f"  combined majority-vote: {a.shape[0]} genomes x {len(cols)} phenotypes, diffs={ndiff}")
assert ndiff == 0, "combined majority-vote mismatch"

print("OK: combined calls reproduce the traitar3 reference exactly.")
PY
echo "PASSED"
