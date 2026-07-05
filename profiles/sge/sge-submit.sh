#!/usr/bin/env bash
# Snakemake cluster-generic submit wrapper for SGE / UGE.
# Translates Snakemake resources into a qsub command and prints the SGE job id.
set -euo pipefail

THREADS=1
MEM_MB=4000
RUNTIME=60          # minutes
RULE="smk"
JOBSCRIPT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --rule)     RULE="$2";    shift 2 ;;
        --threads)  THREADS="$2"; shift 2 ;;
        --mem-mb)   MEM_MB="$2";  shift 2 ;;
        --runtime)  RUNTIME="$2"; shift 2 ;;
        *)          JOBSCRIPT="$1"; shift ;;   # last positional arg = job script
    esac
done

if [[ -z "$JOBSCRIPT" ]]; then
    echo "sge-submit.sh: no job script given" >&2
    exit 1
fi

THREADS=${THREADS:-1}
[[ "$THREADS" -lt 1 ]] && THREADS=1

# SGE h_vmem is requested per slot, so divide total memory by the slot count.
MEM_PER_SLOT=$(( (MEM_MB + THREADS - 1) / THREADS ))

# runtime minutes -> HH:MM:SS
HH=$(( RUNTIME / 60 ))
MM=$(( RUNTIME % 60 ))
HRT=$(printf '%02d:%02d:00' "$HH" "$MM")

PE="${SGE_PE:-smp}"          # parallel environment name (override via SGE_PE)
LOGDIR="logs/sge"
mkdir -p "$LOGDIR"

# -terse prints just the job id; that is what Snakemake reads from stdout.
# Build -pe only when more than one slot is requested (some PEs reject 1).
PE_ARG=()
if [[ "$THREADS" -gt 1 ]]; then
    PE_ARG=(-pe "$PE" "$THREADS")
fi

jobid=$(qsub -terse -cwd -V \
    -v OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,NUMEXPR_NUM_THREADS=1 \
    -N "smk.${RULE}" \
    "${PE_ARG[@]}" \
    -l h_vmem="${MEM_PER_SLOT}M" \
    -l h_rt="${HRT}" \
    -o "$LOGDIR" -e "$LOGDIR" \
    "$JOBSCRIPT")

# strip any whitespace/newline
echo "${jobid//[[:space:]]/}"
