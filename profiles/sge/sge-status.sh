#!/usr/bin/env bash
# Snakemake cluster-generic status wrapper for SGE / UGE.
# Usage: sge-status.sh <jobid>   ->  prints "running" | "success" | "failed"
set -uo pipefail

jobid="$1"

# 1) Still pending/running? qstat -j returns 0 while the job is known to the
#    scheduler. Detect the error state 'E' explicitly.
if qstat -j "$jobid" >/dev/null 2>&1; then
    state=$(qstat 2>/dev/null | awk -v id="$jobid" '$1==id {print $5; exit}')
    case "$state" in
        *E*) echo "failed" ;;     # Eqw / Error
        *)   echo "running" ;;    # qw, r, t, hqw, ...
    esac
    exit 0
fi

# 2) Job has left the queue -> consult accounting for the exit status.
if command -v qacct >/dev/null 2>&1; then
    # qacct can lag a few seconds after a job ends; retry briefly.
    for _ in 1 2 3 4 5; do
        acct=$(qacct -j "$jobid" 2>/dev/null) && [[ -n "$acct" ]] && break
        sleep 2
    done
    if [[ -n "${acct:-}" ]]; then
        exit_status=$(awk '/^exit_status/ {print $2; exit}' <<<"$acct")
        failed=$(awk '/^failed/ {print $2; exit}' <<<"$acct")
        if [[ "${exit_status:-1}" == "0" && "${failed:-1}" == "0" ]]; then
            echo "success"
        else
            echo "failed"
        fi
        exit 0
    fi
fi

# 3) No accounting available: the job is gone from the queue and we can't prove
#    failure. Report success and rely on Snakemake's output-file check (it marks
#    the job failed if the expected outputs are missing).
echo "success"
