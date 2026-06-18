#!/usr/bin/env bash
# Run Primer3 on the prepared input and report design diagnostics.
set -euo pipefail
cd /Users/fvrodriguez/Desktop/primermu

primer3_core work/p3_in.txt > work/p3_out.txt 2> work/p3_err.txt
echo "primer3 exit ok"
echo "--- stderr (if any) ---"
head work/p3_err.txt || true
echo "--- explain / counts ---"
grep -E '^PRIMER_PAIR_NUM_RETURNED|_EXPLAIN' work/p3_out.txt || true
