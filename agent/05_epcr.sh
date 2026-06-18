#!/usr/bin/env bash
# Genome-wide in-silico PCR specificity for the 9 GSP F x R combos.
# A pair is specific if it yields exactly ONE amplicon (at the target locus).
set -euo pipefail
ROOT=/Users/fvrodriguez/Desktop/primermu
cd "$ROOT/data"   # re-PCR finds W22.famap (basename stored in hash) relative to CWD
echo "=== re-PCR: 1 mismatch/primer, 1 gap, size margin 200 ==="
re-PCR -S W22.hash -n 1 -g 1 -m 200 -o "$ROOT/results/epcr_raw.txt" "$ROOT/work/epcr_sts.txt" 2> "$ROOT/work/epcr.log" || true
echo "--- log ---"; tail -3 "$ROOT/work/epcr.log" || true
echo "--- raw hits ---"
cat "$ROOT/results/epcr_raw.txt"
