#!/usr/bin/env bash
# Forward e-PCR: hash the primer pairs, stream the whole W22 genome FASTA.
# Specific pair == exactly ONE amplicon (the target locus).
set -euo pipefail
ROOT=/Users/fvrodriguez/Desktop/primermu
cd "$ROOT"
echo "=== e-PCR (forward): w7, 1 mismatch, 1 gap, margin 200, tabular ==="
# POSIX options MUST precede the stsfile; trailing args are all treated as FASTA inputs.
e-PCR -n 1 -g 1 -m 200 -w 7 -t 3 -o results/epcr_fwd.txt work/epcr_sts.txt data/W22.fa
echo "--- hits ---"
cat results/epcr_fwd.txt
echo "--- amplicon count per primer pair (expect 1 each) ---"
grep -v '^#' results/epcr_fwd.txt | awk '{print $1}' | sort | uniq -c
