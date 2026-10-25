#!/usr/bin/env bash
# Reproduce the Phase-2 results (endomorphism search, threshold test, C3 fit + held-out).
# Uses plain python; the recorded evidence runs in results/phase2/run-00N/ were made with the
# computation-audit bounded runner (see each manifest.json for the exact argv and input hashes).
set -euo pipefail
cd "$(dirname "$0")"
OUT=${1:-results/phase2/reproduce}
python3 -m pytest -q tests
python3 experiments/exp_c_endomorphisms.py --dmax 20 --dmax163 42 --out "$OUT/run-001"
python3 experiments/exp_c2_threshold.py --out "$OUT/run-005"
python3 experiments/exp_c3_small_subgroups.py --out "$OUT/run-006" --label fit
python3 experiments/exp_c3_small_subgroups.py --out "$OUT/run-007" --label heldout \
  --qs 67,79,89,103,113,127,131,137,139,149,151,163,167,173,179,181,191,193,197,199 \
  --curves-per-q 1 --seed 20261029 --cm=-3,-4,-7
