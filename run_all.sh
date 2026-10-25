#!/usr/bin/env bash
# Reproduce every Phase-1 result from scratch (~4-5 minutes on 1 CPU core).
set -euo pipefail
cd "$(dirname "$0")"
python3 -m pytest -q tests
cd experiments
python3 exp_a_geometry.py --perms 200
python3 exp_b_hypotheses.py
python3 validate_curve_level_null.py
