# ECC Mathematical Research Lab — Phase 1: Mathematical Structure Exploration

An experimental-mathematics lab for studying the structure of the scalar-multiplication
sequence G, 2G, 3G, … on **small synthetic elliptic curves over prime fields**
(p ≈ 2^16–2^17). The goal of this phase is understanding which structure exists, which of it is
already explained by known theory, and what deserves further study.

**Scope and safety.** The lab uses toy curves only. It contains no discrete-logarithm solver and
does not touch Bitcoin wallets, networks, addresses or keys. The published secp256k1 domain
parameters appear once, as a big-integer arithmetic test vector in `tests/test_ec.py`.

## Layout

```
src/
  ec.py          curve arithmetic: field helpers, point validity, add/double/negate,
                 scalar multiplication (double-and-add + Montgomery ladder),
                 point counting (two independent methods), point orders, group structure
  curvegen.py    seeded curve generation (generic / j=0 / j=1728) and classification of
                 special properties (supersingular, anomalous, embedding degree, CM field,
                 automorphism group, ...)
  harness.py     seeding, timing, provenance and JSON output
tests/           54 unit tests (arithmetic, counting, group structure, generation, classification)
experiments/
  common.py                     the curve set (fixed by MASTER_SEED) and plotting style
  exp_a_geometry.py             Experiment A: geometry/statistics of x(kG), y(kG)
  exp_b_hypotheses.py           Experiment B: hypotheses H1–H5
  validate_curve_level_null.py  independent re-check of an apparent anomaly from H5
results/
  curves.json                   parameters, seeds, generators and properties of all curves
  exp_a/  summary.md|json, figures A1–A5
  exp_b/  summary.md|json, figures B1–B3, validate_curve_level_null.json
research_log.md   chronological log: hypotheses, methods, failures, conclusions (Phases 1–2)
findings.md       scientific interpretation
next_steps.md     candidate directions (Phase-2 direction 1 done; Phase-3 proposals)
findings_phase2.md  Phase-2 interpretation (endomorphisms; small-subgroup conjecture C3)
run_all.sh        reproduces Phase 1
run_phase2.sh     reproduces Phase 2
```

## Setup

Python 3.13 (any Python 3.10+ should work) with numpy, scipy, matplotlib and pytest:

```bash
pip install -r requirements.txt
```

The versions used for the recorded results were numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2
and pytest 9.1.1. They are also written into each `results/*/summary.json`.

## Reproducing

```bash
./run_all.sh                 # tests + Experiment A + Experiment B + validation, ~3–5 min, 1 core
python3 -m pytest -q tests   # tests only, ~1 s
```

Individual steps (run from `experiments/`):

```bash
python3 exp_a_geometry.py --perms 200
python3 exp_b_hypotheses.py
python3 validate_curve_level_null.py
```

All randomness is seeded: the curve search uses `MASTER_SEED = 20261025`, Experiment A uses seed 1,
Experiment B seed 2, and the validation script uses the seeds listed inside it. Re-running
produces byte-identical `summary.md` files; only timestamps and timings in the JSON change.

## The curve set

All curves are cyclic. G is a certified generator: its order is verified to equal #E exactly.

| name | family | seed | p | a | b | n = #E | trace t | n factored | G | \|Aut\| | CM field d |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GEN-A | generic | 20261025 | 128483 | 24512 | 66006 | 128747 | −263 | prime | (14206, 112558) | 2 | −444763 |
| GEN-B | generic | 20262025 | 99377 | 41730 | 46564 | 99529 | −151 | prime | (64712, 9780) | 2 | −374707 |
| GEN-C6 | generic | 20263040 | 123887 | 65122 | 5691 | 123846 | 42 | 2·3·20641 | (66324, 65101) | 2 | −123446 |
| J0-A | j = 0 | 20261025 | 118423 | 0 | 87528 | 119107 | −683 | prime | (75613, 77925) | 6 | −3 |
| J0-B | j = 0 | 20262025 | 114547 | 0 | 29121 | 115099 | −551 | prime | (10207, 102632) | 6 | −3 |
| J0-SS | j = 0, p ≡ 2 (mod 3) | 20261025 | 72269 | 0 | 44436 | 72270 | 0 | 2·3²·5·11·73 | (54804, 48486) | 2 | −72269 |
| J1728 | j = 1728, p ≡ 1 (mod 4) | 20261025 | 72269 | 44436 | 0 | 72746 | −476 | 2·36373 | (71233, 54818) | 4 | −1 |

Known special properties are flagged automatically (`results/curves.json`, field `notes`):

- J0-A and J0-B are ordinary j = 0 curves, the same family as secp256k1. They have an
  order-3 automorphism, which gives the GLV endomorphism.
- J0-SS is supersingular, with #E = p + 1 and embedding degree 2.
- J1728 has an order-4 automorphism, and (0, 0) is a 2-torsion point.

J0-SS and J1728 share p because the same seed was used with compatible congruence filters.

## Phase 2: blind discovery of small-degree endomorphisms

Additional code:
- `src/cm.py`: class-number-one orders and j-invariants, norm forms, eigenvalues of
  endomorphisms on a subgroup, and polynomial helpers over F_p.
- `experiments/exp_c_endomorphisms.py`: the batched exact search for multipliers m with
  x(mP) = R(x(P)).
- `experiments/exp_c2_threshold.py`: the small-subgroup boundary test.
- `experiments/exp_c3_small_subgroups.py`: exact d_min(m) for every m, and conjecture C3.
- `tests/test_cm.py`, `tests/test_phase2.py`.

```bash
./run_phase2.sh              # ~6–8 min on 1 core; writes results/phase2/reproduce/
```

The current evidence runs are `results/phase2/run-001` (endomorphism search), `run-005`
(threshold test), `run-006` (C3 fit) and `run-007` (C3 held-out). Runs 002–004 are earlier,
byte-identical versions; their manifests are stale because a script was edited afterwards
(see research_log.md P2.6). Each has its contract
(`results/phase2/contract-00N.json`) and a manifest with input hashes, argv, limits and output
hashes, made with the computation-audit bounded runner. There is no git repository, so each
manifest records the commit as "unavailable" and relies on input-file hashes instead.
