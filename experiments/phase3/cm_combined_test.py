"""Test of the combined conjecture C2+C3 on CM curves with a small subgroup:

    d_min(m) = min( C3(m),  min{ N(alpha) : alpha in End(E), alpha acts as +-m on <P> } )

for E with j = j(O) (class number one), P of prime order q with q || #E(F_p) (then <P> is
End-stable and every alpha acts by a scalar; the eigenvalue comes from src/cm.py as in Phase 2).
For q inert in the CM field, q || #E is impossible (q | N(pi - 1) forces q^2 | #E); such curves are
tested with q^2 | #E as a control, where the prediction is C3 alone.

Usage: python experiments/phase3/cm_combined_test.py q [q ...]
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1]))

from exp_b_hypotheses import monomials, nullspace_mod_p  # noqa: E402
from exp_c3_small_subgroups import conjectured_dmin  # noqa: E402
from src.cm import CLASS_NUMBER_ONE, curve_from_j, eigenvalue_of_w, predicted_multipliers  # noqa: E402
from src.ec import Curve, is_probable_prime, legendre  # noqa: E402


def dmin_all(E, G, q):
    x = np.array([P[0] for P in E.multiples(G, q - 1)], dtype=np.int64)
    M = (q - 1) // 2
    d_int = M // 2
    half = np.arange(1, M + 1)
    out = {}
    for m in range(1, M + 1):
        for d in range(1, d_int + 1):
            if nullspace_mod_p(monomials(x[half - 1], x[(m * half) % q - 1], d, 1, E.p), E.p):
                out[m] = d
                break
    return out, d_int


def find_case(D, q, rng, want_exact):
    j = CLASS_NUMBER_ONE[D]
    for _ in range(200000):
        p = rng.randrange(3000, 400000)
        if not is_probable_prime(p):
            continue
        if D == -3:
            cands = [Curve(p, 0, rng.randrange(1, p)) for _ in range(6)]
        elif D == -4:
            cands = [Curve(p, rng.randrange(1, p), 0) for _ in range(4)]
        else:
            if j % p in (0, 1728 % p):
                continue
            nr = next(c for c in range(2, p) if legendre(c, p) == -1)
            cands = [curve_from_j(p, j), curve_from_j(p, j, twist=nr)]
        for E in cands:
            N = E.order()
            t = p + 1 - N
            if t % p == 0 or N % q:
                continue
            exact = N % (q * q) != 0
            if exact != want_exact:
                continue
            while True:
                G = E.mul(N // q if exact else N // (q * q), E.random_point(rng))
                if G is not None and (exact or E.mul(q, G) is None):
                    if E.mul(q, G) is None and G is not None:
                        break
            return E, G, t, exact
    return None


def main():
    rng = random.Random(20261030)
    for q in [int(v) for v in sys.argv[1:]]:
        M = (q - 1) // 2
        for D in sorted(CLASS_NUMBER_ONE, reverse=True)[:8]:
            split = legendre(D % q, q) == 1
            case = find_case(D, q, rng, want_exact=split)
            if case is None:
                print(f"q={q} D={D}: no suitable curve found")
                continue
            E, G, t, exact = case
            obs, d_int = dmin_all(E, G, q)
            c3 = {m: conjectured_dmin(m, q, d_int)[0] for m in range(1, M + 1)}
            pred = dict(c3)
            if exact:
                mu = eigenvalue_of_w(E.p, t, D, q)
                for d in range(1, d_int + 1):
                    for mm in predicted_multipliers(D, d, mu, q):
                        mm = min(mm, q - mm)
                        if mm and d < pred[mm]:
                            pred[mm] = d
            bad = {m: (obs[m], pred[m]) for m in obs if obs[m] != pred[m]}
            lower = {m: obs[m] for m in obs if obs[m] < c3[m]}
            print(f"q={q} D={D} ({'split' if split else 'inert/ramified'}; q{'||' if exact else '^2|'}#E) p={E.p}: "
                  f"combined prediction matches {M - len(bad)}/{M}; below-C3 multipliers {lower}; mismatches (obs, pred) {bad}",
                  flush=True)


if __name__ == "__main__":
    main()
