"""Tests for the Phase-2 search machinery (experiments/exp_c_endomorphisms.py)."""
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))

from exp_b_hypotheses import monomials, nullspace_mod_p  # noqa: E402
from exp_c_endomorphisms import (batch_rank_deficient, doubling_formula, modinv_vec,  # noqa: E402
                                 reference_screen, screen_degree)
from exp_c3_small_subgroups import conjectured_dmin  # noqa: E402
from src.ec import Curve  # noqa: E402


def test_modinv_vec():
    p = 10007
    a = np.arange(1, 2000, dtype=np.int64)
    assert np.all(a * modinv_vec(a, p) % p == 1)


def test_batch_rank_deficient_matches_reference():
    rng = np.random.default_rng(0)
    p = 8191
    for r, c in [(3, 2), (5, 4), (9, 8), (12, 11)]:
        B = rng.integers(0, p, size=(400, r, c))
        # plant rank deficiency in a quarter of the stack: last column = combination of others
        for i in range(0, 400, 4):
            coef = rng.integers(0, p, size=c - 1)
            B[i, :, -1] = (B[i, :, :-1] @ coef) % p
        # and some zero columns / repeated rows to exercise pivoting
        B[1, :, 0] = 0
        B[2, 1] = B[2, 0]
        got = batch_rank_deficient(B, p)
        ref = np.array([len(nullspace_mod_p(B[i], p)) > 0 for i in range(len(B))])
        assert np.array_equal(got, ref)


def test_screen_has_no_false_negatives_on_planted_relation():
    # On any curve, m = 2 has a degree-4 relation (doubling); the screen must report it at d = 4
    E = Curve(10007, 3, 7)
    rng = random.Random(1)
    N = E.order()
    G = E.random_point(rng)
    from src.ec import factorize
    q = max(factorize(N))
    G = E.mul(N // q, G)
    while G is None:
        G = E.mul(N // q, E.random_point(rng))
    x = np.array([P[0] for P in E.multiples(G, q - 1)], dtype=np.int64)
    ms = np.arange(1, q, dtype=np.int64)
    hits4 = set(screen_degree(x, q, E.p, 4, ms).tolist())
    assert {1, 2, q - 1, q - 2} <= hits4
    assert 2 not in set(screen_degree(x, q, E.p, 3, ms).tolist())
    for m in [1, 2, 3, 5, q - 2]:
        assert (m in hits4) == reference_screen(x, q, E.p, 4, m)
    # the closed-form doubling map satisfies the relation on every point
    num, den = doubling_formula(E.a, E.b, E.p)
    p = E.p
    for k in range(1, 30):
        X = int(x[k - 1])
        Y = int(x[(2 * k) % q - 1])
        nv = sum(c * pow(X, i, p) for i, c in enumerate(num)) % p
        dv = sum(c * pow(X, i, p) for i, c in enumerate(den)) % p
        assert (Y * dv - nv) % p == 0


def test_conjectured_dmin_integer_multipliers():
    # b = 1 gives d = m^2 (integer multipliers) whenever that is below threshold
    q = 157
    d_int = next(d for d in range(1, q) if 2 * d + 2 > (q - 1) // 2)
    for m in range(1, 7):
        d, rep = conjectured_dmin(m, q, d_int)
        assert d == min(m * m, d_int) or d < m * m
    assert conjectured_dmin(52, q, d_int)[0] == 26     # m = -1/3 mod 157: (M + 1 - 1)/3 = 26
