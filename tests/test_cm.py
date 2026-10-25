"""Tests for the CM helpers in src/cm.py (hand-checkable cases)."""
import random

import pytest

from src.cm import (CLASS_NUMBER_ONE, curve_from_j, eigenvalue_of_w, elements_of_norm,
                    fundamental_discriminant, frobenius_coordinates, kernel_structure, pdivmod,
                    pgcd, pmul, predicted_multipliers, r_O, trace_for_multiplier)
from src.ec import Curve, is_probable_prime


@pytest.mark.parametrize("D,d,expected", [
    (-4, 1, 4), (-4, 2, 4), (-4, 3, 0), (-4, 5, 8), (-4, 4, 4),       # Z[i]
    (-3, 1, 6), (-3, 2, 0), (-3, 3, 6), (-3, 4, 6), (-3, 7, 12),      # Z[zeta_3]
    (-12, 1, 2), (-12, 3, 2), (-12, 4, 6), (-12, 7, 4),               # Z[sqrt(-3)]
    (-16, 1, 2), (-16, 4, 4), (-16, 5, 4),                            # Z[2i]
    (-8, 1, 2), (-8, 2, 2), (-8, 3, 4),                               # Z[sqrt(-2)]
    (-7, 1, 2), (-7, 2, 4), (-7, 4, 6),                               # Z[(1+sqrt(-7))/2]
    (-163, 1, 2), (-163, 40, 0), (-163, 41, 4),
])
def test_representation_counts(D, d, expected):
    assert r_O(D, d) == expected


def test_elements_of_norm_bruteforce():
    # compare with a wide brute-force box for several orders
    for D in (-3, -4, -7, -11, -12, -27, -28):
        c = (D * D - D) // 4
        for d in range(1, 30):
            box = {(x, y) for x in range(-60, 61) for y in range(-20, 21)
                   if x * x + D * x * y + c * y * y == d}
            assert box == set(elements_of_norm(D, d)), (D, d)


@pytest.mark.parametrize("delta,expected", [
    (-3, (-3, 1)), (-12, (-3, 2)), (-27, (-3, 3)), (-4, (-4, 1)), (-16, (-4, 2)),
    (-7, (-7, 1)), (-28, (-7, 2)), (-8, (-8, 1)), (-32, (-8, 2)), (-63, (-7, 3)), (-20, (-20, 1)),
])
def test_fundamental_discriminant(delta, expected):
    assert fundamental_discriminant(delta) == expected


def test_curve_from_j():
    for p in (10007, 30011):
        for D, j in CLASS_NUMBER_ONE.items():
            if D in (-3, -4):
                continue
            E = curve_from_j(p, j)
            assert E.j_invariant == j % p
            assert curve_from_j(p, j, twist=5).j_invariant == j % p


def test_eigenvalue_on_known_curve():
    # j = 1728: y^2 = x^3 + x over p = 1 mod 4; End = Z[i]; i acts as (x, y) -> (-x, i y)
    rng = random.Random(1)
    p = 10009
    E = Curve(p, 1, 0)
    N = E.order()
    t = p + 1 - N
    q = max(int(k) for k in __import__("src.ec", fromlist=["factorize"]).factorize(N))
    assert is_probable_prime(q) and q > 100
    G = E.mul(N // q, E.random_point(rng))
    while G is None:
        G = E.mul(N // q, E.random_point(rng))
    mu = eigenvalue_of_w(p, t, -4, q)  # w = (-4 + sqrt(-4))/2 = -2 + i, so i = w + 2
    lam = (mu + 2) % q
    i = __import__("src.ec", fromlist=["sqrt_mod"]).sqrt_mod(p - 1, p)
    P = E.mul(lam, G)
    assert P in ((-G[0] % p, i * G[1] % p), (-G[0] % p, (-i) * G[1] % p))
    # the predicted set at d = 1 is {+-1, +-lambda}
    assert predicted_multipliers(-4, 1, mu, q) == {1, q - 1, lam, q - lam}


def test_frobenius_coordinates_rejects_wrong_order():
    p = 10009
    t = p + 1 - Curve(p, 1, 0).order()
    frobenius_coordinates(p, t, -4)
    with pytest.raises(ValueError):
        frobenius_coordinates(p, t, -7)


def test_trace_for_multiplier():
    q = 101
    # m = 2 has m^2 - 4m + 4 = 0 -> T = 4 for d = 4
    assert trace_for_multiplier(2, 4, q) == 4
    assert trace_for_multiplier(q - 1, 1, q) == -2


def test_polynomial_helpers_and_kernel_structure():
    p = 10007
    rng = random.Random(3)
    for _ in range(50):
        a = [rng.randrange(p) for _ in range(rng.randrange(1, 8))] + [1]
        b = [rng.randrange(p) for _ in range(rng.randrange(1, 6))] + [1]
        qq, r = pdivmod(pmul(a, b, p), b, p)
        assert qq == a and r == []
        g = pgcd(pmul(a, b, p), b, p)
        assert g == b
    f = [3, 2, 0, 1]                       # x^3 + 2x + 3
    e = [5, 1]                             # x + 5
    h = pmul([7, 1], [11, 1], p)           # (x + 7)(x + 11)
    den = pmul(e, pmul(h, h, p), p)
    ks = kernel_structure(den, pmul(e, [1, 0, 1], p), p)
    assert ks["exact_split"] and ks["deg_h"] == 2 and ks["deg_s"] == 1 and ks["s_divides_f"]
    assert not kernel_structure(den, f, p)["s_divides_f"]
