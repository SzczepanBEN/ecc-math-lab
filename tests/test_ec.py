"""Verification of the field / group arithmetic in src/ec.py."""
import itertools
import math
import random

import pytest

from src.ec import (INF, Curve, cube_roots_of_unity, factorize, find_generator,
                    is_probable_prime, legendre, sqrt_mod)

SMALL_CURVES = [
    Curve(17, 2, 2),        # textbook, #E = 19
    Curve(23, 1, 1),        # textbook, #E = 28
    Curve(97, 2, 3),
    Curve(101, 0, 7),       # j = 0, p = 2 mod 3 (supersingular)
    Curve(103, 0, 5),       # j = 0, p = 1 mod 3
    Curve(107, 3, 0),       # j = 1728, p = 3 mod 4 (supersingular)
    Curve(109, 5, 0),       # j = 1728, p = 1 mod 4
    Curve(131, 17, 42),
]


def brute_points(E):
    p = E.p
    return [INF] + [(x, y) for x in range(p) for y in range(p) if (y * y - E.rhs(x)) % p == 0]


# ---------------- number theory ----------------

def test_primality_matches_trial_division():
    for n in range(-5, 5000):
        trial = n >= 2 and all(n % d for d in range(2, math.isqrt(n) + 1))
        assert is_probable_prime(n) == trial, n


def test_factorize_roundtrip():
    rng = random.Random(1)
    for _ in range(200):
        n = rng.randrange(2, 10**9)
        assert math.prod(q**e for q, e in factorize(n).items()) == n
        assert all(is_probable_prime(q) for q in factorize(n))


@pytest.mark.parametrize("p", [7, 13, 17, 97, 101, 65537, 1000003, 2**61 - 1])
def test_sqrt_mod(p):
    rng = random.Random(p)
    for _ in range(100):
        a = rng.randrange(p)
        r = sqrt_mod(a, p)
        if legendre(a, p) == -1:
            assert r is None
        else:
            assert r * r % p == a


def test_cube_roots_of_unity():
    for p in [7, 13, 103, 82009]:
        roots = cube_roots_of_unity(p)
        assert len(roots) == 2 and all(pow(r, 3, p) == 1 and r != 1 for r in roots)
    assert cube_roots_of_unity(101) == []


# ---------------- curve construction ----------------

def test_rejects_singular_and_bad_p():
    with pytest.raises(ValueError):
        Curve(17, 0, 0)
    with pytest.raises(ValueError):
        Curve(23, 3 * 23 - 3, 2)  # a=-3, b=2: 4(-27)+27*4 = 0
    with pytest.raises(ValueError):
        Curve(15, 1, 1)


def test_j_invariant_special_values():
    assert Curve(103, 0, 5).j_invariant == 0
    assert Curve(109, 5, 0).j_invariant == 1728 % 109


# ---------------- point counting ----------------

@pytest.mark.parametrize("E", SMALL_CURVES, ids=lambda E: f"p{E.p}a{E.a}b{E.b}")
def test_point_count_three_ways(E):
    brute = brute_points(E)
    assert len(brute) == E.order() == E.order_legendre()
    assert sorted(brute[1:]) == sorted(E.all_points()[1:])
    assert all(E.is_on_curve(P) for P in brute)
    # Hasse bound
    assert abs(E.p + 1 - E.order()) <= 2 * math.sqrt(E.p)


def test_textbook_vectors():
    # Paar & Pelzl, "Understanding Cryptography", Ex. 9.x: y^2 = x^3 + 2x + 2 over F_17
    E = Curve(17, 2, 2)
    G = (5, 1)
    assert E.order() == 19 and E.double(G) == (6, 3) and E.mul(19, G) is INF
    assert E.multiples(G, 19)[-1] is INF
    # y^2 = x^3 + x + 1 over F_23 has 28 points (standard textbook example)
    assert Curve(23, 1, 1).order() == 28


def test_supersingular_families():
    for p in [101, 107, 131, 1013, 120863]:
        if p % 3 == 2:
            assert Curve(p, 0, 7).order() == p + 1
        if p % 4 == 3:
            assert Curve(p, 3, 0).order() == p + 1


# ---------------- group law ----------------

@pytest.mark.parametrize("E", SMALL_CURVES, ids=lambda E: f"p{E.p}a{E.a}b{E.b}")
def test_group_axioms_exhaustive_small(E):
    pts = brute_points(E)
    rng = random.Random(E.p)
    for P in pts:
        assert E.add(P, INF) == P == E.add(INF, P)
        assert E.add(P, E.neg(P)) is INF
        assert E.is_on_curve(E.double(P))
    for P, Q in itertools.product(pts, repeat=2):
        R = E.add(P, Q)
        assert E.is_on_curve(R)
        assert R == E.add(Q, P)
    for _ in range(3000):
        P, Q, R = (rng.choice(pts) for _ in range(3))
        assert E.add(E.add(P, Q), R) == E.add(P, E.add(Q, R))


def test_associativity_random_large():
    E = Curve(1000003, 12345, 6789)
    rng = random.Random(7)
    for _ in range(2000):
        P, Q, R = (E.random_point(rng) for _ in range(3))
        assert E.add(E.add(P, Q), R) == E.add(P, E.add(Q, R))
    # edge cases: doubling via add, P + P + (-P)
    P = E.random_point(rng)
    assert E.add(E.add(P, P), E.neg(P)) == P


def test_two_torsion():
    E = Curve(97, 2, 3)
    for P in brute_points(E)[1:]:
        if P[1] == 0:
            assert E.double(P) is INF and E.point_order(P) == 2


# ---------------- scalar multiplication ----------------

@pytest.mark.parametrize("E", SMALL_CURVES[:4], ids=lambda E: f"p{E.p}")
def test_scalar_mul_matches_repeated_addition(E):
    rng = random.Random(3)
    N = E.order()
    for _ in range(5):
        P = E.random_point(rng)
        seq = E.multiples(P, 2 * N + 3)
        for k in range(1, 2 * N + 4):
            assert E.mul(k, P) == seq[k - 1] == E.mul_ladder(k, P)
        assert E.mul(0, P) is INF and E.mul(-1, P) == E.neg(P)


def test_scalar_mul_homomorphism_large():
    E = Curve(1000003, 12345, 6789)
    rng = random.Random(11)
    N = E.order()
    for _ in range(300):
        P = E.random_point(rng)
        a, b = rng.randrange(-N, N), rng.randrange(-N, N)
        assert E.add(E.mul(a, P), E.mul(b, P)) == E.mul(a + b, P)
        assert E.mul(a, E.mul(b, P)) == E.mul(a * b, P)
        assert E.mul(a, P) == E.mul_ladder(a, P)
        assert E.mul(N, P) is INF  # Lagrange


# ---------------- orders and group structure ----------------

@pytest.mark.parametrize("E", SMALL_CURVES, ids=lambda E: f"p{E.p}a{E.a}b{E.b}")
def test_group_structure_vs_bruteforce(E):
    pts = brute_points(E)
    N = len(pts)
    orders = []
    for P in pts:
        # brute-force order by repeated addition
        k, Q = 1, P
        while Q is not INF:
            Q, k = E.add(Q, P), k + 1
        assert E.point_order(P) == k and N % k == 0
        orders.append(k)
    n1 = max(orders)  # exponent of a finite abelian group = max element order
    n2 = N // n1
    assert E.group_structure(random.Random(0)) == (n1, n2)
    assert n1 % n2 == 0 and (E.p - 1) % n2 == 0
    # Z/n1 x Z/n2 has exactly n2^2 points of order dividing n2 -- check structurally
    assert sum(1 for k in orders if n2 % k == 0) == n2 * n2
    if n2 == 1:
        G = find_generator(E, random.Random(0))
        assert len(set(E.multiples(G, N))) == N


# ---------------- public test vector: secp256k1 ----------------
# Used only as a big-integer arithmetic test vector (published SEC 2 domain parameters).
SECP_P = 2**256 - 2**32 - 977
SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_G = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
          0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)


def test_secp256k1_parameters_as_arithmetic_vector():
    E = Curve(SECP_P, 0, 7)
    assert E.is_on_curve(SECP_G) and E.j_invariant == 0 and SECP_P % 3 == 1
    assert is_probable_prime(SECP_N)
    assert E.mul(SECP_N, SECP_G) is INF
    assert E.mul(SECP_N + 5, SECP_G) == E.mul_ladder(5, SECP_G)
    # GLV endomorphism: some cube root of unity beta mod p and lambda mod n with
    # lambda*G = (beta*x_G, y_G).
    betas = cube_roots_of_unity(SECP_P)
    lams = cube_roots_of_unity(SECP_N)
    hits = [(bt, lm) for bt in betas for lm in lams
            if E.mul(lm, SECP_G) == (bt * SECP_G[0] % SECP_P, SECP_G[1])]
    assert len(hits) == 2  # (beta, lambda) and (beta^2, lambda^2)
