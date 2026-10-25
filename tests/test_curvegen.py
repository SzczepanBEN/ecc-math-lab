"""Curve generation is reproducible and its classifications are correct."""
import random

from src.curvegen import classify, embedding_degree, generate_curve
from src.ec import Curve, is_probable_prime


def test_seed_reproducibility():
    a = generate_curve(5, 14, "generic", require_prime_order=True)
    b = generate_curve(5, 14, "generic", require_prime_order=True)
    assert a.to_json() == b.to_json()
    assert generate_curve(6, 14).to_json() != a.to_json()


def test_generic_curves_are_valid_and_classified():
    for seed in range(15):
        r = generate_curve(seed, 13, "generic")
        E = r.curve()
        assert E.j_invariant not in (0, 1728 % E.p)
        assert r.order == E.order_legendre()
        assert r.properties["hasse_ok"]
        assert r.properties["aut_group_order"] == 2
        n1, n2 = r.structure
        assert n1 * n2 == r.order and (E.p - 1) % n2 == 0


def test_prime_order_and_generator():
    for seed in range(5):
        r = generate_curve(seed, 14, "generic", require_prime_order=True)
        assert is_probable_prime(r.order) and r.properties["prime_order"]
        E = r.curve()
        assert E.is_on_curve(r.generator) and E.mul(r.order, r.generator) is None


def test_j0_classification():
    r = generate_curve(1, 14, "j0", p_mod=(2, 3))
    assert r.properties["supersingular"] and r.order == r.p + 1
    assert r.properties["embedding_degree_of_r"] in (1, 2)
    r = generate_curve(1, 14, "j0", p_mod=(1, 3), require_prime_order=True)
    assert not r.properties["supersingular"] and r.properties["aut_group_order"] == 6
    # for j = 0 ordinary curves the CM field is Q(sqrt(-3)): t^2 - 4p = -3 v^2
    assert r.properties["cm_field_squarefree_d"] == -3


def test_j1728_classification():
    r = generate_curve(2, 14, "j1728", p_mod=(3, 4))
    assert r.properties["supersingular"]
    r = generate_curve(2, 14, "j1728", p_mod=(1, 4))
    assert r.properties["aut_group_order"] == 4 and r.properties["cm_field_squarefree_d"] == -1


def test_embedding_degree():
    assert embedding_degree(11, 3) == 2      # 11 = -1 mod 3
    assert embedding_degree(13, 3) == 1
    assert embedding_degree(101, 7) is not None


def test_anomalous_detection():
    # search a tiny anomalous curve (#E = p) by brute force and check the flag
    rng = random.Random(0)
    for _ in range(20000):
        p = rng.choice([53, 59, 61, 67, 71, 73, 79, 83, 89, 97])
        try:
            E = Curve(p, rng.randrange(p), rng.randrange(p))
        except ValueError:
            continue
        if E.order() == p:
            assert classify(E)["anomalous"]
            return
    raise AssertionError("no anomalous curve found")
