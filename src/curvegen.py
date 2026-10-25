"""Reproducible generation and classification of toy elliptic curves.

Every generated curve carries the seed and search parameters that produced it, and a
classification of the standard special properties so that "structure" which follows from
them is attributed to known theory rather than reported as a discovery.
"""
from __future__ import annotations

import math
import random
from dataclasses import asdict, dataclass

from .ec import Curve, cube_roots_of_unity, factorize, find_generator, is_probable_prime


@dataclass
class CurveRecord:
    name: str
    family: str          # "generic" | "j0" | "j1728"
    p: int
    a: int
    b: int
    seed: int
    order: int
    trace: int
    structure: tuple     # (n1, n2): E = Z/n1 x Z/n2
    generator: tuple | None
    properties: dict

    def curve(self) -> Curve:
        return Curve(self.p, self.a, self.b, self.name)

    def to_json(self) -> dict:
        d = asdict(self)
        d["structure"] = list(self.structure)
        d["generator"] = list(self.generator) if self.generator else None
        return d


def random_prime(rng: random.Random, lo: int, hi: int, mod: tuple[int, int] | None = None) -> int:
    """Uniform-ish random prime in [lo, hi), optionally with p = r (mod m) for mod=(r, m)."""
    while True:
        p = rng.randrange(lo, hi)
        if mod and p % mod[1] != mod[0]:
            continue
        if is_probable_prime(p):
            return p


def embedding_degree(p: int, r: int, kmax: int = 50) -> int | None:
    """Smallest k <= kmax with r | p^k - 1 (MOV/Frey-Rueck degree), None if > kmax."""
    if r == p:
        return None
    x = 1
    for k in range(1, kmax + 1):
        x = x * p % r
        if x == 1:
            return k
    return None


def squarefree_part(n: int) -> int:
    sign = -1 if n < 0 else 1
    s = 1
    for q, e in factorize(abs(n)).items():
        if e % 2:
            s *= q
    return sign * s


def classify(E: Curve) -> dict:
    """Known special properties of E(F_p) relevant to interpreting experiments."""
    p, a, b = E.p, E.a, E.b
    N = E.order()
    t = p + 1 - N
    fac = factorize(N)
    r = max(fac)
    disc_frob = t * t - 4 * p  # discriminant of Frobenius x^2 - t x + p; < 0 by Hasse
    props = {
        "j_invariant": E.j_invariant,
        "j0": a == 0,
        "j1728": b == 0,
        "hasse_ok": abs(t) <= 2 * math.isqrt(p) + 1 and t * t <= 4 * p,
        "supersingular": t % p == 0,
        "anomalous": N == p,
        "prime_order": len(fac) == 1 and fac[r] == 1,
        "order_factorization": {str(k): v for k, v in sorted(fac.items())},
        "largest_prime_factor": r,
        "cofactor": N // r,
        "embedding_degree_of_r": embedding_degree(p, r),
        "twist_order": p + 1 + t,
        "frobenius_discriminant": disc_frob,
        "cm_field_squarefree_d": squarefree_part(disc_frob),
        "p_mod_3": p % 3,
        "p_mod_4": p % 4,
    }
    # Automorphisms beyond +-1 exist only for j = 0 (order 6 when p = 1 mod 3)
    # and j = 1728 (order 4 when p = 1 mod 4).
    if a == 0:
        props["aut_group_order"] = 6 if p % 3 == 1 else 2
        props["cube_roots_of_unity"] = cube_roots_of_unity(p)
    elif b == 0:
        props["aut_group_order"] = 4 if p % 4 == 1 else 2
    else:
        props["aut_group_order"] = 2
    notes = []
    if props["supersingular"]:
        notes.append("supersingular (t=0): #E=p+1, embedding degree <= 2 (MOV/Frey-Rueck)")
    if props["anomalous"]:
        notes.append("anomalous (#E=p): Smart/Satoh-Araki/Semaev p-adic lift applies")
    if a == 0 and p % 3 == 2:
        notes.append("j=0 with p=2 mod 3: x->x^3 is a bijection of F_p, hence supersingular")
    if b == 0 and p % 4 == 3:
        notes.append("j=1728 with p=3 mod 4: supersingular")
    if a == 0 and p % 3 == 1:
        notes.append("j=0 ordinary: order-3 automorphism (x,y)->(beta x, y), GLV endomorphism")
    if b == 0 and p % 4 == 1:
        notes.append("j=1728 ordinary: order-4 automorphism (x,y)->(-x, i y); (0,0) is 2-torsion")
    k = props["embedding_degree_of_r"]
    if k is not None and k <= 6 and not props["supersingular"]:
        notes.append(f"small embedding degree k={k}")
    props["notes"] = notes
    return props


def make_record(E: Curve, family: str, seed: int, rng: random.Random) -> CurveRecord:
    n1, n2 = E.group_structure(rng)
    G = find_generator(E, rng) if n2 == 1 else None
    N = E.order()
    return CurveRecord(E.name, family, E.p, E.a, E.b, seed, N, E.p + 1 - N, (n1, n2), G, classify(E))


def generate_curve(seed: int, bits: int, family: str = "generic", require_prime_order: bool = False,
                   p_mod: tuple[int, int] | None = None, name: str = "",
                   max_tries: int = 100_000) -> CurveRecord:
    """Search for a curve from a seeded RNG.

    family:  "generic" -- random a, b with a != 0, b != 0 (so j != 0, 1728)
             "j0"      -- y^2 = x^3 + b
             "j1728"   -- y^2 = x^3 + a x
    """
    rng = random.Random(seed)
    lo, hi = 1 << (bits - 1), 1 << bits
    for _ in range(max_tries):
        p = random_prime(rng, lo, hi, p_mod)
        a = rng.randrange(1, p) if family in ("generic", "j1728") else 0
        b = rng.randrange(1, p) if family in ("generic", "j0") else 0
        try:
            E = Curve(p, a, b, name or f"{family}-{bits}b-s{seed}")
        except ValueError:
            continue
        if family == "generic" and E.j_invariant in (0, 1728 % p):
            continue
        if require_prime_order and not is_probable_prime(E.order()):
            continue
        return make_record(E, family, seed, rng)
    raise RuntimeError("no curve found")
