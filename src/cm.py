"""Complex-multiplication helpers for Phase 2 (toy CM curves over F_p).

Conventions
-----------
* An imaginary quadratic order O is given by its discriminant D < 0 (D = 0, 1 mod 4), with
  Z-basis {1, w}, w = (D + sqrt(D))/2, so that  N(x + y w) = x^2 + D x y + (D^2 - D)/4 y^2.
* For an ordinary curve E/F_p with trace t, t^2 - 4p = v^2 D_O whenever Z[pi] is contained in
  O = End(E); then pi = (t - v D_O)/2 + v w.
* On a cyclic subgroup <G'> of prime order q fixed by Frobenius, every alpha in O acts as an
  integer multiplier m_alpha mod q ("eigenvalue"); pi acts as 1.

Everything is exact integer / F_p arithmetic. Polynomials over F_p are lists of ints, lowest
degree first, without trailing zeros (the zero polynomial is []).
"""
from __future__ import annotations

import math

from .ec import Curve, factorize, sqrt_mod

# The 13 imaginary quadratic orders of class number one and their j-invariants.
# Values cross-checked against published tables (see research_log.md, Phase 2, literature check);
# the experiment additionally checks each one computationally (CM field of t^2 - 4p).
CLASS_NUMBER_ONE = {
    -3: 0,
    -4: 1728,
    -7: -3375,
    -8: 8000,
    -11: -32768,
    -12: 54000,
    -16: 287496,
    -19: -884736,
    -27: -12288000,
    -28: 16581375,
    -43: -884736000,
    -67: -147197952000,
    -163: -262537412640768000,
}


# ----------------------------------------------------------------- quadratic orders

def squarefree_part(n: int) -> int:
    s = -1 if n < 0 else 1
    for q, e in factorize(abs(n)).items():
        if e % 2:
            s *= q
    return s


def fundamental_discriminant(delta: int) -> tuple[int, int]:
    """Write delta < 0 as f^2 * D_K with D_K fundamental; return (D_K, f)."""
    if delta >= 0:
        raise ValueError("expected a negative discriminant")
    s = squarefree_part(delta)
    dk = s if s % 4 == 1 else 4 * s
    f2, r = divmod(delta, dk)
    f = math.isqrt(f2)
    if r or f * f != f2:
        raise ValueError(f"{delta} is not of the form f^2 * D_K")
    return dk, f


def norm_form(D: int):
    c = (D * D - D) // 4
    return lambda x, y: x * x + D * x * y + c * y * y


def elements_of_norm(D: int, d: int) -> list[tuple[int, int]]:
    """All (x, y) with N(x + y w) = d in the order of discriminant D (exhaustive)."""
    N = norm_form(D)
    ymax = math.isqrt(4 * d // abs(D)) + 1     # N >= |D| y^2 / 4
    out = []
    for y in range(-ymax, ymax + 1):
        # N = (x + yD/2)^2 + |D| y^2 / 4, so |x + y D / 2| <= sqrt(d)
        centre = -y * D // 2
        for x in range(centre - math.isqrt(d) - 2, centre + math.isqrt(d) + 3):
            if N(x, y) == d:
                out.append((x, y))
    return out


def r_O(D: int, d: int) -> int:
    """Number of elements of norm d in the order of discriminant D."""
    return len(elements_of_norm(D, d))


# ----------------------------------------------------------------- curves with given j

def curve_from_j(p: int, j: int, twist: int = 1) -> Curve:
    """y^2 = x^3 + a x + b with j(E) = j (j != 0, 1728 mod p), optionally twisted by c."""
    j %= p
    if j in (0, 1728 % p):
        raise ValueError("use the j = 0 / 1728 families directly")
    k = j * (1728 - j) % p
    a, b = 3 * k % p, 2 * k * (1728 - j) % p
    return Curve(p, a * twist * twist % p, b * pow(twist, 3, p) % p)


# ----------------------------------------------------------------- eigenvalues on <G'>

def frobenius_coordinates(p: int, t: int, D: int) -> tuple[int, int]:
    """(u, v) with pi = u + v w in the order of discriminant D, or raise if Z[pi] is not in it."""
    disc = t * t - 4 * p
    v2, r = divmod(disc, D)
    v = math.isqrt(v2) if v2 >= 0 else -1
    if r or v < 0 or v * v != v2:
        raise ValueError(f"t^2 - 4p = {disc} is not v^2 * {D}")
    u2 = t - v * D
    if u2 % 2:
        raise ValueError("pi not integral in this order")
    return u2 // 2, v


def eigenvalue_of_w(p: int, t: int, D: int, q: int) -> int:
    """The multiplier mu by which w acts on a Frobenius-fixed subgroup of prime order q.

    pi acts as 1, so u + v mu = 1 (mod q). The other root of w's minimal polynomial corresponds
    to the conjugate embedding; the predicted multiplier SETS below do not depend on the choice.
    """
    u, v = frobenius_coordinates(p, t, D)
    if v % q == 0:
        raise ValueError("q divides v: eigenvalue not determined by Frobenius")
    mu = (1 - u) * pow(v, -1, q) % q
    c = (D * D - D) // 4
    assert (mu * mu - D * mu + c) % q == 0, "mu is not a root of w's minimal polynomial"
    return mu


def predicted_multipliers(D: int, d: int, mu: int, q: int) -> set[int]:
    """S_d = { x + y mu mod q : N(x + y w) = d }."""
    return {(x + y * mu) % q for x, y in elements_of_norm(D, d)}


def trace_for_multiplier(m: int, d: int, q: int) -> int | None:
    """The unique integer T, |T| <= 2 sqrt(d), with m^2 - T m + d = 0 mod q (None if absent)."""
    tmax = math.isqrt(4 * d)
    hits = [T for T in range(-tmax, tmax + 1) if (m * m - T * m + d) % q == 0]
    assert len(hits) <= 1
    return hits[0] if hits else None


# ----------------------------------------------------------------- polynomials over F_p

def ptrim(a: list[int]) -> list[int]:
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def pdeg(a: list[int]) -> int:
    return len(ptrim(a)) - 1


def pmonic(a: list[int], p: int) -> list[int]:
    a = ptrim(a)
    if not a:
        return a
    inv = pow(a[-1], -1, p)
    return [c * inv % p for c in a]


def pmul(a, b, p):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return ptrim(out)


def pdivmod(a, b, p):
    a, b = ptrim([c % p for c in a]), ptrim([c % p for c in b])
    if not b:
        raise ZeroDivisionError
    qt = [0] * max(len(a) - len(b) + 1, 1)
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % p
        s = len(a) - len(b)
        qt[s] = c
        for i, y in enumerate(b):
            a[s + i] = (a[s + i] - c * y) % p
        a = ptrim(a)
    return ptrim(qt), a


def pgcd(a, b, p):
    a, b = ptrim([c % p for c in a]), ptrim([c % p for c in b])
    while b:
        a, b = b, pdivmod(a, b, p)[1]
    return pmonic(a, p)


def pderiv(a, p):
    return ptrim([i * c % p for i, c in enumerate(a)][1:])


def kernel_structure(den: list[int], f: list[int], p: int) -> dict:
    """Split a denominator as  den = c * s * h^2  (h = gcd(den, den'), s = den / h^2) and test
    the Velu-type prediction: s | f (2-torsion kernel points), h squarefree, s, h coprime."""
    den = pmonic(den, p)
    h = pgcd(den, pderiv(den, p), p) if pdeg(den) > 0 else [1]
    s, r = pdivmod(den, pmul(h, h, p), p)
    return {
        "deg_den": pdeg(den), "deg_h": pdeg(h), "deg_s": pdeg(s) if s else -1,
        "exact_split": not r,
        "s_divides_f": bool(s) and not pdivmod(f, s, p)[1],
        "h_squarefree": pdeg(h) <= 0 or pdeg(pgcd(h, pderiv(h, p), p)) == 0,
        "s_h_coprime": pdeg(h) <= 0 or pdeg(pgcd(s, h, p)) == 0,
    }
