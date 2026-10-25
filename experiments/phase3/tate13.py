"""Exact arithmetic for the universal elliptic curve with a point of order 13 (Phase 3).

Tate normal form  E_{b,c}: y^2 + (1-c) x y - b y = x^3 - b x^2,  P = (0,0).
13P = O  <=>  x(6P) = x(7P); the numerator of x(6P) - x(7P) is one irreducible F(b,c) (computed).
All computations are in Q(b,c) with PARI (cypari2); reduction modulo F is exact division in b
(F is monic in b).
"""
from __future__ import annotations

import cypari2

pari = cypari2.Pari()
pari.allocatemem(4 * 10**9, silent=True) if "silent" in pari.allocatemem.__doc__ else pari.allocatemem(4 * 10**9)
b, c = pari("b"), pari("c")
A1, A2, A3, A4 = 1 - c, -b, -b, pari(0)


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 + y2 + A1 * x2 + A3 == 0:
            return None
        lam = (3 * x1**2 + 2 * A2 * x1 + A4 - A1 * y1) / (2 * y1 + A1 * x1 + A3)
    else:
        lam = (y2 - y1) / (x2 - x1)
    nu = y1 - lam * x1
    x3 = lam**2 + A1 * lam - A2 - x1 - x2
    return (x3, -(lam + A1) * x3 - nu - A3)


def multiples(n: int):
    P = (pari(0), pari(0))
    pts = [None, P]
    for _ in range(2, n + 1):
        pts.append(add(pts[-1], P))
    return pts


PTS = multiples(7)
X = {k: PTS[k][0] for k in range(1, 7)}                 # x(kP), k = 1..6 (x(-Q) = x(Q))
F = pari.numerator(PTS[6][0] - PTS[7][0])               # order-13 locus X_1(13)


def xk(k: int, q: int = 13):
    k %= q
    return X[min(k, q - k)]


def reduce_mod_F(poly):
    """Remainder of a polynomial in Q[b,c] on division by F (monic in b)."""
    return pari.lift(pari.Mod(poly, F))


def relation_matrix(m: int, d: int):
    rows = []
    for k in range(1, 7):
        Xv, Yv = xk(k), xk(m * k)
        rows.append([Xv**i for i in range(d + 1)] + [Xv**i * Yv for i in range(d + 1)])
    return pari.matrix(6, 2 * d + 2, [e for r in rows for e in r])
