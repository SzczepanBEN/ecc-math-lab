"""Universal elliptic curve with a point of prime order q (Tate normal form), for any q >= 5.

E_{b,c}: y^2 + (1-c)xy - by = x^3 - bx^2, P = (0,0). With h = (q-1)/2, qP = O <=> x(hP) = x((h+1)P);
the order-q locus F_q is the factor of the numerator of x(hP) - x((h+1)P) that the caller selects
(for q = 13 the numerator is already irreducible). Generalises tate13.py.
"""
from __future__ import annotations

import cypari2

pari = cypari2.Pari()
try:
    pari.allocatemem(6 * 10**9)
except Exception:  # pragma: no cover
    pass
b, c = pari("b"), pari("c")
A1, A2, A3, A4 = 1 - c, -b, -b, pari(0)
a1, a2, a3 = A1, A2, A3
b2 = a1**2 + 4 * a2
b4 = a1 * a3
b6 = a3**2
b8 = a2 * a3**2
DELTA = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
C4 = b2**2 - 24 * b4


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


class TateModel:
    def __init__(self, q: int):
        self.q = q
        self.h = (q - 1) // 2
        P = (pari(0), pari(0))
        pts = [None, P]
        for _ in range(2, self.h + 2):
            pts.append(add(pts[-1], P))
        self.pts = pts
        self.X = {k: pts[k][0] for k in range(1, self.h + 1)}
        self.raw = pari.numerator(pts[self.h][0] - pts[self.h + 1][0])
        fa = pari.factor(self.raw)
        self.raw_factors = [(fa[0][i], int(fa[1][i])) for i in range(len(fa[0]))]
        # the order-q component: the factor of largest total degree (checked over F_p by callers)
        self.F = max((f for f, _ in self.raw_factors),
                     key=lambda f: int(pari.poldegree(f, "b")) + int(pari.poldegree(f, "c")))

    def xk(self, k: int):
        k %= self.q
        return self.X[min(k, self.q - k)]

    def zero_in_K(self, expr) -> bool:
        """expr == 0 in the function field Q(b,c)/(F)."""
        num = pari.numerator(expr)
        return pari.divrem(num, self.F, "b")[1] == 0

    def relation_matrix(self, m: int, d: int):
        rows = []
        for k in range(1, self.h + 1):
            Xv, Yv = self.xk(k), self.xk(m * k)
            rows.append([Xv**i for i in range(d + 1)] + [Xv**i * Yv for i in range(d + 1)])
        return pari.matrix(self.h, 2 * d + 2, [e for r in rows for e in r])
