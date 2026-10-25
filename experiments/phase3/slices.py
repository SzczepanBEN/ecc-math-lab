"""Slice certificates for statements "G = 0 in K" / "G != 0 in K", K = Q(X_1(q)) = Q(b,c)/(F).

For a minor of the relation matrix, clear row denominators to get a polynomial matrix P (row
denominators are products of c, b - c, b - c^2 - c and similar factors that are nonzero on the
smooth part of F = 0, checked by the caller), so the minor vanishes in K iff F | G := det P[R, C].

* Nonvanishing: one rational c0 with G(b, c0) != 0 mod F(b, c0) proves F does not divide G (F is
  monic in b, so division commutes with specialising c).
* Vanishing (Bezout): if F does not divide G, the affine curves F = 0 and G = 0 (F irreducible) meet
  in at most deg F * deg G points. Each slice c = c0 with F(b, c0) squarefree on which
  G(b, c0) = 0 mod F(b, c0) gives deg_b F distinct common points. Enough slices prove F | G.
  deg G is bounded by sum_j max_i totdeg(P_ij).
Everything is exact over Q (univariate in b after substituting c = c0).
"""
from __future__ import annotations

import itertools

from tateq import TateModel, pari

t = pari("varlower(\"t\")")


def totdeg(f) -> int:
    if f == 0:
        return -1
    g = pari.substpol(pari.substpol(f, "b", pari("b") * t), "c", pari("c") * t)
    return int(pari.poldegree(g, t))


class SliceProver:
    def __init__(self, T: TateModel):
        self.T = T
        self.F = T.F
        self.nb = int(pari.poldegree(self.F, "b"))
        self.degF = totdeg(self.F)
        self.cache = {}

    def poly_matrix(self, m: int, d: int):
        A = self.T.relation_matrix(m, d)
        rows, cols = int(pari.matsize(A)[0]), int(pari.matsize(A)[1])
        P, dens = [], []
        for i in range(rows):
            den = pari(1)
            for j in range(cols):
                den = pari.lcm(den, pari.denominator(A[i, j]))
            dens.append(den)
            P.append([pari.simplify(A[i, j] * den) for j in range(cols)])
        return P, dens

    def slice_ok(self, c0) -> bool:
        Fc = pari.substpol(self.F, "c", c0)
        return pari.poldisc(Fc, "b") != 0

    def minor_on_slice(self, P, R, C, c0):
        M = pari.matrix(len(R), len(C), [pari.substpol(P[i][j], "c", c0) for i in R for j in C])
        G0 = pari.matdet(M)
        return pari.divrem(G0, pari.substpol(self.F, "c", c0), "b")[1]

    def degree_bound(self, P, R, C) -> int:
        return sum(max(totdeg(P[i][j]) for i in R) for j in C)

    def certify_vanishing(self, P, R, C, start=2):
        """Return (True, slices_used, bound) if F | det P[R,C] is certified, or (False, c0, ...) at a
        slice where the remainder is nonzero (which proves F does not divide it)."""
        bound = self.degF * self.degree_bound(P, R, C)
        need = bound // self.nb + 1
        used, c0 = 0, start
        while used < need:
            if self.slice_ok(c0):
                if self.minor_on_slice(P, R, C, c0) != 0:
                    return False, c0, bound
                used += 1
            c0 += 1
        return True, used, bound

    def certify_nonvanishing(self, P, R, C, tries=10, start=2):
        c0, done = start, 0
        while done < tries:
            if self.slice_ok(c0):
                if self.minor_on_slice(P, R, C, c0) != 0:
                    return True, c0
                done += 1
            c0 += 1
        return False, None


def all_minors(rows: int, cols: int, k: int):
    return [(R, C) for R in itertools.combinations(range(rows), k) for C in itertools.combinations(range(cols), k)]
