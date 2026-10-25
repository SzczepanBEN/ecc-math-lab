"""Enumerate ALL affine smooth points of X_1(q)(F_p) and find where multipliers have a lower
minimal degree than the generic (C3) value; report the j-invariants of those exceptional points.

Prediction (CM theory, Phase 2): exceptional points are CM points whose endomorphism ring contains
an element of small norm with an eigenvalue on <P>. This is evidence about the exceptional locus
in characteristic p, not a proof of the characteristic-0 statement.

Usage: python experiments/phase3/exceptional_fp.py q p1 [p2 ...]
"""
from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1]))

from exp_b_hypotheses import monomials, nullspace_mod_p  # noqa: E402
from exp_c3_small_subgroups import conjectured_dmin  # noqa: E402
from src.cm import CLASS_NUMBER_ONE  # noqa: E402
from tateq import TateModel, pari  # noqa: E402


def tate_multiples(b0, c0, p, n):
    """x(kP), k = 1..n on y^2 + (1-c)xy - by = x^3 - bx^2 over F_p (None if kP = O)."""
    a1, a2, a3 = (1 - c0) % p, (-b0) % p, (-b0) % p

    def add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if (y1 + y2 + a1 * x2 + a3) % p == 0:
                return None
            lam = (3 * x1 * x1 + 2 * a2 * x1 - a1 * y1) * pow(2 * y1 + a1 * x1 + a3, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        nu = (y1 - lam * x1) % p
        x3 = (lam * lam + a1 * lam - a2 - x1 - x2) % p
        return (x3, (-(lam + a1) * x3 - nu - a3) % p)

    P = (0, 0)
    out, R = [], None
    for _ in range(n):
        R = add(R, P)
        out.append(R)
    return out


def j_of(b0, c0, p):
    a1, a2, a3 = (1 - c0) % p, (-b0) % p, (-b0) % p
    b2, b4, b6, b8 = (a1 * a1 + 4 * a2) % p, (a1 * a3) % p, (a3 * a3) % p, (a2 * a3 * a3) % p
    D = (-b2 * b2 * b8 - 8 * b4**3 - 27 * b6 * b6 + 9 * b2 * b4 * b6) % p
    c4 = (b2 * b2 - 24 * b4) % p
    return None if D == 0 else c4**3 * pow(D, -1, p) % p


def main():
    q = int(sys.argv[1])
    primes = [int(v) for v in sys.argv[2:]]
    T = TateModel(q)
    M = (q - 1) // 2
    d_int = M // 2
    generic = {m: conjectured_dmin(m, q, d_int)[0] for m in range(1, M + 1)}
    for p in primes:
        cm_mod_p = defaultdict(list)
        for D, j in CLASS_NUMBER_ONE.items():
            cm_mod_p[j % p].append(D)
        npts = 0
        exc = Counter()
        exc_j = defaultdict(Counter)
        for c0 in range(1, p):
            for r in pari.polrootsmod(pari.substpol(T.F, "c", c0), p):
                b0 = int(pari.lift(r))
                if b0 == 0 or j_of(b0, c0, p) is None:
                    continue
                pts = tate_multiples(b0, c0, p, q)
                assert pts[-1] is None and all(P is not None for P in pts[:-1])
                x = np.array([P[0] for P in pts[:-1]], dtype=np.int64)
                if len(set(x[:M].tolist())) < M:
                    continue
                npts += 1
                half = np.arange(1, M + 1)
                for m in range(2, M + 1):
                    d = generic[m] - 1
                    if d < 1:
                        continue
                    if nullspace_mod_p(monomials(x[half - 1], x[(m * half) % q - 1], d, 1, p), p):
                        exc[m] += 1
                        jj = j_of(b0, c0, p)
                        exc_j[m][(jj, tuple(cm_mod_p.get(jj, [])))] += 1
        print(f"q={q} p={p}: {npts} smooth affine points of X_1({q})(F_p)", flush=True)
        for m in range(2, M + 1):
            if exc[m]:
                print(f"   m={m} (generic degree {generic[m]}): {exc[m]} exceptional points; "
                      f"j (mod p) -> CM discriminants with that j: {dict(exc_j[m])}", flush=True)


if __name__ == "__main__":
    main()
