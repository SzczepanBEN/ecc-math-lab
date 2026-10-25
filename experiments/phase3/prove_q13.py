"""Phase 3 -- exact certificate for conjecture C3 at q = 13 (Theorem T13).

Universal object: Tate normal form E_{b,c}: y^2 + (1-c)xy - by = x^3 - bx^2, P = (0,0), over the
function field K = Q(b,c)/(F) of X_1(13), where F(b,c) is the numerator of x(6P) - x(7P).

For m in 1..6 and d >= 1 let A_{m,d} be the 6 x (2d+2) matrix with rows
    [1, X, .., X^d, Y, XY, .., X^d Y],   (X, Y) = (x(kP), x(mkP)),  k = 1..6.
A relation  Y psi(X) + phi(X) = 0 (deg <= d) holding for all 6 pairs <=> A_{m,d} has a nonzero
kernel vector. Certificate items (all exact, Q(b,c) arithmetic in PARI):

  1. F is irreducible over Q, monic in b, integral; random F_p points of F give curves on which P
     has exact order 13 (checked with the independent Phase-1 arithmetic in src/ec.py).
  2. det A_{4,2} = 0 in K                (existence of a degree-<=2 relation for m = +-4);
     rank A_{4,2} = 5 in K                (it is unique up to scaling);
     some 4x4 minor of A_{4,1} != 0 in K  (no degree-1 relation generically);
     det A_{m,2} != 0 in K for m = 2,3,5,6 (no degree-<=2 relation generically).
     (d = 3: A has 6 rows and 8 columns, so a relation exists for every m.)
  3. Integrality: x(kP) in Z[b,c][1/(c (b-c) (b-c^2-c))] and F divides (the primitive integral
     numerator of) det A_{4,2} in Z[b,c]; the denominators meet F only at c = 0 (forcing b = 0,
     the singular cubic). Hence the existence statement holds for every elliptic curve with a
     point of order 13 over every field.
  4. Exceptional loci: zeros on X_1(13) of det A_{m,2} (m = 2,3,5,6) and the common zeros of all
     4x4 minors of A_{4,1}, classified as cusps (Delta = 0) or elliptic curves with their j-invariants
     (tested against Hilbert class polynomials).

Usage:  python experiments/phase3/prove_q13.py [--out results/phase3/q13]
"""
from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from exceptional import Delta, c4, cm_disc, exceptional_points, y  # noqa: E402
from tate13 import F, X, pari, reduce_mod_F, relation_matrix  # noqa: E402

from src.ec import Curve  # noqa: E402

b, c = pari("b"), pari("c")


def zero_in_K(expr) -> bool:
    return reduce_mod_F(pari.numerator(expr)) == 0


def item1(rng) -> dict:
    fa = pari.factor(F)
    hits = 0
    for p in (10007, 10009, 10037, 10039, 100003):
        for _ in range(30):
            c0 = rng.randrange(1, p)
            for r in pari.polrootsmod(pari.substpol(F, "c", c0), p):
                b0 = int(pari.lift(r))
                if b0 == 0:
                    continue
                a1, a2, a3 = (1 - c0) % p, (-b0) % p, (-b0) % p
                b2, b4, b6 = (a1 * a1 + 4 * a2) % p, (a1 * a3) % p, (a3 * a3) % p
                c4_, c6_ = (b2 * b2 - 24 * b4) % p, (-b2**3 + 36 * b2 * b4 - 216 * b6) % p
                try:
                    E = Curve(p, -27 * c4_ % p, -54 * c6_ % p)
                except ValueError:
                    continue
                Pp = (3 * b2 % p, 108 * a3 % p)          # image of (0,0) under x' = 36x + 3b2, y' = 108(2y + a1x + a3)
                assert E.is_on_curve(Pp)
                assert E.mul(13, Pp) is None and all(E.mul(k, Pp) is not None for k in range(1, 13))
                hits += 1
    return {"F": str(F), "bidegree_b_c": [int(pari.poldegree(F, "b")), int(pari.poldegree(F, "c"))],
            "irreducible_over_Q": len(fa[0]) == 1 and int(fa[1][0]) == 1,
            "monic_in_b": pari.pollead(F, "b") == 1, "integral": pari.denominator(pari.content(F)) == 1,
            "Fp_points_checked_exact_order_13": hits}


def item2() -> dict:
    out = {}
    A42 = relation_matrix(4, 2)
    out["det_A_4_2_is_zero_in_K"] = zero_in_K(pari.matdet(A42))
    rank5 = any(not zero_in_K(pari.matdet(pari.matrix(5, 5, [A42[i, j] for i in r for j in cc])))
                for r in itertools.combinations(range(6), 5) for cc in itertools.combinations(range(6), 5))
    out["rank_A_4_2_is_5"] = rank5
    A41 = relation_matrix(4, 1)
    out["A_4_1_has_nonzero_4x4_minor"] = any(
        not zero_in_K(pari.matdet(pari.matrix(4, 4, [A41[i, j] for i in r for j in range(4)])))
        for r in itertools.combinations(range(6), 4))
    for m in (2, 3, 5, 6):
        out[f"det_A_{m}_2_is_nonzero_in_K"] = not zero_in_K(pari.matdet(relation_matrix(m, 2)))
    out["m1_degree1"] = True  # x(P) = x(P)
    return out


def integral_poly(expr) -> bool:
    """True iff expr is a polynomial in b, c with integer coefficients."""
    if expr == 0 or str(pari.type(expr)) == "t_INT":
        return True
    if str(pari.type(expr)) != "t_POL":
        return False
    for i in range(int(pari.poldegree(expr, "b")) + 1):
        cb = pari.polcoef(expr, i, "b")
        if cb == 0 or str(pari.type(cb)) == "t_INT":
            continue
        if str(pari.type(cb)) != "t_POL":
            return False
        for j in range(int(pari.poldegree(cb, "c")) + 1):
            if str(pari.type(pari.polcoef(cb, j, "c"))) != "t_INT":
                return False
    return True


def item3() -> dict:
    allowed = [c, b - c, b - c**2 - c]
    D = (c * (b - c) * (b - c**2 - c)) ** 8
    den_ok = all(integral_poly(X[k] * D) for k in range(1, 7))
    meets = {str(g): str(pari.factor(pari.polresultant(F, g, "b"))) if pari.poldegree(g, "b") > 0 else "c"
             for g in allowed}
    num4 = pari.numerator(pari.matdet(relation_matrix(4, 2)))
    num4p = num4 / pari.content(num4)                      # primitive integral polynomial
    q, r = pari.divrem(num4p, F, "b")
    return {"x_kP_denominators_only_c_bminusc_bminusc2minusc_and_integral_numerators": den_ok,
            "resultant_F_with_denominator_factors": meets,
            "F_divides_primitive_numerator_of_det_A42_over_Z": r == 0 and pari.denominator(pari.content(q)) == 1,
            "Delta_factorization": str(pari.factor(Delta))}


def item4() -> dict:
    out = {}
    for m in (2, 3, 5, 6):
        out[f"m={m}, d=2"] = exceptional_points(pari.numerator(pari.matdet(relation_matrix(m, 2))))
    A = relation_matrix(4, 1)
    minors = [pari.numerator(pari.matdet(pari.matrix(4, 4, [A[i, j] for i in r for j in range(4)])))
              for r in itertools.combinations(range(6), 4)]
    comps = []
    fa = pari.factor(pari.polresultant(F, minors[0], "b"))
    for i in range(len(fa[0])):
        hc = fa[0][i]
        if hc == c:
            comps.append({"component": "c = 0 (singular cubic)"})
            continue
        hy = pari.substpol(hc, "c", y)
        nf = pari.nfinit(hy)
        for g in pari.nffactor(nf, pari.substpol(F, "c", y))[0]:
            if pari.poldegree(g, "b") != 1:
                continue
            b0 = pari.lift(-pari.polcoef(g, 0, "b") / pari.polcoef(g, 1, "b"))

            def ev(f):
                return pari.Mod(pari.substpol(pari.substpol(f, "b", b0), "c", y), hy)
            if ev(minors[0]) != 0:
                continue
            D0 = ev(Delta)
            entry = {"component": f"c-field degree {int(pari.poldegree(hc, 'c'))}",
                     "all_15_minors_vanish": all(ev(mn) == 0 for mn in minors), "cusp": D0 == 0}
            if D0 != 0:
                mp = pari.minpoly(ev(c4)**3 / D0, "x")
                entry.update({"j_minpoly": str(mp), "CM_discriminant": cm_disc(mp)})
            comps.append(entry)
    out["m=4, d=1 (all 4x4 minors)"] = comps
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/phase3/q13")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    cert = {"theorem": "T13 (C3 at q = 13)", "item1_model": item1(random.Random(13)),
            "item2_ranks": item2(), "item3_integrality": item3(), "item4_exceptional": item4()}
    i2 = cert["item2_ranks"]
    i3 = cert["item3_integrality"]
    cert["all_checks_pass"] = bool(
        cert["item1_model"]["irreducible_over_Q"] and cert["item1_model"]["monic_in_b"]
        and cert["item1_model"]["Fp_points_checked_exact_order_13"] > 0 and all(i2.values())
        and i3["x_kP_denominators_only_c_bminusc_bminusc2minusc_and_integral_numerators"]
        and i3["F_divides_primitive_numerator_of_det_A42_over_Z"])
    cert["runtime_s"] = round(time.time() - t0, 1)
    (out / "certificate.json").write_text(json.dumps(cert, indent=2, default=str))
    print(json.dumps(cert, indent=2, default=str))


if __name__ == "__main__":
    main()
