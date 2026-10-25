"""Exceptional points on X_1(13): where the generic statements of Theorem T13 fail."""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from tate13 import pari, F, relation_matrix

b, c = pari('b'), pari('c')
y = pari('varlower("w")')      # lowest-priority variable for number fields
a1, a2, a3 = 1 - c, -b, -b
b2 = a1**2 + 4*a2; b4 = a1*a3; b6 = a3**2; b8 = a2*a3**2
Delta = -b2**2*b8 - 8*b4**3 - 27*b6**2 + 9*b2*b4*b6
c4 = b2**2 - 24*b4
jfun = c4**3 / Delta
BYH = {}
for D in range(-3, -5001, -1):
    if D % 4 in (0, 1):
        BYH.setdefault(int(pari.qfbclassno(D)), []).append(D)


def cm_disc(minpoly):
    for D in BYH.get(int(pari.poldegree(minpoly)), []):
        if pari.polclass(D) == minpoly:
            return D
    return None


def exceptional_points(num):
    """Irreducible components of {F = 0, num = 0} away from c = 0, with j and CM test."""
    out = []
    fa = pari.factor(pari.polresultant(F, num, 'b'))
    for i in range(len(fa[0])):
        hc = fa[0][i]
        if pari.poldegree(hc, 'c') == 1 and hc == c:
            out.append({"component": "c = 0 (singular cubic, b = 0)", "multiplicity": int(fa[1][i])})
            continue
        hy = pari.substpol(hc, 'c', y)
        nf = pari.nfinit(hy)
        Fy = pari.substpol(F, 'c', y)
        Ny = pari.substpol(num, 'c', y)
        for g in pari.nffactor(nf, Fy)[0]:
            # does g divide num over the number field?  reduce num mod g with coefficients mod hy
            rem = pari.lift(pari.lift(pari.Mod(pari.Mod(Ny, hy), g)))  # remainder in b with coeffs mod hy
            if rem != 0 or pari.poldegree(g, 'b') != 1:
                if rem == 0:
                    out.append({"component": f"deg_c {pari.poldegree(hc,'c')}", "note": f"common factor of b-degree {pari.poldegree(g,'b')}"})
                continue
            b0 = pari.lift(-pari.polcoef(g, 0, 'b') / pari.polcoef(g, 1, 'b'))
            ev = lambda f: pari.Mod(pari.substpol(pari.substpol(f, 'b', b0), 'c', y), hy)
            D0, c40 = ev(Delta), ev(c4)
            if D0 == 0:
                out.append({"component": f"c-field degree {int(pari.poldegree(hc,'c'))}", "multiplicity": int(fa[1][i]),
                            "type": "CUSP (Delta = 0: degenerate curve)", "c4_zero": bool(c40 == 0)})
                continue
            j0 = c40**3 / D0
            mp = pari.minpoly(j0, 'x')
            out.append({"component": f"c-field degree {int(pari.poldegree(hc,'c'))}", "multiplicity": int(fa[1][i]),
                        "type": "elliptic curve", "j_minpoly_degree": int(pari.poldegree(mp)), "j_minpoly": str(mp),
                        "CM_discriminant": cm_disc(mp)})
    return out


if __name__ == "__main__":
    for m in (2, 3, 5, 6):
        num = pari.numerator(pari.matdet(relation_matrix(m, 2)))
        for e in exceptional_points(num):
            print(f"m={m}:", e, flush=True)
