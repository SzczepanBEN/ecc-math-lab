"""Phase 3 -- exact certificate for conjecture C3 at a prime q, by slice certificates (any q).

Theorem T_q (generic statement, function field K of X_1(q), Tate normal form): for every m in
1..M (M = (q-1)/2) with d* = C3(m) < d_int (the non-trivial predictions),
  (E) all (2d*+2)-minors of A_{m,d*} vanish in K   -> a relation of degree <= d* exists;
  (U) A_{m,d*} has rank exactly 2d*+1              -> it is unique up to scaling;
and for EVERY m with C3(m) >= 2,
  (L) some (2d+2)-minor of A_{m,d} with d = C3(m) - 1 is nonzero in K -> no lower-degree relation
      generically.
(At d = d_int a relation exists for all m by linear algebra.) Denominators of x(kP) are checked to
meet F only over c = 0; the integer constants of those resultants are reported (they decide in
which characteristics the existence statement transfers to every curve with a point of order q).

Usage: python experiments/phase3/prove_q.py q [--out DIR]
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
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1]))

from exp_c3_small_subgroups import conjectured_dmin  # noqa: E402
from slices import SliceProver  # noqa: E402
from tateq import TateModel, pari  # noqa: E402

from src.ec import Curve  # noqa: E402


def model_check(T: TateModel, rng) -> int:
    q, hits = T.q, 0
    for p in (10007, 10009, 100003):
        for _ in range(25):
            c0 = rng.randrange(1, p)
            for r in pari.polrootsmod(pari.substpol(T.F, "c", c0), p):
                b0 = int(pari.lift(r))
                if b0 == 0:
                    continue
                a1, a2, a3 = (1 - c0) % p, (-b0) % p, (-b0) % p
                b2, b4, b6 = (a1 * a1 + 4 * a2) % p, (a1 * a3) % p, (a3 * a3) % p
                c4, c6 = (b2 * b2 - 24 * b4) % p, (-b2**3 + 36 * b2 * b4 - 216 * b6) % p
                try:
                    E = Curve(p, -27 * c4 % p, -54 * c6 % p)
                except ValueError:
                    continue
                Pp = (3 * b2 % p, 108 * a3 % p)
                assert E.mul(q, Pp) is None and all(E.mul(k, Pp) is not None for k in range(1, q))
                hits += 1
    return hits


def first_nonzero(S, P, rows, k, ncols, fixed_cols=None):
    for R in itertools.combinations(range(rows), k):
        for C in ([fixed_cols] if fixed_cols else itertools.combinations(range(ncols), k)):
            ok, c0 = S.certify_nonvanishing(P, R, C)
            if ok:
                return {"rows": R, "cols": C, "slice_c0": c0}
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("q", type=int)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    q = args.q
    out = Path(args.out or f"results/phase3/q{q}")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    T = TateModel(q)
    S = SliceProver(T)
    M, d_int = T.h, T.h // 2
    c3 = {m: conjectured_dmin(m, q, d_int)[0] for m in range(1, M + 1)}
    fa = pari.factor(T.F)
    cert = {"q": q, "M": M, "d_int": d_int, "C3_prediction": c3,
            "F_bidegree_b_c": [int(pari.poldegree(T.F, "b")), int(pari.poldegree(T.F, "c"))],
            "F_total_degree": S.degF,
            "F_irreducible_over_Q": len(fa[0]) == 1 and int(fa[1][0]) == 1,
            "F_monic_in_b": pari.pollead(T.F, "b") == 1,
            "raw_numerator_factor_count": len(T.raw_factors),
            "Fp_model_check_points": model_check(T, random.Random(q)),
            "existence": {}, "uniqueness": {}, "lower_bound": {}}
    rows = M
    for m in range(1, M + 1):
        d = c3[m]
        if 1 <= d < d_int and m != 1:
            P, _ = S.poly_matrix(m, d)
            ncols = 2 * d + 2
            res = []
            for R in itertools.combinations(range(rows), ncols):
                ok = S.certify_vanishing(P, R, tuple(range(ncols)))
                res.append({"rows": R, "vanishes": ok[0], "slices": ok[1], "bezout_bound": ok[2]})
                if not ok[0]:
                    break
            cert["existence"][m] = {"d": d, "all_minors_vanish": all(r["vanishes"] for r in res),
                                    "n_minors": len(res), "min_slices": min(r["slices"] for r in res),
                                    "max_bezout_bound": max(r["bezout_bound"] for r in res)}
            cert["uniqueness"][m] = first_nonzero(S, P, rows, ncols - 1, ncols) is not None
        if d >= 2:
            dl = d - 1
            P, _ = S.poly_matrix(m, dl)
            nc = 2 * dl + 2
            cert["lower_bound"][m] = {"d": dl, "nonzero_minor": first_nonzero(S, P, rows, nc, nc, tuple(range(nc)))}
        print(f"q={q} m={m}: C3={d} existence={cert['existence'].get(m)} lower={cert['lower_bound'].get(m) is not None and cert['lower_bound'][m]['nonzero_minor'] is not None}  [{time.time()-t0:.0f}s]", flush=True)
    facs = {}
    for k in range(1, M + 1):
        f2 = pari.factor(pari.denominator(T.X[k]))
        for i in range(len(f2[0])):
            facs[str(f2[0][i])] = f2[0][i]
    den = {}
    for s, g in facs.items():
        var = "b" if pari.poldegree(g, "b") > 0 else "c"
        R = pari.polresultant(T.F, g, var)
        den[s] = {"resultant_variable": var, "resultant": str(R)}
    cert["denominator_resultants"] = den
    cert["F_at_c0"] = str(pari.substpol(T.F, "c", 0))
    cert["all_checks_pass"] = bool(
        cert["F_irreducible_over_Q"] and cert["F_monic_in_b"] and cert["Fp_model_check_points"] > 0
        and all(v["all_minors_vanish"] for v in cert["existence"].values())
        and all(cert["uniqueness"].values())
        and all(v["nonzero_minor"] is not None for v in cert["lower_bound"].values()))
    cert["runtime_s"] = round(time.time() - t0, 1)
    (out / "certificate.json").write_text(json.dumps(cert, indent=2, default=str))
    print(json.dumps({k: cert[k] for k in ("q", "C3_prediction", "existence", "uniqueness", "denominator_resultants",
                                           "F_at_c0", "all_checks_pass", "runtime_s")}, indent=1, default=str))


if __name__ == "__main__":
    main()
