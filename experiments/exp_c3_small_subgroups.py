"""Phase 2 follow-up: non-endomorphism multipliers on small subgroups (unexpected result of
exp_c2_threshold.py).

For generic curves with a subgroup of small prime order q we compute, for EVERY multiplier m, the
minimal degree d_min(m) of a rational map R with x(mP) = R(x(P)) on <G'> (exact rank over all
M = (q-1)/2 distinct pairs; monotone in d, so a binary search over d in [1, d_int] is exact).
Two different curves are used for each q, to see whether d_min depends on the curve.

Usage:  python experiments/exp_c3_small_subgroups.py [--out results/phase2/small_subgroups]
"""
from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

import numpy as np

from common import save_json
from exp_b_hypotheses import monomials, nullspace_mod_p
from src.cm import fundamental_discriminant
from src.ec import Curve, factorize, is_probable_prime
from src.harness import provenance

MASTER_SEED = 20261028
BITS = 14
Q_LIST = [61, 73, 97, 109, 157, 59, 71, 83, 101, 107]      # five q = 1 mod 6, five q = 5 mod 6
CURVES_PER_Q = 2


def find_curve_with_q(q: int, seed: int) -> dict:
    rng = random.Random(seed)
    while True:
        p = rng.randrange(1 << (BITS - 1), 1 << BITS)
        if not is_probable_prime(p):
            continue
        try:
            E = Curve(p, rng.randrange(1, p), rng.randrange(1, p))
        except ValueError:
            continue
        if E.j_invariant in (0, 1728 % p):
            continue
        N = E.order()
        t = p + 1 - N
        if N % q or N % (q * q) == 0 or t % p == 0:    # q || N exactly
            continue
        dk, f = fundamental_discriminant(t * t - 4 * p)
        if abs(dk) <= 4 * q:          # endomorphisms of norm < q must be integers
            continue
        pr = random.Random(seed + 1)
        while True:
            G = E.mul(N // q, E.random_point(pr))
            if G is not None:
                break
        return {"p": p, "a": E.a, "b": E.b, "N": N, "t": t, "q": q, "G": list(G), "seed": seed,
                "D_K": dk, "conductor_Z_pi": f}


def min_degrees(c: dict) -> tuple[dict[int, int], int]:
    p, q = c["p"], c["q"]
    E = Curve(p, c["a"], c["b"])
    x = np.array([P[0] for P in E.multiples(tuple(c["G"]), q - 1)], dtype=np.int64)
    M = (q - 1) // 2
    half = np.arange(1, M + 1, dtype=np.int64)
    d_int = next(d for d in range(1, q) if 2 * d + 2 > M)

    def has(m, d):
        return bool(nullspace_mod_p(monomials(x[half - 1], x[(m * half) % q - 1], d, 1, p), p))

    out = {}
    for m in range(1, M + 1):                 # m and q - m behave identically (x is even)
        lo, hi = 1, d_int                     # has(m, d_int) is always True
        assert has(m, d_int)
        while lo < hi:
            mid = (lo + hi) // 2
            if has(m, mid):
                hi = mid
            else:
                lo = mid + 1
        out[m] = lo
    return out, d_int


def conjectured_dmin(m: int, q: int, d_int: int) -> tuple[int, tuple[int, int] | None]:
    """Conjecture C3 (fitted on Q_LIST, frozen before the held-out run):

        d_min(m) = min(d_int, min over (a, b): b odd, gcd(a, b) = 1, b m = +-a (mod q),
                                q = 4a^2 (mod b)  of  ((b-1) q + 4a^2 - eps_b) / (4b)),
        eps_b = 0 if b = 1 (mod 4), 2b if b = 3 (mod 4).
    """
    best, arg = d_int, None
    for b in range(1, q, 2):
        r = (b * m) % q
        for a in (r, r - q):
            if not a or math.gcd(abs(a), b) != 1:
                continue
            num = (b - 1) * q + 4 * a * a - (2 * b if b % 4 == 3 else 0)
            if num % (4 * b) == 0 and num // (4 * b) < best:
                best, arg = num // (4 * b), (a, b)
    return best, arg


def small_fractions(m: int, q: int, H: int = 12) -> list[tuple[int, int]]:
    """All (a, b) with |a| <= H, 1 <= b <= H, gcd(a, b) = 1, a = b m (mod q) up to sign of m."""
    res = []
    for b in range(1, H + 1):
        for a in range(-H, H + 1):
            if a and math.gcd(a, b) == 1 and ((a - b * m) % q == 0 or (a + b * m) % q == 0):
                res.append((a, b))
    return sorted(res, key=lambda ab: (ab[0] ** 2 + ab[1] ** 2, ab[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/phase2/small_subgroups")
    ap.add_argument("--qs", default=",".join(map(str, Q_LIST)))
    ap.add_argument("--curves-per-q", type=int, default=CURVES_PER_Q)
    ap.add_argument("--seed", type=int, default=MASTER_SEED)
    ap.add_argument("--label", default="fit")
    ap.add_argument("--cm", default="", help="comma list of CM discriminants for a side check (not part of C3)")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    q_list = [int(v) for v in args.qs.split(",")]
    cases = []
    for i, q in enumerate(q_list):
        for j in range(args.curves_per_q):
            c = find_curve_with_q(q, args.seed + 1000 * i + 100 * j)
            dmin, d_int = min_degrees(c)
            M = (q - 1) // 2
            extras = [{"m": m, "d": d, "fractions": small_fractions(m, q)[:3]}
                      for m, d in sorted(dmin.items(), key=lambda kv: kv[1])
                      if d < d_int and not (d == m * m)]
            pred = {m: conjectured_dmin(m, q, d_int) for m in dmin}
            n_match = sum(1 for m in dmin if pred[m][0] == dmin[m])
            n_nontriv = sum(1 for m in dmin if pred[m][0] < d_int)
            cases.append({"curve": c, "M": M, "d_int": d_int, "dmin": dmin, "extras": extras,
                          "conjecture_C3": {"n_multipliers": len(dmin), "n_match": n_match,
                                            "n_predicted_below_threshold": n_nontriv,
                                            "n_observed_below_threshold": sum(1 for d in dmin.values() if d < d_int),
                                            "mismatches": [{"m": m, "observed": dmin[m], "predicted": pred[m][0],
                                                            "rep": pred[m][1]} for m in dmin if pred[m][0] != dmin[m]]}})
            print(f"q={q} (q mod 6 = {q % 6}) p={c['p']}: d_int={d_int}, below-threshold extras: "
                  + ", ".join(f"m={e['m']}:{e['d']}{tuple(e['fractions'][0]) if e['fractions'] else ''}"
                              for e in extras[:8]), flush=True)
    # curve-(in)dependence: same q, several curves
    same = {}
    k = args.curves_per_q
    for i, q in enumerate(q_list):
        profs = [cases[k * i + j]["dmin"] for j in range(k)]
        same[q] = {"identical_dmin_profile": all(pp == profs[0] for pp in profs),
                   "n_m_differing": sum(1 for m in profs[0] if any(pp[m] != profs[0][m] for pp in profs))}
    # the b = 3 formula: d(+-a/3) = (M + a^2 - 1)/3 for q = 1 mod 6
    checks = []
    for cs in cases:
        q, M, d_int = cs["curve"]["q"], cs["M"], cs["d_int"]
        for a in range(1, 9):
            if a % 3 == 0:
                continue
            m = a * pow(3, -1, q) % q
            m = min(m, q - m)
            pred = (M + a * a - 1) / 3
            checks.append({"q": q, "p": cs["curve"]["p"], "a": a, "m": m, "observed": cs["dmin"][m],
                           "formula": pred, "formula_below_threshold": pred < d_int,
                           "match": cs["dmin"][m] == pred})
    cm_side = []
    if args.cm:
        from exp_c2_threshold import find_case
        for i, D in enumerate(int(v) for v in args.cm.split(",")):
            c = find_case(f"D={D}", args.seed + 50000 + 100 * i)
            dmin, d_int = min_degrees(c)
            pred = {m: conjectured_dmin(m, c["q"], d_int)[0] for m in dmin}
            cm_side.append({"D": D, "p": c["p"], "q": c["q"], "d_int": d_int,
                            "n_multipliers": len(dmin), "n_match_C3": sum(pred[m] == dmin[m] for m in dmin),
                            "n_observed_lower_than_C3": sum(dmin[m] < pred[m] for m in dmin),
                            "n_observed_higher_than_C3": sum(dmin[m] > pred[m] for m in dmin)})
            print(f"CM side check D={D} q={c['q']}: C3 matches {cm_side[-1]['n_match_C3']}/{len(dmin)}, "
                  f"observed lower {cm_side[-1]['n_observed_lower_than_C3']}, higher {cm_side[-1]['n_observed_higher_than_C3']}", flush=True)
    tot = sum(cs["conjecture_C3"]["n_multipliers"] for cs in cases)
    ok = sum(cs["conjecture_C3"]["n_match"] for cs in cases)
    nontriv = sum(cs["conjecture_C3"]["n_observed_below_threshold"] for cs in cases)
    summary = {"provenance": provenance(args.seed), "label": args.label, "q_list": q_list,
               "curves_per_q": k, "cases": cases, "curve_dependence": same, "b3_formula_checks": checks,
               "cm_side_check": cm_side,
               "conjecture_C3_total": {"multipliers": tot, "exact_matches": ok,
                                       "observed_below_threshold": nontriv}}
    print(f"C3: {ok}/{tot} exact matches; {nontriv} multipliers observed below threshold", flush=True)
    save_json(summary, out / "summary.json")
    L = ["# Small-subgroup multipliers (auto-generated)", "",
         f"Run label: {args.label}. Conjecture C3: {ok}/{tot} multipliers match exactly "
         f"({nontriv} of them have a relation below the interpolation threshold).", "",
         "| q | q mod 6 | p (curve 1 / curve 2) | d_int | identical d_min profile for both curves | below-threshold extras (curve 1): m:d (a/b) |",
         "|---|---|---|---|---|---|"]
    for i, q in enumerate(q_list):
        c1, c2 = cases[k * i], cases[k * i + k - 1]
        ex = ", ".join(f"{e['m']}:{e['d']} ({e['fractions'][0][0]}/{e['fractions'][0][1]})" if e["fractions"]
                       else f"{e['m']}:{e['d']}" for e in c1["extras"])
        L.append(f"| {q} | {q % 6} | {c1['curve']['p']} / {c2['curve']['p']} | {c1['d_int']} | "
                 f"{same[q]['identical_dmin_profile'] if k > 1 else 'n/a (1 curve)'} ({same[q]['n_m_differing']} m differ) | {ex} |")
    L += ["", "## b = 3 formula d(+-a/3) = (M + a^2 - 1)/3", "",
          "| q | p | a | m | observed d_min | formula | formula < d_int | match |", "|---|---|---|---|---|---|---|---|"]
    for ck in checks:
        if ck["formula_below_threshold"] or ck["observed"] < (ck["q"] - 1) // 4:
            L.append(f"| {ck['q']} | {ck['p']} | {ck['a']} | {ck['m']} | {ck['observed']} | {ck['formula']:.2f} | "
                     f"{ck['formula_below_threshold']} | {ck['match']} |")
    (out / "summary.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
