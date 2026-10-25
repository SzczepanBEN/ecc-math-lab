"""Phase 2, boundary test for conjecture C2 (small subgroups, interpolation regime).

Heuristic behind C2: a relation  Y*psi(X) + phi(X) = 0  of bidegree (d,1) has 2d+2 unknown
coefficients, and the subgroup <G'> of prime order q supplies (q-1)/2 distinct pairs
(x_k, x_{mk}) (k and q-k give the same pair). Hence:

  * if 2d+2 > (q-1)/2, EVERY m has a relation (linear algebra; a theorem);
  * if 2d+2 < (q-1)/2, a non-endomorphism m should have one only by coincidence (C2);
  * at 2d+2 = (q-1)/2, a square system: expected ~ q/p coincidences.

This script measures, for small q, the multipliers with a relation at each d (exact rank over all
(q-1)/2 distinct pairs -- an implementation independent of the batched screen in
exp_c_endomorphisms.py) and compares them with the cumulative endomorphism prediction
U_d = { alpha mod q : N(alpha) <= d } \\ {0}.

Usage:  python experiments/exp_c2_threshold.py [--out results/phase2/threshold]
"""
from __future__ import annotations

import argparse
import random
from pathlib import Path

import numpy as np

from common import save_json
from exp_b_hypotheses import monomials, nullspace_mod_p
from src.cm import CLASS_NUMBER_ONE, curve_from_j, eigenvalue_of_w, fundamental_discriminant, predicted_multipliers
from src.ec import Curve, factorize, is_probable_prime, legendre
from src.harness import provenance

MASTER_SEED = 20261027
BITS = 14
Q_RANGE = (40, 160)


def find_case(family: str, seed: int) -> dict:
    """First curve (seeded) whose group order has a prime factor q in Q_RANGE."""
    rng = random.Random(seed)
    while True:
        p = rng.randrange(1 << (BITS - 1), 1 << BITS)
        if not is_probable_prime(p):
            continue
        if family == "D=-3":
            E, D = Curve(p, 0, rng.randrange(1, p)), -3
        elif family == "D=-4":
            E, D = Curve(p, rng.randrange(1, p), 0), -4
        elif family.startswith("D="):
            D = int(family[2:])
            if CLASS_NUMBER_ONE[D] % p in (0, 1728 % p):
                continue
            c = next(c for c in range(2, p) if legendre(c, p) == -1)
            E = curve_from_j(p, CLASS_NUMBER_ONE[D], twist=rng.choice([1, c]))
        else:
            try:
                E = Curve(p, rng.randrange(1, p), rng.randrange(1, p))
            except ValueError:
                continue
            D = None
        N = E.order()
        t = p + 1 - N
        if t % p == 0:
            continue
        qs = [q for q in factorize(N) if Q_RANGE[0] <= q <= Q_RANGE[1] and N % (q * q)]
        if not qs:
            continue
        q = qs[0]
        if D is None:
            dk, _ = fundamental_discriminant(t * t - 4 * p)
            if abs(dk) <= q:        # keep the generic prediction (integers only) unambiguous
                continue
        pt_rng = random.Random(seed + 1)
        while True:
            G = E.mul(N // q, E.random_point(pt_rng))
            if G is not None:
                break
        return {"family": family, "D": D, "p": p, "a": E.a, "b": E.b, "N": N, "t": t, "q": q,
                "G": list(G), "seed": seed}


def run_case(c: dict) -> dict:
    p, q, t = c["p"], c["q"], c["t"]
    E = Curve(p, c["a"], c["b"])
    pts = E.multiples(tuple(c["G"]), q - 1)
    x = np.array([P[0] for P in pts], dtype=np.int64)
    half = np.arange(1, (q - 1) // 2 + 1, dtype=np.int64)      # one k per distinct pair
    assert len(set(x[half - 1].tolist())) == len(half)
    d_int = next(d for d in range(1, q) if 2 * d + 2 > (q - 1) // 2)
    mu = eigenvalue_of_w(p, t, c["D"] if c["D"] is not None else t * t - 4 * p, q)
    Dpred = c["D"] if c["D"] is not None else t * t - 4 * p
    rows = []
    for d in range(1, d_int + 2):
        U = set()
        for dd in range(1, d + 1):
            U |= predicted_multipliers(Dpred, dd, mu, q)
        U.discard(0)
        found = []
        for m in range(1, q):
            X, Y = x[half - 1], x[(m * half) % q - 1]
            if nullspace_mod_p(monomials(X, Y, d, 1, p), p):
                found.append(m)
        fs = set(found)
        rows.append({"d": d, "unknowns": 2 * d + 2, "constraints": (q - 1) // 2,
                     "n_found": len(fs), "n_predicted_endomorphisms": len(U),
                     "extras": sorted(fs - U), "missing": sorted(U - fs)})
    first_extra = next((r["d"] for r in rows if r["extras"]), None)
    return {"case": c, "d_interpolation": d_int, "first_d_with_extras": first_extra,
            "extras_only_at_or_above_threshold": first_extra is None or first_extra >= d_int,
            "all_m_at_threshold": next(r for r in rows if r["d"] == d_int)["n_found"] == q - 1,
            "never_missing": all(not r["missing"] for r in rows),
            "rows": rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/phase2/threshold")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    families = ["D=-3", "D=-4", "D=-7", "D=-8", "D=-11", "generic", "generic", "generic"]
    res = []
    for i, fam in enumerate(families):
        c = find_case(fam, MASTER_SEED + 100 * i)
        r = run_case(c)
        res.append(r)
        print(f"{fam:>8} p={c['p']} q={c['q']}: d_int={r['d_interpolation']} first extras at d={r['first_d_with_extras']}"
              f"  extras only >= threshold: {r['extras_only_at_or_above_threshold']}  never missing: {r['never_missing']}",
              flush=True)
    summary = {"provenance": provenance(MASTER_SEED), "q_range": Q_RANGE, "bits": BITS, "cases": res,
               "all_extras_only_at_or_above_threshold": all(r["extras_only_at_or_above_threshold"] for r in res),
               "all_never_missing": all(r["never_missing"] for r in res),
               "all_m_at_threshold": all(r["all_m_at_threshold"] for r in res)}
    save_json(summary, out / "summary.json")
    L = ["# Phase 2 threshold test (auto-generated)", "",
         "| case | p | q | d_int (2d+2 > (q-1)/2) | first d with non-endomorphism multipliers | all m at d_int | endomorphisms never missing |",
         "|---|---|---|---|---|---|---|"]
    for r in res:
        c = r["case"]
        L.append(f"| {c['family']} | {c['p']} | {c['q']} | {r['d_interpolation']} | {r['first_d_with_extras']} | "
                 f"{r['all_m_at_threshold']} | {r['never_missing']} |")
    L += ["", "Per-degree counts (found / predicted endomorphism multipliers):", ""]
    for r in res:
        c = r["case"]
        L.append(f"- {c['family']} q={c['q']}: " + ", ".join(
            f"d={row['d']}: {row['n_found']}/{row['n_predicted_endomorphisms']}" for row in r["rows"]))
    (out / "summary.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
