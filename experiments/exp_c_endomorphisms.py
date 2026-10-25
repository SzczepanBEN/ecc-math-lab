"""Phase 2, direction 1 -- blind discovery of small-degree endomorphisms.

For each curve E/F_p and a point G' of prime order q we take x_k = x(kG'), k = 1..q-1, and for
every multiplier m in [1, q-1] and every degree d <= DMAX we decide exactly (over F_p) whether

        x_{mk} = R(x_k)   for all k,  for one rational function R of degree <= d,

i.e. whether a relation  Y*psi(X) + phi(X) = 0  with deg phi, deg psi <= d holds on every pair
(x_k, x_{mk}). We record the minimal such d for every m that has one ("found multipliers").

Prediction (see research_log.md, Phase 2): if End(E) = O (an imaginary quadratic order of
discriminant D) then the multipliers of minimal degree exactly d are
        S_d = { alpha mod q : alpha in O, N(alpha) = d },   |S_d| = r_O(d),
where the inclusion "S_d has degree <= d" is a theorem and "nothing else" is conjecture C2.
The prediction is computed from (p, t, q, D) alone, independently of the search.

Curves: one curve for each of the 13 class-number-one orders (j = j(O) mod p) plus generic
controls. Output: results/phase2/.

Usage:  python experiments/exp_c_endomorphisms.py [--dmax 20] [--dmax163 42]
"""
from __future__ import annotations

import argparse
import math
import random
import time

import numpy as np

from common import RESULTS, save_json
from exp_b_hypotheses import evaluate, monomials, normalize, nullspace_mod_p
from src.cm import (CLASS_NUMBER_ONE, curve_from_j, eigenvalue_of_w, fundamental_discriminant,
                    kernel_structure, pdeg, pgcd, pmul, predicted_multipliers, r_O,
                    trace_for_multiplier)
from src.ec import Curve, factorize, is_probable_prime, legendre
from src.harness import provenance

OUT = RESULTS / "phase2"
MASTER_SEED = 20261026
BITS = 14          # p in [2^13, 2^14)
MAX_COFACTOR = 8


# ------------------------------------------------------------------ curve construction


def _accept(E: Curve):
    N = E.order()
    t = E.p + 1 - N
    if t % E.p == 0:
        return None  # supersingular: p does not split in K
    q = max(factorize(N))
    h = N // q
    if h > MAX_COFACTOR or not is_probable_prime(q):
        return None
    return {"p": E.p, "a": E.a, "b": E.b, "N": N, "t": t, "q": q, "h": h}


def build_cm_curve(D: int, seed: int) -> dict:
    j = CLASS_NUMBER_ONE[D]
    rng = random.Random(seed)
    lo, hi = 1 << (BITS - 1), 1 << BITS
    while True:
        p = rng.randrange(lo, hi)
        if not is_probable_prime(p) or p < 5:
            continue
        cands = []
        if D == -3:
            cands = [Curve(p, 0, rng.randrange(1, p)) for _ in range(12)]
        elif D == -4:
            cands = [Curve(p, rng.randrange(1, p), 0) for _ in range(12)]
        else:
            if j % p in (0, 1728 % p):
                continue
            c = next(c for c in range(2, p) if legendre(c, p) == -1)
            cands = [curve_from_j(p, j), curve_from_j(p, j, twist=c)]
        for E in cands:
            r = _accept(E)
            if r:
                r.update({"family": f"CM D={D}", "D_table": D, "j_table": j, "seed": seed,
                          "j_mod_p": E.j_invariant})
                return r


def build_generic_curve(seed: int, dmax: int) -> dict:
    rng = random.Random(seed)
    lo, hi = 1 << (BITS - 1), 1 << BITS
    while True:
        p = rng.randrange(lo, hi)
        if not is_probable_prime(p):
            continue
        try:
            E = Curve(p, rng.randrange(1, p), rng.randrange(1, p))
        except ValueError:
            continue
        if E.j_invariant in (0, 1728 % p):
            continue
        r = _accept(E)
        if not r:
            continue
        dk, _ = fundamental_discriminant(r["t"] ** 2 - 4 * p)
        if abs(dk) <= 4 * dmax:     # keep only controls whose prediction is unambiguous
            continue
        r.update({"family": "generic", "D_table": None, "j_table": None, "seed": seed,
                  "j_mod_p": E.j_invariant, "D_K_of_frobenius": dk})
        return r


def subgroup_sequence(c: dict, seed: int) -> np.ndarray:
    E = Curve(c["p"], c["a"], c["b"])
    rng = random.Random(seed)
    while True:
        G = E.mul(c["h"], E.random_point(rng))
        if G is not None:
            break
    assert E.mul(c["q"], G) is None
    c["G"] = list(G)
    pts = E.multiples(G, c["q"] - 1)
    return np.array([P[0] for P in pts], dtype=np.int64)


# ------------------------------------------------------------------ batched exact search


def modinv_vec(a: np.ndarray, p: int) -> np.ndarray:
    res = np.ones_like(a)
    base = a % p
    e = p - 2
    while e:
        if e & 1:
            res = res * base % p
        base = base * base % p
        e >>= 1
    return res


def batch_rank_deficient(B: np.ndarray, p: int) -> np.ndarray:
    """For a stack B (M, r, c) with r >= c over F_p: True where rank(B[m]) < c.

    Column-by-column elimination: if column j has no pivot among rows >= j after eliminating
    columns < j, column j lies in the span of the previous columns, so rank < c.
    """
    B = B.copy()
    M, r, c = B.shape
    ar = np.arange(M)
    deficient = np.zeros(M, dtype=bool)
    for j in range(c):
        nz = B[:, j:, j] != 0
        has = nz.any(axis=1)
        deficient |= ~has
        piv = j + nz.argmax(axis=1)
        rj = B[ar, j, :].copy()
        B[ar, j, :] = B[ar, piv, :]
        B[ar, piv, :] = rj
        pv = np.where(has, B[:, j, j], 1)
        B[:, j, j:] = B[:, j, j:] * modinv_vec(pv, p)[:, None] % p
        if j + 1 < r:
            f = B[:, j + 1:, j]
            B[:, j + 1:, j:] = (B[:, j + 1:, j:] - f[:, :, None] * B[:, j, None, j:]) % p
    return deficient


def screen_degree(x: np.ndarray, q: int, p: int, d: int, ms: np.ndarray) -> np.ndarray:
    """Return the m in `ms` for which pairs (x_k, x_{mk}), k = 1..2d+3, admit a (d,1) relation."""
    s = 2 * d + 3
    k = np.arange(1, s + 1, dtype=np.int64)
    xs = x[k - 1]
    V = monomials(xs, np.ones_like(xs), d, 0, p)               # (s, d+1): x^0..x^d
    Nl = np.array(nullspace_mod_p(V.T.copy(), p), dtype=np.int64)  # (d+2, s), Nl @ V = 0
    assert Nl.shape == (d + 2, s) and not np.any(Nl @ V % p)
    C = (Nl.T[:, :, None] * V[:, None, :]) % p                 # (s, d+2, d+1)
    hits = []
    chunk = max(256, int(3e7 // ((d + 2) * (d + 1) * 8)))
    for lo in range(0, len(ms), chunk):
        mm = ms[lo:lo + chunk]
        Y = x[(mm[:, None] * k[None, :]) % q - 1]               # (M, s)
        B = np.tensordot(Y, C, axes=(1, 0)) % p                # (M, d+2, d+1)
        hits.append(mm[batch_rank_deficient(B, p)])
    return np.concatenate(hits) if hits else np.array([], dtype=np.int64)


def reference_screen(x, q, p, d, m) -> bool:
    """Independent per-m decision: rank of the s x (2d+2) monomial matrix < 2d+2."""
    s = 2 * d + 3
    k = np.arange(1, s + 1, dtype=np.int64)
    X, Y = x[k - 1], x[(m * k) % q - 1]
    return len(nullspace_mod_p(monomials(X, Y, d, 1, p), p)) > 0


def verify(x, q, p, d, m, rng) -> np.ndarray | None:
    """Fit a (d,1) relation on random pairs; return it if it vanishes on ALL q-1 pairs."""
    k_all = np.arange(1, q, dtype=np.int64)
    X, Y = x, x[(m * k_all) % q - 1]
    idx = rng.choice(q - 1, size=min(q - 1, 4 * d + 30), replace=False)
    ns = nullspace_mod_p(monomials(X[idx], Y[idx], d, 1, p), p)
    for v in ns:
        if np.all(evaluate(v, X, Y, d, 1, p) == 0):
            return v
    return None


def search(x, q, p, dmax, rng, log) -> dict[int, dict[int, np.ndarray]]:
    """found[d] = {m: relation coefficients} for every verified m at degree d (cumulative)."""
    ms = np.arange(1, q, dtype=np.int64)
    found = {}
    for d in range(1, dmax + 1):
        t0 = time.perf_counter()
        cand = screen_degree(x, q, p, d, ms)
        rel = {}
        for m in cand.tolist():
            v = verify(x, q, p, d, m, rng)
            if v is not None:
                rel[m] = v
        found[d] = rel
        log.append({"d": d, "screen_candidates": int(len(cand)), "verified": len(rel),
                    "seconds": round(time.perf_counter() - t0, 2)})
    return found


# ------------------------------------------------------------------ analysis


def split_relation(v: np.ndarray, d: int, p: int):
    """Relation sum c_{i,j} X^i Y^j (j in {0,1}) -> R = num/den in lowest terms."""
    phi = [int(v[2 * i]) for i in range(d + 1)]
    psi = [int(v[2 * i + 1]) for i in range(d + 1)]
    num = [(-c) % p for c in phi]
    g = pgcd(num, psi, p)
    from src.cm import pdivmod
    num, den = pdivmod(num, g, p)[0], pdivmod(psi, g, p)[0]
    return num, den


def doubling_formula(a, b, p):
    """x(2P) = (x^4 - 2a x^2 - 8b x + a^2) / (4 (x^3 + a x + b))  (closed form)."""
    num = [a * a % p, (-8 * b) % p, (-2 * a) % p, 0, 1]
    den = [4 * b % p, 4 * a % p, 0, 4]
    return num, den


def analyse_curve(c: dict, x: np.ndarray, found, dmax: int, rng) -> dict:
    p, q, t, a, b = c["p"], c["q"], c["t"], c["a"], c["b"]
    f_poly = [b, a, 0, 1]
    # --- observed: minimal degree of every found multiplier
    mindeg = {}
    for d in range(1, dmax + 1):
        for m in found[d]:
            mindeg.setdefault(m, d)
    monotone = all(set(found[d]) <= set(found[d + 1]) for d in range(1, dmax))
    new = {d: sorted(m for m, dd in mindeg.items() if dd == d) for d in range(1, dmax + 1)}

    # --- orthogonal check of the j-table: CM field of Frobenius
    dk_frob, f_frob = fundamental_discriminant(t * t - 4 * p)
    res = {"curve": {k: v for k, v in c.items()}, "D_K_of_frobenius": dk_frob,
           "conductor_of_Z_pi": f_frob, "monotone_in_d": monotone,
           "observed_new_by_degree": {d: v for d, v in new.items() if v},
           "observed_counts": {d: len(v) for d, v in new.items()}}

    # --- predictions under competing hypotheses for End(E)
    hyps = {}
    if c["D_table"] is not None:
        hyps["End = O(table)"] = c["D_table"]
    hyps["End = O_K"] = dk_frob
    hyps["End = Z[pi]"] = t * t - 4 * p
    comp = {}
    for name, D in hyps.items():
        try:
            mu = eigenvalue_of_w(p, t, D, q)
        except (ValueError, AssertionError) as e:
            comp[name] = {"D": D, "applicable": False, "reason": str(e)}
            continue
        pred = {d: sorted(predicted_multipliers(D, d, mu, q)) for d in range(1, dmax + 1)}
        mism = {d: {"predicted": pred[d], "observed": new[d]}
                for d in range(1, dmax + 1) if pred[d] != new[d]}
        comp[name] = {"D": D, "applicable": True, "mu": mu,
                      "predicted_counts": {d: len(v) for d, v in pred.items()},
                      "r_O_check": all(len(pred[d]) == r_O(D, d) for d in pred),
                      "exact_set_match_all_d": not mism, "mismatches": mism}
    res["hypotheses"] = comp

    # --- blind inference of the CM structure from the found multipliers alone
    inferred = []
    for m, d in sorted(mindeg.items(), key=lambda kv: (kv[1], kv[0])):
        T = trace_for_multiplier(m, d, q)
        entry = {"m": m, "d": d, "trace": T}
        if T is not None and T * T - 4 * d < 0:
            entry["delta"] = T * T - 4 * d
            entry["D_K"], entry["f"] = fundamental_discriminant(entry["delta"])
        inferred.append(entry)
    nonint = [e for e in inferred if "D_K" in e]
    res["inference"] = {
        "all_have_trace": all(e["trace"] is not None for e in inferred),
        "n_integer_multipliers": sum(1 for e in inferred if e.get("trace") is not None
                                     and e["trace"] ** 2 == 4 * e["d"]),
        "n_nonintegral": len(nonint),
        "D_K_values": sorted({e["D_K"] for e in nonint}),
        "inferred_order_disc": (nonint[0]["D_K"] * math.gcd(*[e["f"] for e in nonint]) ** 2
                                if nonint else None),
    }

    # --- structure of each found map at its minimal degree
    struct = []
    for m, d in sorted(mindeg.items(), key=lambda kv: (kv[1], kv[0])):
        num, den = split_relation(found[d][m], d, p)
        ks = kernel_structure(den, f_poly, p) if pdeg(den) > 0 else {"deg_den": pdeg(den)}
        ok = (pdeg(num) == d and pdeg(den) == d - 1 and
              (d == 1 or (ks["exact_split"] and ks["h_squarefree"] and ks["s_h_coprime"]
                          and (ks["deg_s"] <= 0 or ks["s_divides_f"]))))
        struct.append({"m": m, "d": d, "deg_num": pdeg(num), "deg_den": pdeg(den),
                       "kernel": ks, "ok": ok})
    res["map_structure_all_ok"] = all(s["ok"] for s in struct)
    res["map_structure"] = struct
    if 2 in found.get(4, {}):
        num, den = split_relation(found[4][2], 4, p)
        n0, d0 = doubling_formula(a, b, p)
        res["doubling_formula_match"] = pmul(num, d0, p) == pmul(den, n0, p)
    return res


def audit_screen(x, q, p, d, rng, n_random=150) -> dict:
    """Compare the batched screen with the independent per-m reference on a random subset."""
    ms = np.arange(1, q, dtype=np.int64)
    batched = set(screen_degree(x, q, p, d, ms).tolist())
    sample = set(rng.choice(ms, size=min(n_random, q - 1), replace=False).tolist()) | batched
    disagreements = [m for m in sorted(sample) if (m in batched) != reference_screen(x, q, p, d, m)]
    return {"d": d, "checked": len(sample), "batched_positive": len(batched),
            "disagreements": disagreements}


def main():
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--dmax", type=int, default=20)
    ap.add_argument("--dmax163", type=int, default=42)
    ap.add_argument("--generic", type=int, default=4)
    ap.add_argument("--out", default=str(OUT), help="output directory (created if missing)")
    args = ap.parse_args()
    from pathlib import Path
    OUT = Path(args.out)
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(MASTER_SEED)
    curves = []
    for i, D in enumerate(sorted(CLASS_NUMBER_ONE, reverse=True)):
        curves.append((build_cm_curve(D, MASTER_SEED + 100 * i), args.dmax163 if D == -163 else args.dmax))
    for g in range(args.generic):
        curves.append((build_generic_curve(MASTER_SEED + 5000 + 100 * g, args.dmax), args.dmax))
    results, timing = [], {}
    for idx, (c, dmax) in enumerate(curves):
        t0 = time.perf_counter()
        x = subgroup_sequence(c, MASTER_SEED + idx)
        log = []
        found = search(x, c["q"], c["p"], dmax, rng, log)
        r = analyse_curve(c, x, found, dmax, rng)
        r["dmax"] = dmax
        r["search_log"] = log
        if idx in (0, 1, 3):  # D = -3, -4, -8: audit the batched screen against the reference
            r["screen_audit"] = [audit_screen(x, c["q"], c["p"], d, rng) for d in (1, 2, 3, 5, 8)]
        name = c["family"] if c["D_table"] is None else f"D={c['D_table']}"
        timing[f"{name}#{idx}"] = round(time.perf_counter() - t0, 2)
        results.append(r)
        h = r["hypotheses"]
        print(f"{name:>8} p={c['p']} q={c['q']} t={c['t']}: counts "
              f"{ {d: n for d, n in r['observed_counts'].items() if n} }  "
              + "  ".join(f"{k}: {v.get('exact_set_match_all_d')}" for k, v in h.items())
              + f"  [{timing[f'{name}#{idx}']}s]", flush=True)
    summary = {"provenance": provenance(MASTER_SEED), "bits": BITS, "max_cofactor": MAX_COFACTOR,
               "dmax": args.dmax, "dmax163": args.dmax163, "timing_s": timing, "curves": results}
    save_json(summary, OUT / "summary.json")
    write_markdown(summary)
    figure(summary)


def figure(s):
    """Observed number of multipliers of minimal degree d, per curve; cells that disagree with the
    table-order prediction (CM curves) or the integer-only prediction (generic) are outlined."""
    from common import SERIES, TEXT, TEXT2, setup_matplotlib
    from matplotlib.colors import LinearSegmentedColormap
    plt = setup_matplotlib()
    dmax = s["dmax"]
    rows, labels, bad = [], [], []
    for r in s["curves"]:
        c = r["curve"]
        key = "End = O(table)" if c["D_table"] is not None else "End = O_K"
        h = r["hypotheses"][key]
        obs = [r["observed_counts"].get(str(d), r["observed_counts"].get(d, 0)) for d in range(1, dmax + 1)]
        rows.append(obs)
        labels.append(f"D={c['D_table']}" if c["D_table"] is not None else f"generic (D_K={r['D_K_of_frobenius']})")
        mism = {int(k) for k in h.get("mismatches", {})}
        bad.append([d in mism for d in range(1, dmax + 1)])
    A = np.array(rows, dtype=float)
    cmap = LinearSegmentedColormap.from_list("seq", ["#fcfcfb", SERIES[0]])
    fig, ax = plt.subplots(figsize=(11, 0.42 * len(rows) + 1.6))
    ax.imshow(np.sqrt(A), cmap=cmap, aspect="auto")
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if A[i, j]:
                ax.text(j, i, int(A[i, j]), ha="center", va="center", fontsize=7,
                        color=TEXT if A[i, j] < 9 else "#ffffff")
            if bad[i][j]:
                ax.add_patch(__import__("matplotlib").patches.Rectangle((j - .5, i - .5), 1, 1,
                             fill=False, ec="#e34948", lw=2))
    ax.set_xticks(range(dmax), [str(d) for d in range(1, dmax + 1)])
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ax.set_xlabel("degree d of the rational map x(mP) = R(x(P))")
    ax.grid(False)
    nbad = int(np.sum(bad))
    ax.set_title(f"Number of multipliers m of minimal degree d (blind search over all m).\n"
                 f"Predicted value r_O(d) from CM theory: {nbad} mismatching cells (outlined) "
                 f"out of {A.size}", fontsize=9)
    fig.savefig(OUT / "C1_multipliers_by_degree.png")
    plt.close(fig)


def write_markdown(s):
    L = ["# Phase 2 summary (auto-generated)", "",
         f"p in [2^{s['bits'] - 1}, 2^{s['bits']}), cofactor <= {s['max_cofactor']}, "
         f"dmax = {s['dmax']} (D=-163: {s['dmax163']}). Seed {s['provenance']['seed']}.", "",
         "## Observed multipliers by minimal degree vs predictions", "",
         "| curve | p | q | t | D_K(t^2-4p) | observed #new at d (nonzero only) | End=O(table) | End=O_K | End=Z[pi] | inferred disc | maps ok | [2] formula |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in s["curves"]:
        c = r["curve"]
        name = c["family"]
        h = r["hypotheses"]

        def hv(k):
            if k not in h:
                return "n/a"
            return "not applicable" if not h[k]["applicable"] else ("match" if h[k]["exact_set_match_all_d"] else "MISMATCH")
        L.append(f"| {name} | {c['p']} | {c['q']} | {c['t']} | {r['D_K_of_frobenius']} | "
                 f"{ {d: n for d, n in r['observed_counts'].items() if n} } | {hv('End = O(table)')} | "
                 f"{hv('End = O_K')} | {hv('End = Z[pi]')} | {r['inference']['inferred_order_disc']} | "
                 f"{r['map_structure_all_ok']} | {r.get('doubling_formula_match', '-')} |")
    L += ["", "## Screen audit (batched vs independent per-m reference)", ""]
    for r in s["curves"]:
        if "screen_audit" in r:
            L.append(f"- {r['curve']['family']}: " + "; ".join(
                f"d={a['d']}: {a['checked']} m checked, {len(a['disagreements'])} disagreements"
                for a in r["screen_audit"]))
    L += ["", "## Screening statistics", ""]
    tot_c = sum(e["screen_candidates"] for r in s["curves"] for e in r["search_log"])
    tot_v = sum(e["verified"] for r in s["curves"] for e in r["search_log"])
    L.append(f"Screen candidates {tot_c}, verified {tot_v} (unverified screen positives: {tot_c - tot_v}).")
    (OUT / "summary.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
