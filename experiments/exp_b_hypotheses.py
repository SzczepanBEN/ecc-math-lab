"""Experiment B -- testing and classifying the structure seen in Experiment A.

  H1  Translation relation: (x_k, x_{k+l}) always satisfy one fixed polynomial of
      bidegree (2,2) -- predicted to be Semaev's 3rd summation polynomial S3(X, Y, x(lG)) --
      dropping to bidegree (1,1) exactly when lG is 2-torsion.
  H2  Index-scaling symmetries: the multipliers m for which x_{mk} = M(x_k) for one Mobius map
      M are exactly the automorphisms (|Aut E| of them: 2 generic, 6 for j=0, 4 for j=1728).
  H3  No low-complexity description of k -> x_k: no polynomial F(k, x_k) = 0 of small
      bidegree, and linear complexity ~ N/2 as for a random sequence. Positive controls included.
  H4  "Consecutive x's are predictable" by a plug-in conditional entropy: expected to be a
      statistical artifact; contrasted with the exact algebraic test from H1.
  H5  Validation of the nominal Experiment-A hits (J0-B diff chi^2; GEN-B x-spread) using
      bigger nulls, other generators and a curve-level null.
Outputs go to results/exp_b/.

Usage:  python experiments/exp_b_hypotheses.py
"""
from __future__ import annotations

import math
import random

import numpy as np

from common import (BASELINE, RESULTS, SERIES, TEXT2, curve_set, save_json, sequence,
                    setup_matplotlib)
from exp_a_geometry import chi2_counts, mc_pvalue, ordering_stats
from src.curvegen import generate_curve
from src.ec import Curve
from src.harness import provenance, timed

OUT = RESULTS / "exp_b"
SEED = 2


# ------------------------------------------------------------------ linear algebra mod p

def nullspace_mod_p(A: np.ndarray, p: int) -> list[np.ndarray]:
    """Basis of {v : A v = 0 mod p}. Entries must be < p < 2^20 (int64 products stay exact)."""
    A = (A % p).astype(np.int64)
    r, c = A.shape
    pivots, row = [], 0
    for col in range(c):
        if row == r:
            break
        nz = np.flatnonzero(A[row:, col])
        if nz.size == 0:
            continue
        pr = row + nz[0]
        A[[row, pr]] = A[[pr, row]]
        A[row] = A[row] * pow(int(A[row, col]), -1, p) % p
        f = A[:, col].copy()
        f[row] = 0
        A = (A - np.outer(f, A[row]) % p) % p
        pivots.append(col)
        row += 1
    basis = []
    for fc in (j for j in range(c) if j not in pivots):
        v = np.zeros(c, dtype=np.int64)
        v[fc] = 1
        for i, pc in enumerate(pivots):
            v[pc] = (-A[i, fc]) % p
        basis.append(v)
    return basis


def monomials(X: np.ndarray, Y: np.ndarray, d1: int, d2: int, p: int) -> np.ndarray:
    """Columns X^i Y^j (0<=i<=d1, 0<=j<=d2) mod p, ordered i-major."""
    px = [np.ones_like(X)]
    for _ in range(d1):
        px.append(px[-1] * X % p)
    py = [np.ones_like(Y)]
    for _ in range(d2):
        py.append(py[-1] * Y % p)
    return np.stack([px[i] * py[j] % p for i in range(d1 + 1) for j in range(d2 + 1)], axis=1)


def normalize(v: np.ndarray, p: int) -> np.ndarray:
    nz = np.flatnonzero(v)
    return v * pow(int(v[nz[-1]]), -1, p) % p


def evaluate(coeffs: np.ndarray, X, Y, d1, d2, p) -> np.ndarray:
    M = monomials(X, Y, d1, d2, p)
    return (M * coeffs[None, :] % p).sum(axis=1) % p


def relation_search(X, Y, p, bidegrees, rng, samples=400):
    """For each bidegree: dimension of the space of F with F(X_i, Y_i) = 0 on a random sample,
    and for a 1-dim space, the fraction of ALL pairs on which that F vanishes."""
    idx = rng.choice(len(X), size=min(samples, len(X)), replace=False)
    out = {}
    for d1, d2 in bidegrees:
        ns = nullspace_mod_p(monomials(X[idx], Y[idx], d1, d2, p), p)
        rec = {"nullspace_dim": len(ns), "n_monomials": (d1 + 1) * (d2 + 1)}
        if len(ns) == 1:
            F = normalize(ns[0], p)
            rec["coeffs"] = F.tolist()
            rec["frac_all_pairs_vanishing"] = float(np.mean(evaluate(F, X, Y, d1, d2, p) == 0))
        out[f"({d1},{d2})"] = rec
    return out


def semaev_s3_coeffs(x3: int, a: int, b: int, p: int) -> np.ndarray:
    """S3(X, Y, x3) = (X-Y)^2 x3^2 - 2((X+Y)(XY+a) + 2b) x3 + (XY-a)^2 - 4b(X+Y),
    as coefficients of X^i Y^j, i-major over 0..2 x 0..2 (Semaev 2004)."""
    c = np.zeros((3, 3), dtype=object)
    c[2, 2] = 1
    c[2, 1] = c[1, 2] = -2 * x3
    c[2, 0] = c[0, 2] = x3 * x3
    c[1, 1] = -2 * x3 * x3 - 2 * a
    c[1, 0] = c[0, 1] = -2 * a * x3 - 4 * b
    c[0, 0] = a * a - 4 * b * x3
    return np.array([int(v) % p for v in c.flatten()], dtype=np.int64)


# ------------------------------------------------------------------ H1


def h1(rng) -> dict:
    res = {}
    bideg = [(1, 1), (1, 2), (2, 1), (2, 2)]
    for rec in curve_set():
        p, a, b, n = rec.p, rec.a, rec.b, rec.order
        x, _ = sequence(rec.name)
        lags = [1, 2, 17, 1000, (n - 1) // 3]
        if n % 2 == 0:
            lags.append(n // 2)
        if n % 3 == 0:
            lags.append(n // 3)
        per = {}
        for l in lags:
            # pairs (x_k, x_{k+l}) for k, k+l in 1..n-1  (array index i <-> k = i+1)
            X, Y = x[: n - 1 - l], x[l:]
            r = relation_search(X, Y, p, bideg, rng)
            xq = int(x[l - 1])
            s3 = normalize(semaev_s3_coeffs(xq, a, b, p), p)
            r22 = r["(2,2)"]
            entry = {"search": r, "x(lG)": xq}
            if "coeffs" in r22:
                entry["matches_semaev_S3"] = bool(np.array_equal(np.array(r22["coeffs"]), s3))
            if r["(1,1)"]["nullspace_dim"] == 1:
                # predicted Mobius relation for a 2-torsion translate T=(t,0): (Y-t)(X-t) = 3t^2+a
                t = xq
                pred = normalize(np.array([(t * t - (3 * t * t + a)) % p, (-t) % p, (-t) % p, 1],
                                          dtype=np.int64), p)
                entry["mobius_matches_(Y-t)(X-t)=3t^2+a"] = bool(np.array_equal(
                    np.array(r["(1,1)"]["coeffs"]), pred))
            # baseline: random permutation of the x values
            Xs, Ys = rng.permutation(X), rng.permutation(Y)
            base = relation_search(Xs, Ys, p, [(1, 1), (2, 2)], rng)
            entry["baseline_nullspace_dims"] = {k: v["nullspace_dim"] for k, v in base.items()}
            entry["baseline_frac_S3_vanishing"] = float(np.mean(evaluate(s3, Xs, Ys, 2, 2, p) == 0))
            per[str(l)] = entry
        res[rec.name] = per
    return res


# ------------------------------------------------------------------ H2


def mobius_multipliers(x: np.ndarray, n: int, p: int) -> list[int]:
    """All m in [1, n-1] with x_{mk} = M(x_k) for a single Mobius M (blind search).

    Screen: Mobius maps preserve cross-ratios, so CR(x_1..x_4) = CR(x_m..x_4m) and
    CR(x_2..x_5) = CR(x_2m..x_5m) (cross-multiplied, no inverses). Survivors are then fitted
    with a (1,1) relation and verified on every k.
    """
    def xs(j):  # x_j for an integer array j (mod n); -1 marks the identity
        j = j % n
        out = np.full(j.shape, -1, dtype=np.int64)
        ok = j != 0
        out[ok] = x[j[ok] - 1]
        return out

    m = np.arange(1, n, dtype=np.int64)
    cols = [xs(m * t) for t in range(1, 6)]
    ref = [int(x[t - 1]) for t in range(1, 6)]
    alive = np.all(np.stack(cols) >= 0, axis=0)

    def cr_parts(a, b, c, d):
        return ((a - c) % p) * ((b - d) % p) % p, ((a - d) % p) * ((b - c) % p) % p

    for (i, j, k, l) in [(0, 1, 2, 3), (1, 2, 3, 4)]:
        rn, rd = cr_parts(ref[i], ref[j], ref[k], ref[l])
        cn, cd = cr_parts(cols[i], cols[j], cols[k], cols[l])
        alive &= (cd != 0) & ((cn * rd - rn * cd) % p == 0)
    found = []
    k = np.arange(1, n, dtype=np.int64)
    for mm in m[alive]:
        X, Y = x, xs(mm * k)
        okk = Y >= 0
        ns = nullspace_mod_p(monomials(X[okk][:60], Y[okk][:60], 1, 1, p), p)
        if len(ns) != 1:
            continue
        if np.all(evaluate(ns[0], X[okk], Y[okk], 1, 1, p) == 0) and okk.sum() == n - 1:
            found.append(int(mm))
    return found


def h2() -> dict:
    res = {"main_curves": {}, "extra_curves": []}
    for rec in curve_set():
        x, _ = sequence(rec.name)
        ms = mobius_multipliers(x, rec.order, rec.p)
        res["main_curves"][rec.name] = {
            "family": rec.family, "n": rec.order, "multipliers": ms, "count": len(ms),
            "aut_group_order": rec.properties["aut_group_order"],
            "match": len(ms) == rec.properties["aut_group_order"]}
    # generalisation: 24 further cyclic curves with 14-bit p, 8 per family
    specs = [("generic", None, True), ("j0", (1, 3), True), ("j1728", (1, 4), False)]
    for family, mod, prime in specs:
        seed, got = 500, 0
        while got < 8:
            seed += 1
            r = generate_curve(seed, 14, family, require_prime_order=prime, p_mod=mod)
            if r.structure[1] != 1:
                continue
            E = r.curve()
            pts = E.multiples(r.generator, r.order - 1)
            x = np.array([P[0] for P in pts], dtype=np.int64)
            ms = mobius_multipliers(x, r.order, r.p)
            res["extra_curves"].append({
                "family": family, "seed": seed, "p": r.p, "a": r.a, "b": r.b, "n": r.order,
                "count": len(ms), "aut_group_order": r.properties["aut_group_order"],
                "match": len(ms) == r.properties["aut_group_order"]})
            got += 1
    res["all_match"] = all(v["match"] for v in res["main_curves"].values()) and \
        all(v["match"] for v in res["extra_curves"])
    return res


# ------------------------------------------------------------------ H3


def berlekamp_massey_profile(s: list[int], p: int) -> list[int]:
    """Linear complexity L_N of s_0..s_{N-1} over F_p for every N (Berlekamp-Massey)."""
    C, B = [1], [1]
    L, m, bb = 0, 1, 1
    prof = []
    for i, si in enumerate(s):
        d = si
        for j in range(1, L + 1):
            d = (d + C[j] * s[i - j]) % p
        if d == 0:
            m += 1
        else:
            coef = d * pow(bb, -1, p) % p
            T = C[:]
            C = C + [0] * (len(B) + m - len(C))
            for j, bj in enumerate(B):
                C[j + m] = (C[j + m] - coef * bj) % p
            if 2 * L <= i:
                L, B, bb, m = i + 1 - L, T, d, 1
            else:
                m += 1
        prof.append(L)
    return prof


def h3(rng) -> dict:
    res = {"polynomial_index_relation": {}, "linear_complexity": {}}
    bideg = [(1, 1), (2, 2), (3, 1), (4, 4), (8, 2), (16, 2), (32, 1), (10, 10), (60, 1)]
    for rec in curve_set():
        p, n = rec.p, rec.order
        x, y = sequence(rec.name)
        K = min(n - 1, p - 1)
        k = np.arange(1, K + 1, dtype=np.int64)
        xk = x[:K]
        entry = {}
        for d1, d2 in bideg:
            M = (d1 + 1) * (d2 + 1)
            idx = rng.choice(K, size=2 * M + 40, replace=False)
            dim_real = len(nullspace_mod_p(monomials(k[idx], xk[idx], d1, d2, p), p))
            dim_perm = len(nullspace_mod_p(monomials(k[idx], rng.permutation(xk)[idx], d1, d2, p), p))
            planted = (k ** 3 % p * k % p + 5 * k + 1) % p  # positive control: Y = k^4 + 5k + 1
            dim_ctrl = len(nullspace_mod_p(monomials(k[idx], planted[idx], d1, d2, p), p))
            entry[f"({d1},{d2})"] = {"monomials": M, "rows": 2 * M + 40, "dim_real": dim_real,
                                     "dim_permuted": dim_perm, "dim_planted_control": dim_ctrl}
        res["polynomial_index_relation"][rec.name] = entry
        if rec.name in ("GEN-A", "J0-A", "J0-SS"):
            N = 800
            res["linear_complexity"][rec.name] = {
                "x_profile": berlekamp_massey_profile([int(v) for v in x[:N]], p),
                "y_profile": berlekamp_massey_profile([int(v) for v in y[:N]], p),
                "perm_profile": berlekamp_massey_profile([int(v) for v in rng.permutation(x)[:N]], p),
            }
    p = curve_set()[0].p
    res["linear_complexity"]["control_k^4+5k+1"] = berlekamp_massey_profile(
        [(kk**4 + 5 * kk + 1) % p for kk in range(1, 801)], p)
    return res


# ------------------------------------------------------------------ H4


def cond_entropy_bits(a: np.ndarray, b: np.ndarray) -> float:
    """Plug-in H(B | A) in bits from paired samples."""
    pairs = a.astype(np.int64) * (1 << 21) + b
    _, cab = np.unique(pairs, return_counts=True)
    _, ca = np.unique(a, return_counts=True)
    N = len(a)
    hab = -np.sum(cab / N * np.log2(cab / N))
    ha = -np.sum(ca / N * np.log2(ca / N))
    return float(hab - ha)


def h4(rng) -> dict:
    res = {}
    for rec in curve_set():
        p, n = rec.p, rec.order
        x, _ = sequence(rec.name)
        xs = rng.permutation(x)
        H = (n - 1) // 2
        xh = rng.permutation(x[:H])
        s3 = semaev_s3_coeffs(int(x[0]), rec.a, rec.b, p)
        res[rec.name] = {
            "log2_p": math.log2(p),
            "plugin_H_next_given_curr_full_period": cond_entropy_bits(x[:-1], x[1:]),
            "plugin_H_next_given_curr_full_period_PERMUTED": cond_entropy_bits(xs[:-1], xs[1:]),
            "plugin_H_next_given_curr_half_period": cond_entropy_bits(x[: H - 1], x[1:H]),
            "plugin_H_next_given_curr_half_period_PERMUTED": cond_entropy_bits(xh[:-1], xh[1:]),
            "frac_S3_vanishing_real": float(np.mean(evaluate(s3, x[:-1], x[1:], 2, 2, p) == 0)),
            "frac_S3_vanishing_PERMUTED": float(np.mean(evaluate(s3, xs[:-1], xs[1:], 2, 2, p) == 0)),
            # the two roots of S3(x_k, Y, x_G) are curve x-values, each occurring twice in the
            # shuffled multiset of ~p values, so the chance rate is ~4/p (not 2/p)
            "expected_S3_vanishing_if_shuffled_~4/p": 4 / p,
        }
    return res


# ------------------------------------------------------------------ H5


def x_spread_chi2_from_curve(p, a, b, table) -> float:
    xx = np.arange(p, dtype=np.int64)
    occ = table[((xx * xx % p) * xx % p + a * xx + b) % p]
    return chi2_counts(np.bincount(xx * 64 // p, weights=occ, minlength=64))


def diff_chi2(x: np.ndarray, p: int) -> float:
    """Same statistic as ordering_stats()['diff_chi2'] (64-bin chi^2 of x_{k+1} - x_k mod p)."""
    return chi2_counts(np.bincount((x[1:] - x[:-1]) % p * 64 // p, minlength=64))


def h5(rng) -> dict:
    res = {}
    recs = {r.name: r for r in curve_set()}
    # (a) J0-B consecutive-difference chi^2: bigger null + 30 other generators
    r = recs["J0-B"]
    x, _ = sequence(r.name)
    E = r.curve()
    H = (r.order - 1) // 2
    obs = ordering_stats(x[:H], r.p)["diff_chi2"]
    assert abs(diff_chi2(x[:H], r.p) - obs) < 1e-9
    null = np.array([diff_chi2(rng.permutation(x[:H]), r.p) for _ in range(2000)])
    pr = random.Random(SEED)
    gen_p = []
    for _ in range(30):
        c = pr.randrange(2, r.order - 1)
        G2 = E.mul(c, r.generator)
        xs2 = np.array([P[0] for P in E.multiples(G2, H)], dtype=np.int64)
        gen_p.append(mc_pvalue(diff_chi2(xs2, r.p), null))
    res["J0-B_diff_chi2"] = {
        "observed": obs, "p_mc_2000": mc_pvalue(obs, null),
        "other_generators_p_values": gen_p,
        "other_generators_frac_below_0.05": float(np.mean(np.array(gen_p) < 0.05)),
    }
    # (b) GEN-B x-spread: curve-level null (3000 random curves over the same p) vs relabelling null.
    # A first run used 300 curves (seed 3) and showed a spurious KS difference; see research_log.md
    # and experiments/validate_curve_level_null.py.
    r = recs["GEN-B"]
    p = r.p
    E = r.curve()
    table = E._square_table
    obs = x_spread_chi2_from_curve(p, r.a, r.b, table)
    cr = random.Random(SEED + 1)
    curve_null = []
    while len(curve_null) < 3000:
        a, b = cr.randrange(1, p), cr.randrange(1, p)
        if (4 * a**3 + 27 * b * b) % p:
            curve_null.append(x_spread_chi2_from_curve(p, a, b, table))
    curve_null = np.array(curve_null)
    xx = np.arange(p, dtype=np.int64)
    occ = table[((xx * xx % p) * xx % p + r.a * xx + r.b) % p]
    bins = xx * 64 // p
    relabel_null = np.array([chi2_counts(np.bincount(bins[rng.permutation(p)], weights=occ,
                                                     minlength=64)) for _ in range(2000)])
    from scipy import stats
    res["GEN-B_x_spread"] = {
        "observed": obs, "p_relabel_2000": mc_pvalue(obs, relabel_null),
        "p_curve_level_3000": mc_pvalue(obs, curve_null),
        "curve_null_mean": float(curve_null.mean()), "relabel_null_mean": float(relabel_null.mean()),
        "curve_null_var": float(curve_null.var()), "relabel_null_var": float(relabel_null.var()),
        "KS_curve_vs_relabel_pvalue": float(stats.ks_2samp(curve_null, relabel_null).pvalue),
        "frac_random_curves_with_p_below_0.05_vs_relabel":
            float(np.mean([mc_pvalue(c, relabel_null) < 0.05 for c in curve_null])),
    }
    res["_nulls_for_plot"] = {"curve": curve_null.tolist(), "relabel": relabel_null.tolist(),
                              "obs": obs}
    return res


# ------------------------------------------------------------------ figures


def figures(plt, R):
    # B1: H3 linear complexity
    lc = R["H3"]["linear_complexity"]
    fig, ax = plt.subplots(figsize=(7, 3.8))
    N = np.arange(1, 801)
    ax.plot(N, N / 2, color=BASELINE, lw=1, ls="--", label="N/2 (expected for random)")
    for name, c in zip(["GEN-A", "J0-A", "J0-SS"], SERIES):
        ax.plot(N, lc[name]["x_profile"], color=c, lw=1.4, label=f"{name}: x_k")
    ax.plot(N, lc["GEN-A"]["perm_profile"], color=BASELINE, lw=1, label="GEN-A permuted")
    ax.plot(N, lc["control_k^4+5k+1"], color=TEXT2, lw=2, label="control k^4+5k+1 (L=5)")
    ax.set_xlabel("prefix length N")
    ax.set_ylabel("linear complexity L_N over F_p")
    ax.set_title("H3: linear complexity of x(kG) -- all curves and the permuted baseline\n"
                 "coincide exactly with ceil(N/2); only the planted control is low")
    ax.legend(fontsize=8, loc="upper left")
    fig.savefig(OUT / "B1_linear_complexity.png")
    plt.close(fig)

    # B2: H5 curve-level vs relabelling null for GEN-B x-spread
    d = R["H5"]["_nulls_for_plot"]
    fig, ax = plt.subplots(figsize=(7, 3.6))
    bins = np.linspace(min(d["curve"] + d["relabel"]), max(d["curve"] + d["relabel"] + [d["obs"]]), 40)
    ax.hist(d["relabel"], bins=bins, density=True, color=BASELINE, alpha=0.55,
            label="relabelling null (2000)")
    ax.hist(d["curve"], bins=bins, density=True, histtype="step", color=SERIES[0], lw=2,
            label="3000 random curves, same p")
    ax.axvline(d["obs"], color=SERIES[1], lw=2, label="GEN-B observed")
    ax.set_xlabel("chi^2 of x-occupancy over 64 bins")
    ax.set_ylabel("density")
    ax.set_title("H5: is GEN-B's x-spread unusual among curves over the same field?")
    ax.legend(fontsize=8)
    fig.savefig(OUT / "B2_curve_level_null.png")
    plt.close(fig)

    # B3: H4 entropy artifact vs algebraic test
    h = R["H4"]
    names = list(h)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    w = 0.38
    i = np.arange(len(names))
    axes[0].bar(i - w / 2, [h[k]["plugin_H_next_given_curr_full_period"] for k in names], w * 0.92,
                color=SERIES[0], label="real sequence")
    axes[0].bar(i + w / 2, [h[k]["plugin_H_next_given_curr_full_period_PERMUTED"] for k in names],
                w * 0.92, color=BASELINE, label="random permutation")
    axes[0].set_ylabel("plug-in H(x_k+1 | x_k)  [bits]")
    axes[0].set_title("Plug-in conditional entropy: no difference")
    axes[1].bar(i - w / 2, [h[k]["frac_S3_vanishing_real"] for k in names], w * 0.92,
                color=SERIES[0], label="real sequence")
    axes[1].bar(i + w / 2, [h[k]["frac_S3_vanishing_PERMUTED"] for k in names], w * 0.92,
                color=BASELINE, label="random permutation")
    axes[1].set_ylabel("fraction of k with S3(x_k, x_k+1, x_G) = 0")
    axes[1].set_title("Algebraic test: total separation")
    for ax in axes:
        ax.set_xticks(i, names, rotation=30, fontsize=8)
    fig.legend(*axes[1].get_legend_handles_labels(), loc="lower center", ncol=2, fontsize=8,
               bbox_to_anchor=(0.5, -0.06))
    fig.savefig(OUT / "B3_entropy_artifact.png")
    plt.close(fig)


# ------------------------------------------------------------------ report


def write_markdown(R):
    L = ["# Experiment B summary (auto-generated)", ""]
    L += ["## H1 translation relations (nullspace dimension on 400 sampled pairs)", "",
          "| curve | lag l | (1,1) | (1,2) | (2,1) | (2,2) | (2,2) vanishes on all pairs | = Semaev S3 | Mobius = predicted | baseline dims (1,1)/(2,2) |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for c, per in R["H1"].items():
        for l, e in per.items():
            s = e["search"]
            L.append(f"| {c} | {l} | " + " | ".join(str(s[b]["nullspace_dim"]) for b in s) +
                     f" | {s['(2,2)'].get('frac_all_pairs_vanishing', '-')} | {e.get('matches_semaev_S3', '-')} | "
                     f"{e.get('mobius_matches_(Y-t)(X-t)=3t^2+a', '-')} | "
                     f"{e['baseline_nullspace_dims']['(1,1)']}/{e['baseline_nullspace_dims']['(2,2)']} |")
    L += ["", "## H2 Mobius multipliers vs |Aut(E)|", "",
          "| curve | family | n | multipliers found | |Aut| | match |", "|---|---|---|---|---|---|"]
    for c, e in R["H2"]["main_curves"].items():
        L.append(f"| {c} | {e['family']} | {e['n']} | {e['multipliers']} | {e['aut_group_order']} | {e['match']} |")
    ex = R["H2"]["extra_curves"]
    L += ["", f"Extra 14-bit curves: {sum(e['match'] for e in ex)}/{len(ex)} match "
          f"(counts: {[(e['family'], e['count']) for e in ex]}). All match: {R['H2']['all_match']}.", ""]
    L += ["## H3 polynomial relations F(k, x_k) = 0 (nullspace dims: real / permuted / planted control)", "",
          "| curve | " + " | ".join(next(iter(R["H3"]["polynomial_index_relation"].values()))) + " |",
          "|---|" + "---|" * len(next(iter(R["H3"]["polynomial_index_relation"].values())))]
    for c, e in R["H3"]["polynomial_index_relation"].items():
        L.append(f"| {c} | " + " | ".join(f"{v['dim_real']}/{v['dim_permuted']}/{v['dim_planted_control']}"
                                         for v in e.values()) + " |")
    L += ["", "Linear complexity at N = 800 (random expectation ~400):", ""]
    for c, e in R["H3"]["linear_complexity"].items():
        if isinstance(e, dict):
            L.append(f"- {c}: x {e['x_profile'][-1]}, y {e['y_profile'][-1]}, permuted {e['perm_profile'][-1]}")
        else:
            L.append(f"- {c}: {e[-1]}")
    L += ["", "## H4 entropy artifact", "", "| curve | plug-in H full | permuted | plug-in H half | permuted | S3 vanishing real | permuted |",
          "|---|---|---|---|---|---|---|"]
    for c, e in R["H4"].items():
        L.append(f"| {c} | {e['plugin_H_next_given_curr_full_period']:.4f} | {e['plugin_H_next_given_curr_full_period_PERMUTED']:.4f} | "
                 f"{e['plugin_H_next_given_curr_half_period']:.4f} | {e['plugin_H_next_given_curr_half_period_PERMUTED']:.4f} | "
                 f"{e['frac_S3_vanishing_real']:.6f} | {e['frac_S3_vanishing_PERMUTED']:.6f} |")
    h5 = {k: v for k, v in R["H5"].items() if not k.startswith("_")}
    L += ["", "## H5 validation of nominal Experiment-A hits", "", "```", str(h5), "```"]
    (OUT / "summary.md").write_text("\n".join(L) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    plt = setup_matplotlib()
    timing = {}
    R = {}
    with timed("H1", timing):
        R["H1"] = h1(rng)
    with timed("H2", timing):
        R["H2"] = h2()
    with timed("H3", timing):
        R["H3"] = h3(rng)
    with timed("H4", timing):
        R["H4"] = h4(rng)
    with timed("H5", timing):
        R["H5"] = h5(rng)
    with timed("figures", timing):
        figures(plt, R)
    R["provenance"] = provenance(SEED)
    R["timing_s"] = timing
    save_json(R, OUT / "summary.json")
    write_markdown(R)


if __name__ == "__main__":
    main()
