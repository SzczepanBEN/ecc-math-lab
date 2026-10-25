"""Experiment A -- geometry of the scalar-multiplication sequence G, 2G, ..., (n-1)G.

For every curve in experiments/common.py this script measures
  A1  value multiplicities of x, y (and x^3 for j=0) over the full period, vs theory;
  A2  spread of the occupied coordinate values in [0, p), vs random relabelling of F_p;
  A3  exact negation symmetry x_{n-k} = x_k, y_{n-k} = -y_k;
  A4  ordering statistics on the half period k = 1..(n-1)/2 (autocorrelation, Ljung-Box,
      consecutive-difference histogram, 2-D serial test, rank correlation with k,
      initial-segment discrepancy, periodogram peak) vs random-permutation baselines;
  A5  positions of small-order points (group structure) and automorphism checks.
Outputs go to results/exp_a/.

Usage:  python experiments/exp_a_geometry.py [--perms 200]
"""
from __future__ import annotations

import argparse
import math
from collections import Counter

import numpy as np
from scipy import stats

from common import (BASELINE, RESULTS, SERIES, TEXT2, curve_set, save_json, sequence,
                    setup_matplotlib, write_curve_table)
from src.curvegen import generate_curve
from src.ec import sqrt_mod
from src.harness import provenance, timed

OUT = RESULTS / "exp_a"
SEED = 1
LAGS = 100

# ------------------------------------------------------------------ statistics


def acf(u: np.ndarray, lags: int) -> np.ndarray:
    v = u - u.mean()
    n = len(v)
    f = np.fft.rfft(v, 2 * n)
    r = np.fft.irfft(f * np.conj(f))[: lags + 1]
    return r[1:] / r[0]


def chi2_counts(c: np.ndarray) -> float:
    e = c.sum() / c.size
    return float(((c - e) ** 2 / e).sum())


def ordering_stats(x: np.ndarray, p: int) -> dict:
    """Statistics that depend on the ORDER of the sequence (multiset held fixed)."""
    u = x / p
    H = len(x)
    r = acf(u, LAGS)
    d = (x[1:] - x[:-1]) % p
    nb = 64
    diff_counts = np.bincount(d * nb // p, minlength=nb)
    g = 16
    cells = (x[:-1] * g // p) * g + (x[1:] * g // p)
    serial_counts = np.bincount(cells, minlength=g * g)
    rank_corr = np.corrcoef(np.arange(H), np.argsort(np.argsort(x)))[0, 1]
    pg = np.abs(np.fft.rfft(u - u.mean()))[1:] ** 2
    return {
        "acf_max_abs": float(np.abs(r).max()),
        "ljung_box_Q": float(H * (H + 2) * np.sum(r**2 / (H - np.arange(1, LAGS + 1)))),
        "diff_chi2": chi2_counts(diff_counts),
        "serial2d_chi2": chi2_counts(serial_counts),
        "rank_corr_k_abs": float(abs(rank_corr)),
        "ks_first_1000": float(stats.kstest(u[:1000], "uniform").statistic),
        "periodogram_peak_ratio": float(pg.max() / pg.mean()),
    }


def mc_pvalue(obs: float, null: np.ndarray) -> float:
    """Upper-tail Monte-Carlo p-value with the +1 correction."""
    return float((1 + np.sum(null >= obs)) / (1 + len(null)))


def multiplicity_profile(values: np.ndarray, p: int) -> dict[int, int]:
    """#{v in F_p : v occurs exactly m times}, for each m."""
    occ = np.bincount(values, minlength=p)
    return {int(m): int(c) for m, c in sorted(Counter(occ.tolist()).items())}


def _symmetric_relabelling(p: int, rng: np.random.Generator) -> np.ndarray:
    """Random bijection of F_p commuting with v -> -v: permutes the pairs {v, -v} and flips signs."""
    reps = np.arange(1, (p - 1) // 2 + 1)
    img = rng.permutation(reps) * rng.choice([1, -1], size=reps.size) % p
    perm = np.zeros(p, dtype=np.int64)
    perm[reps] = img
    perm[p - reps] = (p - img) % p
    return perm


def spread_test(values: np.ndarray, p: int, rng: np.random.Generator, B: int,
                symmetric: bool = False) -> dict:
    """Is the occupied multiset spread uniformly over [0,p)? Null: random relabelling of F_p.

    If the occupancy is invariant under v -> -v (true for y always, and for x when
    x -> -x is an automorphism), the binned counts are mirror-symmetric and the naive
    relabelling null has the wrong variance; `symmetric=True` uses relabellings that commute
    with negation. Both p-values are reported so the artifact stays visible.
    """
    nb = 64
    occ = np.bincount(values, minlength=p)
    bins = np.arange(p) * nb // p
    obs = chi2_counts(np.bincount(bins, weights=occ, minlength=nb))

    def null(sym):
        out = np.empty(B)
        for i in range(B):
            perm = _symmetric_relabelling(p, rng) if sym else rng.permutation(p)
            # value v moves to perm[v]; bin the moved occupancy
            out[i] = chi2_counts(np.bincount(bins[perm], weights=occ, minlength=nb))
        return out

    res = {"chi2": obs, "occupancy_negation_invariant": bool(np.array_equal(occ[1:], occ[1:][::-1]))}
    for label, sym in [("naive", False), ("symmetric", True)] if symmetric else [("naive", False)]:
        nl = null(sym)
        res[label] = {"null_mean": float(nl.mean()), "null_var": float(nl.var()),
                      "p_upper": mc_pvalue(obs, nl), "p_lower": float((1 + np.sum(nl <= obs)) / (1 + B))}
    return res


# ------------------------------------------------------------------ per-curve analysis


def theory_y_profile(rec) -> str:
    if rec.family == "j0":
        return ("each y occurs exactly once (x -> x^3 bijective)" if rec.p % 3 == 2 else
                "y occurs 0 or 3 times (x^3 = y^2 - b has 0 or 3 roots), once if y^2 = b")
    return "y occurs 0/1/2/3 times = #roots of a cubic; Chebotarev (S3): 1/3, 1/2, 0, 1/6"


def automorphism_check(rec, x, y) -> dict | None:
    """For j=0 / j=1728 find lambda with lambda*G = aut(G) and test it on every k."""
    p, n = rec.p, rec.order
    if rec.family == "j0" and p % 3 == 1:
        s = sqrt_mod(-3 % p, p)
        zetas = sorted({(-1 + s) * pow(2, -1, p) % p, (-1 - s) * pow(2, -1, p) % p})
        maps = [("x->beta x, y->y", lambda xx, yy, z=z: (z * xx % p, yy)) for z in zetas]
    elif rec.family == "j1728" and p % 4 == 1:
        i = sqrt_mod(p - 1, p)
        maps = [("x->-x, y->i y", lambda xx, yy, i=i: ((-xx) % p, i * yy % p))]
    else:
        return None
    k = np.arange(1, n, dtype=np.int64)
    out = []
    for label, phi in maps:
        target = phi(int(x[0]), int(y[0]))
        hits = np.flatnonzero((x == target[0]) & (y == target[1]))
        lam = int(hits[0]) + 1
        idx = (lam * k) % n
        ok = idx != 0
        px, py = phi(x[ok], y[ok])
        frac = float(np.mean((x[idx[ok] - 1] == px) & (y[idx[ok] - 1] == py)))
        out.append({"map": label, "lambda": lam, "lambda_minimal_poly_check":
                    {"l^2+l+1 mod n": (lam * lam + lam + 1) % n, "l^2+1 mod n": (lam * lam + 1) % n},
                    "fraction_of_k_satisfied": frac})
    return out


def torsion_check(rec, x, y) -> dict:
    p, a, b, n = rec.p, rec.a, rec.b, rec.order
    res = {}
    if n % 2 == 0:
        res["k=n/2: y==0"] = bool(y[n // 2 - 1] == 0)
        res["#k with y_k=0"] = int(np.sum(y == 0))
    if n % 3 == 0:
        def psi3(t):
            return (3 * t**4 + 6 * a * t * t + 12 * b * t - a * a) % p
        res["k=n/3,2n/3: psi_3(x)==0"] = all(psi3(int(x[k - 1])) == 0 for k in (n // 3, 2 * n // 3))
        res["#k with psi_3(x_k)=0"] = int(sum(psi3(int(t)) == 0 for t in set(x.tolist())))
    return res


def analyse(rec, B: int, rng: np.random.Generator) -> dict:
    p, n = rec.p, rec.order
    x, y = sequence(rec.name)
    m = n - 1
    k = np.arange(1, n)
    out = {"curve": rec.name, "family": rec.family, "p": p, "a": rec.a, "b": rec.b, "n": n,
           "generator": list(rec.generator), "notes": rec.properties["notes"]}

    # A1 value multiplicities
    out["x_multiplicity_profile"] = multiplicity_profile(x, p)
    out["y_multiplicity_profile"] = multiplicity_profile(y, p)
    out["y_profile_theory"] = theory_y_profile(rec)
    if rec.family == "generic":
        prof = out["y_multiplicity_profile"]
        obs = np.array([prof.get(0, 0), prof.get(1, 0), prof.get(3, 0)])
        exp = p * np.array([1 / 3, 1 / 2, 1 / 6])
        exp *= obs.sum() / exp.sum()
        out["y_profile_vs_chebotarev"] = {
            "observed_0_1_3": obs.tolist(), "expected": exp.round(1).tolist(),
            "chi2_p_value": float(stats.chisquare(obs, exp).pvalue),
            "count_mult_2": prof.get(2, 0)}
    if rec.family == "j0":
        out["x_cubed_multiplicity_profile"] = multiplicity_profile(x**3 % p, p)

    # A2 spread of occupied values
    occx = np.bincount(x, minlength=p)
    out["spread_x"] = spread_test(x, p, rng, B, symmetric=bool(np.array_equal(occx[1:], occx[1:][::-1])))
    out["spread_y"] = spread_test(y, p, rng, B, symmetric=True)

    # A3 negation symmetry
    out["negation"] = {
        "frac_x_{n-k}=x_k": float(np.mean(x[::-1] == x)),
        "frac_y_{n-k}=-y_k": float(np.mean((x[::-1] + 0 == x) & ((y[::-1] + y) % p == 0))),
        "full_period_corr_x_k_x_{n-k}": float(np.corrcoef(x, x[::-1])[0, 1]),
        "full_period_lag1_acf": float(acf(x / p, 1)[0]),
    }

    # A4 ordering statistics on the half period
    H = m // 2
    xh = x[:H]
    obs = ordering_stats(xh, p)
    null = {s: [] for s in obs}
    for _ in range(B):
        for s, v in ordering_stats(rng.permutation(xh), p).items():
            null[s].append(v)
    out["ordering"] = {s: {"observed": obs[s], "null_median": float(np.median(null[s])),
                           "null_q99": float(np.quantile(null[s], 0.99)),
                           "p_mc": mc_pvalue(obs[s], np.array(null[s]))} for s in obs}
    out["acf_half_period"] = acf(xh / p, LAGS).tolist()
    out["acf_null_q005_q995"] = [None, None]  # filled by caller for the plot

    # A5 group structure and automorphisms
    out["torsion"] = torsion_check(rec, x, y)
    out["automorphisms"] = automorphism_check(rec, x, y)
    return out


# ------------------------------------------------------------------ figures


def fig_index_vs_x(plt):
    """Small (9-bit) toy curves so individual points are visible."""
    toys = [generate_curve(3, 9, "generic", require_prime_order=True, name="toy generic"),
            generate_curve(3, 9, "j0", p_mod=(1, 3), require_prime_order=True, name="toy j=0 ordinary"),
            generate_curve(4, 9, "j0", p_mod=(2, 3), name="toy j=0 supersingular")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=False)
    for ax, t, c in zip(axes, toys, SERIES):
        E = t.curve()
        G = t.generator
        if G is None:  # supersingular toy may be non-cyclic; use a point of maximal order
            import random
            rr = random.Random(0)
            G = max((E.random_point(rr) for _ in range(50)), key=lambda P: E.point_order(P))
        n = E.point_order(G)
        pts = E.multiples(G, n - 1)
        ax.scatter(range(1, n), [P[0] for P in pts], s=6, color=c, linewidths=0)
        ax.axvline(n / 2, color=BASELINE, lw=1, ls="--")
        ax.set_title(f"{t.name}\np={t.p}, a={t.a}, b={t.b}, ord(G)={n}", fontsize=9)
        ax.set_xlabel("k")
    axes[0].set_ylabel("x(kG)")
    fig.suptitle("x(kG) against k: mirror symmetry about k = n/2 is the only visible structure",
                 fontsize=10, y=1.06)
    fig.savefig(OUT / "A1_index_vs_x_toy.png")
    plt.close(fig)


def fig_pairs(plt, rng):
    rec = curve_set()[0]
    x, _ = sequence(rec.name)
    H = (rec.order - 1) // 2
    xh = x[:H]
    sh = rng.permutation(xh)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharex=True, sharey=True)
    m = 15000
    axes[0].scatter(xh[:m], xh[1:m + 1], s=0.6, color=SERIES[0], linewidths=0)
    axes[0].set_title(f"{rec.name}: (x_k, x_(k+1)), first {m} k")
    axes[1].scatter(sh[:m], sh[1:m + 1], s=0.6, color=BASELINE, linewidths=0)
    axes[1].set_title("random permutation of the same values")
    for ax in axes:
        ax.set_xlabel("x_k")
        ax.set_aspect("equal")
    axes[0].set_ylabel("x_(k+1)")
    fig.suptitle("Visually indistinguishable, yet every left-hand point lies on one algebraic "
                 "curve of bidegree (2,2) (Experiment B, H2)", fontsize=9)
    fig.savefig(OUT / "A2_consecutive_pairs.png")
    plt.close(fig)


def fig_acf(plt, results, envelopes):
    fig, ax = plt.subplots(figsize=(8, 3.6))
    lags = np.arange(1, LAGS + 1)
    lo, hi = envelopes
    ax.fill_between(lags, lo, hi, color=BASELINE, alpha=0.25, lw=0,
                    label="99% envelope, random permutations (GEN-A)")
    for r, c in zip([results[0], results[3], results[5]], SERIES):
        ax.plot(lags, r["acf_half_period"], color=c, lw=1.2, label=r["curve"])
    ax.axhline(0, color=TEXT2, lw=0.6)
    ax.set_xlabel("lag")
    ax.set_ylabel("autocorrelation of x_k / p")
    ax.set_title("Autocorrelation over the half period k = 1..(n-1)/2")
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(OUT / "A3_autocorrelation.png")
    plt.close(fig)


def fig_multiplicity(plt, results):
    sel = [results[0], results[3], results[5]]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    for ax, key, title in [(axes[0], "x_multiplicity_profile", "x-values"),
                           (axes[1], "y_multiplicity_profile", "y-values")]:
        ms = sorted({int(m) for r in sel for m in r[key]})
        w = 0.8 / len(sel)
        for i, (r, c) in enumerate(zip(sel, SERIES)):
            prof = {int(k): v for k, v in r[key].items()}
            frac = [prof.get(m, 0) / r["p"] for m in ms]
            ax.bar(np.arange(len(ms)) + (i - 1) * w, frac, width=w * 0.92, color=c, label=r["curve"])
        ax.set_xticks(range(len(ms)), [str(m) for m in ms])
        ax.set_xlabel("number of times a value of F_p occurs in the full period")
        ax.set_ylabel("fraction of F_p")
        ax.set_title(f"Occurrence profile of {title}")
    axes[1].legend(fontsize=8)
    fig.savefig(OUT / "A4_value_multiplicities.png")
    plt.close(fig)


def fig_pvalue_qq(plt, results):
    ps = np.sort([v["p_mc"] for r in results for v in r["ordering"].values()])
    q = (np.arange(1, len(ps) + 1) - 0.5) / len(ps)
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    ax.plot([0, 1], [0, 1], color=BASELINE, lw=1, label="Uniform(0,1): no signal")
    ax.scatter(q, ps, s=14, color=SERIES[0], zorder=3, label=f"{len(ps)} ordering tests")
    ax.set_xlabel("expected p-value quantile")
    ax.set_ylabel("observed Monte-Carlo p-value")
    ax.set_title("Ordering tests vs permutation null")
    ax.legend(fontsize=8, loc="upper left")
    ax.set_aspect("equal")
    fig.savefig(OUT / "A5_ordering_pvalue_qq.png")
    plt.close(fig)


# ------------------------------------------------------------------ main


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=200)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    plt = setup_matplotlib()
    timing = {}
    write_curve_table()
    results = []
    for rec in curve_set():
        with timed(f"A {rec.name}", timing):
            results.append(analyse(rec, args.perms, rng))
    # acf null envelope for the plot (GEN-A)
    rec = curve_set()[0]
    xh = sequence(rec.name)[0][: (rec.order - 1) // 2]
    null_acf = np.array([acf(rng.permutation(xh) / rec.p, LAGS) for _ in range(args.perms)])
    env = (np.quantile(null_acf, 0.005, axis=0), np.quantile(null_acf, 0.995, axis=0))
    results[0]["acf_null_q005_q995"] = [env[0].tolist(), env[1].tolist()]
    with timed("A figures", timing):
        fig_index_vs_x(plt)
        fig_pairs(plt, rng)
        fig_acf(plt, results, env)
        fig_multiplicity(plt, results)
        fig_pvalue_qq(plt, results)
    allp = [v["p_mc"] for r in results for v in r["ordering"].values()]
    summary = {
        "provenance": provenance(SEED), "permutations": args.perms, "timing_s": timing,
        "n_ordering_tests": len(allp), "n_p_below_0.05": int(sum(q < 0.05 for q in allp)),
        "expected_below_0.05_under_null": 0.05 * len(allp),
        "min_ordering_p": min(allp), "bonferroni_threshold": 0.05 / len(allp),
        "curves": results,
    }
    save_json(summary, OUT / "summary.json")
    write_markdown(summary)
    print(f"ordering tests: {len(allp)}, p<0.05: {summary['n_p_below_0.05']}, min p={min(allp):.4f}")


def fmt(sp, key):
    if key not in sp:
        return "n/a"
    return f"{sp[key]['p_upper']:.3f}/{sp[key]['p_lower']:.3f}"


def write_markdown(s):
    L = ["# Experiment A summary (auto-generated)", "",
         f"Permutations per null: {s['permutations']}; seed {s['provenance']['seed']}.", "",
         "## Ordering statistics (half period, Monte-Carlo p-values vs random permutation)", "",
         "| curve | " + " | ".join(s["curves"][0]["ordering"]) + " |",
         "|---|" + "---|" * len(s["curves"][0]["ordering"])]
    for r in s["curves"]:
        L.append(f"| {r['curve']} | " + " | ".join(f"{v['p_mc']:.3f}" for v in r["ordering"].values()) + " |")
    L += ["", f"{s['n_p_below_0.05']} of {s['n_ordering_tests']} p-values < 0.05 "
          f"(expected under the null: {s['expected_below_0.05_under_null']:.1f}); "
          f"min p = {s['min_ordering_p']:.4f}; Bonferroni threshold {s['bonferroni_threshold']:.5f}.", "",
          "## Value distribution and exact symmetries", "",
          "Spread p-values are (upper/lower) tails; 'naive' = arbitrary relabelling of F_p, "
          "'sym' = relabelling commuting with v -> -v (correct null when the occupancy is negation-invariant).", "",
          "| curve | x occurrence profile | y occurrence profile | spread-x naive | spread-x sym | spread-y naive | spread-y sym | x_{n-k}=x_k | y_{n-k}=-y_k |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in s["curves"]:
        L.append(f"| {r['curve']} | {r['x_multiplicity_profile']} | {r['y_multiplicity_profile']} | "
                 f"{fmt(r['spread_x'], 'naive')} | {fmt(r['spread_x'], 'symmetric')} | "
                 f"{fmt(r['spread_y'], 'naive')} | {fmt(r['spread_y'], 'symmetric')} | "
                 f"{r['negation']['frac_x_{n-k}=x_k']:.3f} | {r['negation']['frac_y_{n-k}=-y_k']:.3f} |")
    L += ["", "## Group-structure and automorphism checks", ""]
    for r in s["curves"]:
        L.append(f"- **{r['curve']}** torsion: {r['torsion']}; automorphisms: {r['automorphisms']}")
        if "y_profile_vs_chebotarev" in r:
            L.append(f"  - y-profile vs Chebotarev: {r['y_profile_vs_chebotarev']}")
        if "x_cubed_multiplicity_profile" in r:
            L.append(f"  - x^3 occurrence profile: {r['x_cubed_multiplicity_profile']}")
    (OUT / "summary.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
