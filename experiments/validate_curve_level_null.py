"""Independent validation of an apparent anomaly found in Experiment B / H5 (first run).

With 300 random curves over GEN-B's prime (Random(3)), the 64-bin chi^2 of x-occupancy looked
smaller than under a random-relabelling null (KS p = 0.0006). This script checks whether that is a
property of elliptic curves or a fluctuation of the null sample:
  1. 120 independent seeds x 300 curves: distribution of the 300-curve mean;
  2. pooled curves vs the chi^2_{63} reference;
  3. other primes and bin counts, comparing curves with relabelled versions of themselves.
Runtime ~1-2 minutes.  Output: results/exp_b/validate_curve_level_null.json
"""
from __future__ import annotations

import random

import numpy as np

from common import RESULTS, curve_set, save_json
from src.ec import Curve


def chi2_spread_factory(p: int, nb: int):
    T = Curve(p, 1, 1)._square_table if p != 0 else None
    xx = np.arange(p, dtype=np.int64)
    x3 = xx * xx % p * xx % p
    bins = xx * nb // p

    def occ(a, b):
        return T[(x3 + a * xx + b) % p]

    def chi(o, relabel=None):
        bb = bins if relabel is None else bins[relabel]
        c = np.bincount(bb, weights=o, minlength=nb)
        e = c.sum() / nb
        return float(((c - e) ** 2 / e).sum())
    return occ, chi


def main():
    out = {}
    p = next(c for c in curve_set() if c.name == "GEN-B").p
    occ, chi = chi2_spread_factory(p, 64)
    means, pooled = {}, []
    for seed in [3] + list(range(100, 220)):
        cr = random.Random(seed)
        vals = []
        while len(vals) < 300:
            a, b = cr.randrange(1, p), cr.randrange(1, p)
            if (4 * a**3 + 27 * b * b) % p:
                vals.append(chi(occ(a, b)))
        means[seed] = float(np.mean(vals))
        if seed != 3:
            pooled += vals
    m = np.array([v for s, v in means.items() if s != 3])
    pooled = np.array(pooled)
    out["seed_sweep"] = {
        "p": p, "seed3_mean": means[3], "other_seeds_mean_of_means": float(m.mean()),
        "other_seeds_sd_of_means": float(m.std()), "predicted_sd": float(pooled.std() / np.sqrt(300)),
        "min_other_seed_mean": float(m.min()),
        "n_other_seeds_at_or_below_seed3": int((m <= means[3]).sum()), "n_other_seeds": len(m),
        "pooled_n": len(pooled), "pooled_mean": float(pooled.mean()), "pooled_var": float(pooled.var()),
        "chi2_63_reference": {"mean": 63, "var": 126},
    }
    rows = []
    nrng = np.random.default_rng(7)
    for q in [10007, 99377, 1000003]:
        for nb in [16, 64]:
            occ_q, chi_q = chi2_spread_factory(q, nb)
            cr = random.Random(7)
            cur, rel = [], []
            ncur = 2000 if q < 200000 else 400
            while len(cur) < ncur:
                a, b = cr.randrange(q), cr.randrange(q)
                if (4 * a**3 + 27 * b * b) % q == 0:
                    continue
                o = occ_q(a, b)
                cur.append(chi_q(o))
                if len(rel) < ncur // 3:
                    rel.append(chi_q(o, nrng.permutation(q)))
            cur, rel = np.array(cur), np.array(rel)
            rows.append({"p": q, "bins": nb, "curves_mean": float(cur.mean()),
                         "curves_se": float(cur.std() / np.sqrt(len(cur))),
                         "relabelled_mean": float(rel.mean()), "relabelled_se": float(rel.std() / np.sqrt(len(rel))),
                         "z": float((cur.mean() - rel.mean()) /
                                    np.hypot(cur.std() / np.sqrt(len(cur)), rel.std() / np.sqrt(len(rel))))})
            print(rows[-1])
    out["other_primes"] = rows
    save_json(out, RESULTS / "exp_b" / "validate_curve_level_null.json")
    print(out["seed_sweep"])


if __name__ == "__main__":
    main()
