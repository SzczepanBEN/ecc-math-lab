"""Curve set, sequence computation and plotting style shared by Experiments A and B.

The curve set is fully determined by MASTER_SEED and the search rules below; the chosen
parameters are written to results/curves.json for the record.
"""
from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.curvegen import CurveRecord, generate_curve  # noqa: E402
from src.harness import RESULTS, save_json  # noqa: E402

MASTER_SEED = 20261025
BITS = 17  # p in [2^16, 2^17): every full period is ~1e5 points, ~0.3 s to enumerate


def _search(family, *, start, p_mod=None, prime=False, want=lambda r: True, name):
    """First seed >= start whose curve satisfies `want` (deterministic)."""
    for seed in range(start, start + 5000):
        r = generate_curve(seed, BITS, family, require_prime_order=prime, p_mod=p_mod, name=name)
        if r.structure[1] == 1 and want(r):
            return r
    raise RuntimeError(name)


@lru_cache(maxsize=1)
def curve_set() -> tuple[CurveRecord, ...]:
    s = MASTER_SEED
    recs = (
        # two generic prime-order curves (j != 0, 1728)
        _search("generic", start=s, prime=True, name="GEN-A"),
        _search("generic", start=s + 1000, prime=True, name="GEN-B"),
        # generic cyclic curve with composite order divisible by 6 (small subgroups present)
        _search("generic", start=s + 2000, name="GEN-C6",
                want=lambda r: r.order % 6 == 0 and r.properties["largest_prime_factor"] > 1000),
        # j = 0, p = 1 mod 3: ordinary, extra automorphism of order 3 (secp256k1's family)
        _search("j0", start=s, p_mod=(1, 3), prime=True, name="J0-A"),
        _search("j0", start=s + 1000, p_mod=(1, 3), prime=True, name="J0-B"),
        # j = 0, p = 2 mod 3: supersingular, #E = p + 1
        _search("j0", start=s, p_mod=(2, 3), name="J0-SS"),
        # j = 1728, p = 1 mod 4: ordinary, automorphism of order 4 (x,y) -> (-x, iy)
        _search("j1728", start=s, p_mod=(1, 4), name="J1728"),
    )
    return recs


def write_curve_table() -> Path:
    return save_json({"master_seed": MASTER_SEED, "bits": BITS,
                      "curves": [r.to_json() for r in curve_set()]}, RESULTS / "curves.json")


@lru_cache(maxsize=None)
def sequence(name: str) -> tuple[np.ndarray, np.ndarray]:
    """x_k, y_k for k = 1..n-1 where n = ord(G) = #E (the identity at k = n is dropped).

    Index i of the returned arrays corresponds to k = i + 1.
    """
    r = next(c for c in curve_set() if c.name == name)
    E = r.curve()
    pts = E.multiples(r.generator, r.order - 1)
    assert E.add(pts[-1], r.generator) is None, "generator order check failed"
    x = np.fromiter((P[0] for P in pts), dtype=np.int64, count=len(pts))
    y = np.fromiter((P[1] for P in pts), dtype=np.int64, count=len(pts))
    return x, y


# ---- plotting style (reference palette from the dataviz skill; light mode, static PNG) ----
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
BASELINE = "#8a8985"
TEXT, TEXT2, SURFACE, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"


def setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID, "axes.labelcolor": TEXT2, "xtick.color": TEXT2,
        "ytick.color": TEXT2, "text.color": TEXT, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
        "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "lines.linewidth": 2,
        "legend.frameon": False, "axes.axisbelow": True, "savefig.dpi": 130, "savefig.bbox": "tight",
    })
    return plt
