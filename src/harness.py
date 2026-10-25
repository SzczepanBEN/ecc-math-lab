"""Minimal experiment harness: seeding, timing, provenance and JSON output."""
from __future__ import annotations

import json
import platform
import random
import sys
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def provenance(seed: int) -> dict:
    import matplotlib
    import scipy
    return {
        "seed": seed,
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "platform": platform.platform(),
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def rngs(seed: int) -> tuple[random.Random, np.random.Generator]:
    return random.Random(seed), np.random.default_rng(seed)


def _default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, frozenset)):
        return sorted(o)
    raise TypeError(type(o))


def save_json(obj, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=_default))
    return path


@contextmanager
def timed(label: str, log: dict):
    t0 = time.perf_counter()
    yield
    log[label] = round(time.perf_counter() - t0, 3)
    print(f"[{label}] {log[label]:.2f}s", flush=True)
