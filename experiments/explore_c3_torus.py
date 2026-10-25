# EXPLORATORY (not a recorded evidence run). Usage: python3 experiments/explore_c3_torus.py
# Degenerate (nodal / multiplicative) analogue of C3: x_k = zeta^k + zeta^-k in F_l, l = 1 mod q.
import sys; sys.path.insert(0, '/home/user/experiments'); sys.path.insert(0, '/home/user')
import numpy as np
from exp_b_hypotheses import monomials, nullspace_mod_p
from exp_c3_small_subgroups import conjectured_dmin
from src.ec import is_probable_prime
for q in (19, 31, 37, 43, 61, 73, 97, 59, 71, 101):
    l = next(l for l in range(8193, 10**6) if l % q == 1 and is_probable_prime(l))
    g = next(g for g in range(2, l) if pow(g, (l - 1) // q, l) != 1)
    z = pow(g, (l - 1) // q, l)
    x = np.array([(pow(z, k, l) + pow(z, -k % q, l)) % l for k in range(1, q)], dtype=np.int64)
    M = (q - 1) // 2; d_int = M // 2; half = np.arange(1, M + 1)
    def has(m, d): return bool(nullspace_mod_p(monomials(x[half - 1], x[(m * half) % q - 1], d, 1, l), l))
    agree, diffs = 0, []
    for m in range(1, M + 1):
        lo, hi = 1, d_int
        while lo < hi:
            mid = (lo + hi) // 2
            lo, hi = (lo, mid) if has(m, mid) else (mid + 1, hi)
        pred = conjectured_dmin(m, q, d_int)[0]
        agree += lo == pred
        if lo != pred: diffs.append((m, lo, pred))
    print(f"torus q={q} (l={l}): C3 matches {agree}/{M}; (m, observed, C3) differences: {diffs[:8]}", flush=True)
