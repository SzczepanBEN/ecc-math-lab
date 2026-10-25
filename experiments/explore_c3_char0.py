# EXPLORATORY (not a recorded evidence run). Usage: python3 experiments/explore_c3_char0.py [tau]
# Exploratory: does C3 hold in characteristic 0?  x_k = wp(k/q) on C/(Z + tau Z), generic tau.
import sys, cypari2, mpmath as mp
sys.path.insert(0, '/home/user/experiments'); sys.path.insert(0, '/home/user')
from exp_c3_small_subgroups import conjectured_dmin
pari = cypari2.Pari(); pari.set_real_precision(130); mp.mp.dps = 120
def wp_values(q, tau):
    L = pari('[1, %s]' % tau)
    return [mp.mpc(str(pari.real(v)).replace(' E', 'e'), str(pari.imag(v)).replace(' E', 'e'))
            for v in (pari.ellwp(L, pari(f'{k}/{q}')) for k in range(1, q))]
def sv_ratio(X, Y, d):
    rows = [[x**i for i in range(d + 1)] + [x**i * y for i in range(d + 1)] for x, y in zip(X, Y)]
    A = mp.matrix(rows)
    for j in range(A.cols):                       # column scaling
        s = mp.sqrt(sum(abs(A[i, j])**2 for i in range(A.rows)))
        for i in range(A.rows): A[i, j] /= s
    S = mp.svd_c(A, compute_uv=False)
    return mp.nstr(min(S) / max(S), 3)
def dmin_numeric(x, q, m, d_int, tol=mp.mpf(10)**-60):
    M = (q - 1) // 2
    X = [x[k - 1] for k in range(1, M + 1)]; Y = [x[(m * k) % q - 1] for k in range(1, M + 1)]
    for d in range(1, d_int + 1):
        if 2 * d + 2 > M:          # wide matrix: a null vector exists by linear algebra
            return d
        rows = mp.matrix([[xx**i for i in range(d + 1)] + [xx**i * yy for i in range(d + 1)] for xx, yy in zip(X, Y)])
        for j in range(rows.cols):
            s = mp.sqrt(sum(abs(rows[i, j])**2 for i in range(rows.rows)))
            for i in range(rows.rows): rows[i, j] /= s
        S = mp.svd_c(rows, compute_uv=False)
        if min(S) / max(S) < tol: return d
    return None
import sys as _s
tau = _s.argv[1] if len(_s.argv) > 1 else '0.3141592 + 1.2718281*I'
for q in (19, 31, 37, 43):
    x = wp_values(q, tau); M = (q - 1) // 2; d_int = M // 2
    agree = tot = 0; bad = []
    for m in range(1, M + 1):
        pred = conjectured_dmin(m, q, d_int)[0]
        obs = dmin_numeric(x, q, m, d_int)
        tot += 1; agree += (obs == pred)
        if obs != pred: bad.append((m, obs, pred))
        below = sum(1 for mm in range(1, M + 1) if conjectured_dmin(mm, q, d_int)[0] < d_int) if m == M else None
    print(f"char 0, q={q}, tau={tau}: C3 matches {agree}/{tot} ({below} predicted below threshold); mismatches {bad[:6]}", flush=True)
# sharpness for one case: singular-value ratio just below / at the predicted degree
q = 37; x = wp_values(q, tau); m = 12            # m = -1/3 mod 37 (3*12 = 36)
X = [x[k-1] for k in range(1, 19)]; Y = [x[(m*k) % q - 1] for k in range(1, 19)]
print("q=37, m=-1/3: sigma_min/sigma_max at d=5:", sv_ratio(X, Y, 5), " at d=6 (C3):", sv_ratio(X, Y, 6))
