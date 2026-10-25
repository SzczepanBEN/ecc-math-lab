# Research log — Phase 1

Chronological record of what was done, what was expected, what happened, and what went wrong.
All numbers below are copied from `results/` files produced by `./run_all.sh`.

---

## 0. Inspection of the existing project

**Expected:** `src/ec.py`, `src/curvegen.py`, `src/harness.py` and an incomplete `src/dlog.py`
from a previous session.

**Observed:** `/home/user` was empty. A search of the whole filesystem (`find / -xdev -name ec.py
-o -name curvegen.py -o -name harness.py -o -name dlog.py`, excluding system packages) found no
project files. There was nothing to reuse or audit, and nothing that could break imports.

**Action:** I wrote `src/ec.py`, `src/curvegen.py` and `src/harness.py` from scratch and
verified them independently (section 1). I did **not** create `src/dlog.py` or any other
discrete-logarithm code, since that is outside Phase-1 scope. numpy was installed; scipy,
matplotlib and pytest were missing and were installed with pip.

## 1. Building and verifying the laboratory

Design choices:

- Affine short-Weierstrass arithmetic with exact Python integers; the identity O is `None`.
- Two independent scalar multiplications (double-and-add and Montgomery ladder) that
  cross-check each other.
- Two independent point counts: an O(p) numpy count of #{y : y² = f(x)}, and the Legendre-sum
  formula. For tiny p, a third count by brute-force double loop.
- The group structure Z/n1 × Z/n2 is computed from lcm of point orders. The tests check it
  against brute-force orders of every point and against the structural fact that exactly n2²
  points have order dividing n2.
- Curves used in experiments are cyclic, with a *certified* generator: ord(G) = #E, checked via
  the prime factorisation of #E.

Tests (`tests/`, 54 tests, all pass in about 1 s):

- Primality vs trial division up to 5000.
- Tonelli–Shanks, including p = 2^61 − 1.
- Rejection of singular curves.
- Point counts three ways on 8 small curves, including j=0 and j=1728 in both residue classes.
- Exhaustive closure, commutativity, identity and inverses on those curves, and associativity on
  3000 random triples each plus 2000 triples over p = 10^6 + 3.
- Scalar multiplication vs repeated addition up to 2N + 3; homomorphism (a+b)P = aP + bP and
  a(bP) = (ab)P; Lagrange.
- Group structure vs brute force.
- Seeded reproducibility of curve generation.
- Classification of supersingular, anomalous, j=0 and j=1728 curves, and CM field Q(√−3) / Q(i).
- External test vectors:
  - y² = x³ + 2x + 2 over F_17, with #E = 19 and 2(5,1) = (6,3) (Paar–Pelzl textbook example).
  - y² = x³ + x + 1 over F_23, with #E = 28.
  - The published secp256k1 parameters: G on the curve, n·G = O, and exactly two (β, λ) pairs
    from the cube roots of unity satisfying λG = (βx_G, y_G). This is the GLV endomorphism.

## 2. Curve set

`experiments/common.py`, `MASTER_SEED = 20261025`, 17-bit p. Each curve is the first seed
(from a fixed starting offset) satisfying the family's rule. Parameters are in `README.md` and
`results/curves.json`. There are 7 curves:

- Two generic prime-order curves (GEN-A, GEN-B).
- One generic cyclic curve with n = 2·3·20641 (GEN-C6), so that small subgroups are present.
- Two ordinary j = 0 prime-order curves (J0-A, J0-B), the secp256k1 family.
- One supersingular j = 0 curve (J0-SS, p ≡ 2 mod 3).
- One ordinary j = 1728 curve (J1728).

## 3. Literature check (before designing tests)

- The sequence P_k = kG = P_{k−1} + G is exactly the **elliptic-curve linear congruential
  generator (EC-LCG)**. Hess & Shparlinski (Des. Codes Cryptogr. 35, 2005) prove high linear
  complexity and good multidimensional distribution. Related work includes Kohel & Shparlinski
  (ANTS-IV) on exponential sums over curve points and El Mahassni & Shparlinski on the elliptic
  curve power generator.
- Algebraic relations between x-coordinates of points summing to O are **Semaev's summation
  polynomials** (ePrint 2004/031). S3 has degree 2 in each variable.
- The j=0 endomorphism (x, y) ↦ (βx, y) is the basis of the **GLV method** (Gallant–Lambert–
  Vanstone, CRYPTO 2001).
- Polynomial interpolation of the EC discrete logarithm has proven degree lower bounds
  (**Lange & Winterhof**, COCOON 2002; with improved bounds for curves with endomorphisms of
  degree 2 or 3 in follow-up work).

Consequence: nothing I could find "by pattern" about this sequence is likely to be new.
The value of this phase is a careful, reproducible map of which structure exists and why.

## 4. Experiment A — geometry of the sequence

Design (`experiments/exp_a_geometry.py`). The key decision was to separate two different
questions, each with its own null model:

1. **Which values occur** (value multiset over the full period). Over a full period the sequence
   visits every point exactly once, so its coordinate distribution is the distribution of the
   point set and says nothing about the scalar multiplication order. Null model: the same
   multiset under a random relabelling of F_p.
2. **In what order they occur** (k ↦ x_k). Because x_{n−k} = x_k, I restrict to the half period
   k = 1..(n−1)/2, where every x-value appears once. Null model: random permutations of the same
   half-period values (200 per curve), which holds the multiset fixed exactly.

Ordering statistics:

- max |ACF| over lags 1–100
- Ljung–Box Q(100)
- χ² of consecutive differences (x_{k+1} − x_k) mod p, 64 bins
- 2-D serial χ² of (x_k, x_{k+1}), 16 × 16
- |rank correlation(k, x_k)|
- KS distance of the first 1000 terms vs uniform
- periodogram peak/mean

### Results

- **Exact symmetries.** x_{n−k} = x_k and y_{n−k} = −y_k for 100% of k on all 7 curves.
  Over the full period, corr(x_k, x_{n−k}) = 1.000 while the lag-1 autocorrelation is about 0
  (|r| ≤ 0.0055). Autocorrelation is blind to a perfect reflection symmetry.
- **Value multiplicities.** These match theory exactly:
  - Generic curves: x occurs 0 or 2 times (once for 2-torsion). y occurs 0/1/2/3 times, equal to
    the number of roots of the cubic f(x) − y². The frequencies match the Chebotarev (S3)
    prediction 1/3 : 1/2 : 1/6 (χ² p = 0.39, 0.88, 0.96). The rare count-2 y-values (2 on GEN-B)
    are the critical values of f, where the cubic has a double root.
  - j=0 ordinary: y occurs only 0 or 3 times. x³ occurs 0 or 6 times, which is the orbit of the
    order-6 automorphism group.
  - J0-SS: **every y ∈ F_p occurs exactly once**, because x ↦ x³ is a bijection when
    p ≡ 2 (mod 3). This is perfectly uniform and hence "too uniform" for a random-set test.
  - J1728: the x-occupancy is invariant under x ↦ −x.
- **Ordering tests.** 1 of 49 Monte-Carlo p-values is below 0.05 (expected under the null:
  2.5). The minimum is p = 0.015 (J0-B, difference χ²), against a Bonferroni threshold of
  0.001. The p-value QQ plot (A5) lies near the diagonal, slightly on the conservative side.
  **No ordering structure is detected on any curve family.**
- **Group structure.** On even-order curves the unique point with y = 0 sits at k = n/2. On
  curves with 3 | n, the points at k = n/3 and 2n/3 are the zeros of the 3-division polynomial
  ψ3 (inflection points).
- **Automorphisms.** On J0-A and J0-B, x_{λk} = βx_k for 100% of k, with λ² + λ + 1 ≡ 0
  (mod n). On J1728, x_{λk} = −x_k with λ² ≡ −1 (mod n).

### Failure 1 — mis-specified null (method artifact, caught and fixed)

The first run gave y-spread p-values of 0.025, 0.010 (lower tail) and 0.050 on three curves:
too many extreme values. **Diagnosis:** y and −y always have the same occupancy, so the 64 bin
counts are mirror-symmetric (bin j ↔ bin 63 − j). An arbitrary relabelling of F_p destroys that
symmetry, which gives the null roughly half the true variance. **Fix:** a relabelling null that
commutes with v ↦ −v (permute the pairs {v, −v} and flip signs). Null variances roughly
doubled (GEN-A: 132 → 305), and the extreme p-values disappeared (0.065, 0.060, 0.080). Both
p-values are still reported in `results/exp_a/summary.md`, so the artifact stays visible.

### Failure 2 — a rendering artifact

In the first version of figure A1, gridlines were drawn above the scatter markers. Points lying
on a gridline value looked hollow and seemed to line up at x ∈ {0, 50, 100, …}, an apparent
"pattern" produced by the plotting library. Fixed with `axes.axisbelow = True`.

## 5. Experiment B — hypotheses

### H1 — translation relations (Semaev S3)

- **Relationship:** for every lag ℓ, every pair (x_k, x_{k+ℓ}) satisfies one fixed polynomial
  F_ℓ(X, Y) = 0 of bidegree (2, 2). F_ℓ = S3(X, Y, x(ℓG)). The bidegree drops to (1, 1) when ℓG
  is a 2-torsion point T = (t, 0).
- **Motivation:** P_k + ℓG − P_{k+ℓ} = O, so Semaev's S3 vanishes. Given x(P), x(P + Q) is one of
  the two values x(P ± Q), so the degree is 2 in each variable. When 2Q = O, P + Q and P − Q
  coincide and the map becomes a Möbius map. A direct computation gives
  x(P + T) = t + (3t² + a)/(x − t).
- **Prediction:** a 1-dimensional nullspace at (2, 2), equal to S3, and none at (1,1), (1,2) or
  (2,1). At ℓ = n/2, a (1, 1) relation (Y − t)(X − t) = 3t² + a. Lag n/3 (3-torsion) is *not*
  special. Shuffled baselines have no relation at all.
- **Test:** blind linear algebra over F_p. Build monomials X^i Y^j for 400 sampled pairs,
  compute the nullspace, then verify the found polynomial on **every** pair. Compare the
  coefficients with an independent expansion of S3.
- **Result:** confirmed on all 7 curves. Five generic lags (1, 2, 17, 1000, ≈ n/3) were tested,
  plus n/2 and n/3 where applicable. Every (2, 2) relation vanishes on 100% of pairs and equals
  S3 coefficient for coefficient. At ℓ = n/2 (GEN-C6, J0-SS, J1728) the (1, 1) relation appears
  and equals the predicted Möbius form. At ℓ = n/3 the bidegree stays (2, 2). Baselines: nullspace
  dimension 0 at both (1, 1) and (2, 2).
- **Classification:** a **known mathematical property** (group law / summation polynomials).
  The structure does not depend on the curve family.

### H2 — index-scaling symmetries equal the automorphism group

- **Relationship:** the set of m ∈ [1, n−1] for which x_{mk} = M(x_k) for a single Möbius map M
  has size exactly |Aut(E)|: 2 for generic curves, 6 for j = 0 with p ≡ 1 (mod 3), 4 for j = 1728
  with p ≡ 1 (mod 4), and 2 for supersingular j = 0 with p ≡ 2 (mod 3), where β ∉ F_p.
- **Motivation:** automorphisms act on x by x ↦ ζx (ζ a root of unity), so they give such m.
  The converse ("no other m") is my **conjecture/heuristic**, not a theorem I have proved. It
  matches the idea that a degree-1 x-map on all of E(F_p) should come from a degree-1 isogeny,
  but [m] and [m + n] agree on E(F_p) while having different degrees, so the naive degree
  argument does not apply directly.
- **Test:** a blind search over *all* m. Screen with two cross-ratio invariants, which are exact
  for Möbius maps (false-positive rate about 1/p² per m). Fit survivors with a (1, 1) relation and
  verify it on every k. Run on the 7 main curves plus 24 further 14-bit curves (8 per family).
- **Result:** 31/31 curves match |Aut(E)| exactly.
  - J0-A: m ∈ {1, 4868, 4869, 114238, 114239, 119106} = {±1, ±λ, ±λ²}.
  - J1728: m ∈ {±1, ±19785}.
  - Generic curves and J0-SS: only ±1.
- **Classification:** **known mathematical property** (automorphisms, GLV). The "only
  automorphisms" part has strong experimental support but remains a labelled conjecture.

### H3 — no low-complexity description of k ↦ x_k

- **Relationship (null hypothesis):** no polynomial F(k, x_k) = 0 of small bidegree exists, and
  the linear complexity of (x_k) is about N/2, the typical value for random sequences.
- **Motivation:** Lange–Winterhof-type lower bounds on interpolation degree; Hess–Shparlinski
  linear-complexity bounds for the EC-LCG.
- **Test:** the nullspace method on k^i x_k^j with 2M + 40 rows for M monomials, at bidegrees
  (1,1) through (60,1) and (10,10). Positive control: the planted sequence k⁴ + 5k + 1.
  Baseline: shuffled x. Berlekamp–Massey linear-complexity profile on the first 800 terms of
  x_k and y_k. Control: the polynomial sequence (L = 5).
- **Result:** real and shuffled sequences give nullspace dimension 0 at every bidegree on every
  curve. The control is detected exactly when d1 ≥ 4 and is correctly *not* detected at (1,1),
  (2,2) and (3,1). Linear complexity: L_800 = 400 for x and y on GEN-A, J0-A and J0-SS, identical
  to the shuffled baseline; the control gives L = 5.
- **Caveat found during analysis:** over a large field a random sequence has the "perfect"
  profile L_N = ⌈N/2⌉ with probability close to 1. The test therefore excludes low-complexity
  structure, but it cannot rank sequences that are already random-looking. The degrees tested
  (≤ 60) are tiny compared with the proven lower bounds, so this is a sanity check, not a
  sharp test.
- **Classification:** **negative result**, consistent with known theory. There is no
  low-degree link between the index and the coordinate.

### H4 — "consecutive x's are predictable" by entropy

- **Relationship:** the plug-in estimate of H(x_{k+1} | x_k) is about 1 bit, much less than
  log₂ p ≈ 17 bits, suggesting that x_{k+1} is nearly determined by x_k.
- **Prediction:** this is a **statistical artifact**. Over the full period every x-value occurs
  exactly twice, so *any* ordering of the same multiset gives ≤ 1 bit. Over the half period every
  value occurs once, giving 0 bits for anything. The genuine dependence is algebraic (H1) and
  needs x_G.
- **Result:** the plug-in estimate is 1.0000 bits for the real sequence vs 1.0000 (0.9999 on
  J0-SS) for the permuted one. Over the half period both give 0.0000. The S3 test separates them
  completely: 100% of k for the real sequence vs 0.0014%–0.0061% for the permuted one
  (1–7 coincidences per curve). The chance rate is ≈ 4/p, not 2/p as I first wrote: the two roots
  of S3(x_k, ·, x_G) are curve x-values, and each appears twice in the shuffled multiset.
- **Classification:** the entropy signal is a **statistical artifact** of the estimator; the
  underlying algebraic dependence is the known H1 relation.

### H5 — are the nominal Experiment-A hits real?

- **(a) J0-B difference χ² (p = 0.015 in A).** With 2000 permutations, p = 0.017. If this were a
  property of the curve it should persist for other generators G' = cG, which reorder the same
  multiset. Over 30 random generators only 1/30 has p < 0.05, and the p-values look uniform.
  **Classification: statistical artifact** (multiple testing; one of 49 tests).
- **(b) GEN-B x-spread (p ≈ 0.01–0.045 in A, depending on the null draw).** I compared it with a
  *curve-level* null: random curves over the same p.
  - **First run** (300 curves, `random.Random(3)`): the curve-level null was *narrower* than the
    relabelling null (mean χ² 60.8 vs 63.2, variance 101 vs 128, KS p = 0.0006). That would mean
    real curves spread x more evenly than random sets, contradicting my calculation that
    E_{a,b}[χ(f(x))χ(f(x'))] = 0 for x ≠ x'. I recorded it as unexpected and validated it
    independently (`experiments/validate_curve_level_null.py`):
    - 120 other seeds × 300 curves: mean of means 62.90, sd 0.644 (predicted 0.644). None of the
      120 means is as low as 60.80 (minimum 60.87).
    - Pooled 36,000 curves: mean 62.90, variance 124.3, against 63 and 126 for χ²₆₃.
    - Three primes × two bin counts, curves vs relabelled versions of themselves: all |z| ≤ 1.14.
  - **Conclusion:** the seed-3 sample of the *null itself* was an approximately 3σ outlier.
    **Statistical artifact** caused by an undersized Monte-Carlo null. H5 now uses 3000 curves:
    p_curve = 0.031, p_relabel = 0.031, KS(curve vs relabel) p = 0.36. GEN-B's value is an
    ordinary tail draw, and 5.2% of random curves fall below 0.05 vs the relabelling null.
  - **Lesson:** a Monte-Carlo null is itself a random sample. Report its size, and do not compare
    two Monte-Carlo samples with a KS test without replicating.

## 6. Reproducibility

`./run_all.sh` regenerates everything in about 3 minutes on one core. Re-running Experiments A
and B produced byte-identical `summary.md` files.

## 7. Hypotheses considered but not pursued in Phase 1

- Short-interval discrepancy of x_k for k ∈ [1, K] vs proven character-sum bounds. Only K = 1000
  was tested (no signal); a systematic K-sweep is a Phase-2 candidate.
- Relations x_{mk} = R(x_k) with R of degree > 1, i.e. endomorphisms of degree 2 or 3 on curves
  with small CM discriminant. Phase-2 candidate.
- Anything requiring the index of a given point (DLP-adjacent). Out of scope.

## 8. Conclusions of Phase 1

All structure found is explained by known mathematics. Of the apparent "signals", four were
artifacts of my own methods: a symmetry-breaking null, a plotting z-order, a plug-in entropy
estimator, and an undersized Monte-Carlo null. Linear and statistical tests find no ordering
structure, while exact algebraic tests find strong structure at every lag. See `findings.md`.

---

# Phase 2 — direction 1: blind discovery of small-degree endomorphisms

Approved by the project lead after the Phase-1 review. All runs that support a claim were
executed with the computation-audit bounded runner. Each has a contract, input hashes, resource
caps and a validated manifest in `results/phase2/run-00N/`.

## P2.0 Literature check (bounded by the network policy)

Direct fetches of primary sources (arxiv.org, math.mit.edu, uwaterloo.ca) were refused by the
container's network policy. Every citation below was therefore checked **at search-result level
only**: the source was identified and its statement paraphrased by the search results, but
theorem text was not read.

| item | source identified | checked | how the project depends on it |
|---|---|---|---|
| Deuring reduction: End(Ẽ) = O for ordinary reduction at p ∤ conductor | Lang, *Elliptic Functions*, Ch. 13 §4 Thm 12, as paraphrased in arXiv math/0409501; Silverman ATAEC II.4 | statement only | interpretation of "End = O(table)"; the data independently rejects End = O_K and End = Z[π] |
| class-number-1 orders and their j-invariants | LMFDB CM-discriminant tables; class-number-one notes | values matched for all 13 | `src/cm.py`; independently checked: the CM field of t² − 4p matches on every curve |
| isogeny standard form (u/v, (s/t)·y), deg = max(deg u, deg v), kernel x-coordinates = roots of v | Sutherland, MIT 18.783 lecture notes 5–6 | statement only | the (d,1) search and the kernel-structure check |
| computing End(E) of ordinary curves | Kohel, PhD thesis (Berkeley 1996); Bisson–Sutherland, J. Number Theory 131 (2011) 815–831 | bibliographic | prior work: End(E) is computable by isogeny-volcano methods; our method is a naive small-norm detector, not a competitor |
| GLV curves with degree-2/3 endomorphisms | GLV (CRYPTO 2001); Birkner–Longa–Sica / Longa–Sica (degree-3 families) | bibliographic | context only |
| interpolation of x(P) ↦ x(mP) on subgroups | Lange–Winterhof, AAECC-15 (LNCS 2643, 2003); Kiltz–Winterhof (2004) | bibliographic, theorems **not** seen | possible prior work for P2.4 (open) |

## P2.1 Contract (claim C2)

Let E/F_p be ordinary with End(E) = O (discriminant D), let G′ have prime order q, and let
x_k = x(kG′). For d ≤ dmax define S_d = {α mod q : α ∈ O, N(α) = d}. Here α acts on ⟨G′⟩ by an
eigenvalue computed from π = (t − vD)/2 + vw ≡ 1.

- **Theorem side:** each m ∈ S_d admits a rational R of degree d with x(mP) = R(x(P)). This
  follows from the standard form of the isogeny α.
- **Conjecture C2:** no other m admits such an R of degree ≤ d (for q large relative to d).
- **Finite surrogate:**
  - an exact rank test over F_p for every m ∈ [1, q−1] and d ≤ dmax, followed by verification
    on all q−1 points;
  - comparison with S_d, which is computed only from (p, t, q, D);
  - competing hypotheses End = O(table), End = O_K and End = Z[π], whose predictions differ.
- **Exactness:** the screen has **no false negatives**. A true relation holds at the sample
  points, so the rank test must flag it. All 5424 positives were verified on all points, so the
  observed sets are exact, not estimates.

## P2.2 Main run (run-001): 13 class-number-1 orders + 4 generic controls

- p ∈ [2^13, 2^14), cofactor ≤ 8, dmax = 20 (42 for D = −163); 128 s.
- **Observed = predicted S_d for every d on all 17 curves** (exact set equality).
- Non-maximal orders D = −12, −16, −27, −28: the data rejects End = O_K and matches the table
  order, so the conductor is determined by the data.
- End = Z[π] is rejected on all 13 CM curves.
- **Blind inference:** from (m, d) alone, the unique T with m² − Tm + d ≡ 0 (mod q) gives
  T² − 4d = f²D_K. Taking the gcd of f recovers the exact order discriminant on all 13 CM
  curves; the 4 generic controls give only integers.
- **Orthogonal checks:**
  - the CM field of t² − 4p equals the table's field on every curve;
  - every found map has deg num = d, deg den = d − 1, and a denominator of the form s·h² with
    s | x³ + ax + b;
  - m = 2 reproduces the closed-form doubling formula on all 17 curves;
  - the batched screen agrees with an independent per-m rank test (0 disagreements on 2,404
    audited cases);
  - the validator accepts the manifest.
- D = −163: only integer multipliers (degrees 1, 4, 9, 16, 25, 36) up to d = 40, then exactly
  4 at d = 41 = N((1+√−163)/2).

## P2.3 Boundary test (run-002): C2 fails below the interpolation threshold

The parameter count behind C2 has 2d+2 unknowns against M = (q−1)/2 distinct pairs. It predicts
"every m" once 2d+2 > M, and "only endomorphisms" below that. On small subgroups (q ∈ [40, 160])
the first part holds, but **non-endomorphism multipliers appear well below the threshold** in
6 of 8 cases. For example, the generic curve with q = 61 has extras at d = 10 while
d_int = 15; with q = 157 they appear at 26 vs 39.

**Diagnosis:**
- An independent pure-Python rank computation agrees with the numpy result.
- The found map is a genuine degree-10 rational function in lowest terms.
- Its zeros and poles are not in S, so the coset-polynomial hypothesis is rejected.
- The extras are small-height fractions mod q: ±1/3, ±2/3 and so on.

**Bug caught along the way:** the first version of the follow-up script inverted the q² ∤ N
filter (`N % (q*q)` used where `N % (q*q) == 0` was meant). It accepted only q² | N curves and
looped forever for q = 157. It was found because the run stalled; no output from it was used.
A second slip ran a background job from the wrong directory; it failed immediately and its
non-output was not used.

## P2.4 Characterising the extras (exploratory, then run-003)

- **Rational reconstruction:** every extra multiplier is a small-height fraction a/b mod q.
- **b = 3 fit:** the first fit was d(a/3) = (M + a² − 1)/3 for q ≡ 1 (mod 6), from q = 61 and
  q = 157 (which share the earlier threshold-test curves).
- **More denominators:** from b = 5, 7, 9, 11, 19 the numerator offsets follow the pattern
  −(b+1)/4 for b ≡ 3 (mod 4) and +(b−1)/4 for b ≡ 1 (mod 4).
- **A failed first form:** taking the minimum over *all* representations matched only
  306–552/908, depending on the rounding rule. That rejected the unrestricted version.
- **The fix:** a representation counts only when the formula is an integer, i.e. when
  q ≡ 4a² (mod b). With that restriction the rule matched **908/908**. Since about 15 hand
  cases informed it, this is fitting data, not evidence.
- **Curve independence:** two curves per q over different primes p gave identical d_min
  profiles for all 10 primes q. This was not anticipated.
- **A rejected mechanism:** a coset-polynomial construction (index-3 subgroups of
  (Z/q)*/±1). The found maps have no zeros or poles in S.

## P2.5 Pre-registered held-out test (run-004)

- Formula frozen in `conjectured_dmin` before the run.
- 20 new primes q ∈ [67, 199], new seed, one fresh generic curve each.
- **Result: 1454/1454 exact matches, 253 of them below threshold.**
- CM side check (D = −3, −4, −7, small q), not part of C3: C3 matched 6/21, 11/26 and 13/39
  multipliers. The observed degree was never higher than C3, only lower or equal. This is
  consistent with C3 relations existing on every curve and CM adding more.
- **Recording slip:** the CM side-check numbers were omitted from `summary.json` because an
  edit's search string did not match and `str.replace` silently did nothing. They are in the
  hashed `stdout.txt` of run-004.
- Fixing the script changed a pinned input, so the validator marks runs 002–004 as
  "input changed". All three were therefore re-run with clean provenance as runs 005–007
  (P2.6).

## P2.5b Consistency of C3 with the main run

For b ≥ 3 the C3 degree is at least (q − 3)/6, which is ≥ 182 for every run-001 curve
(dmax ≤ 42). So C3 predicts no extras in run-001, as observed. Conclusion: C2 holds when d is
well below q/6, and C3 describes the regime beyond that.

## P2.5c Literature and novelty check for C3 (bounded)

- Four web searches with phrasing variants: degree of R with x(mP) = R(x(P)) on a subgroup;
  Fricke-function relations on X₁(q); Lange–Winterhof interpolation theorems.
- No statement of C3 found. The closest candidates are Lange–Winterhof (AAECC 2003, LNCS
  2643), "On degrees of polynomial interpolations related to ECC" (WCC 2005), and gonality
  bounds for modular curves (Abramovich; Derickx–Stoll).
- All were seen only as bibliographic records or abstracts (network policy), and the same
  agent did both passes.
- **Classification: apparently new within this narrow scope.** It is not claimed as new.

## P2.6 Clean re-runs and determinism (runs 005–007)

Runs 002–004 were re-run with the fixed script as runs 005 (threshold), 006 (C3 fit) and 007
(C3 held-out, including the CM side check in JSON).
- `summary.md` and `stdout.txt` are **byte-identical** to the originals for all three pairs, so
  the computation is deterministic.
- All three manifests validate against the current inputs.
- Runs 002–004 are kept as history; the validator reports them as stale.
- **Current evidence: run-001 (C2), run-005 (threshold), run-006 (C3 fit), run-007 (C3
  held-out).**

## P2.7 Conclusions of Phase 2 (direction 1)

1. **C2 holds in the tested regime.** Blind search over all multipliers recovers End(E)
   exactly, including the conductor of non-maximal orders, on 17/17 curves for d ≤ 20 (42).
2. **The heuristic behind C2 is wrong for small subgroups.** Curve-independent,
   non-endomorphism relations appear at degrees ≈ q/6 and above.
3. **C3** gives the exact minimal degree of every multiplier on generic curves: 908/908 on the
   fit set and 1454/1454 held out. It is unproved, its mechanism is unknown, and its novelty is
   unchecked beyond a narrow search.
4. **Mistakes caught, none affecting reported results:**
   - an inverted q² filter (stalled run, discarded);
   - a background job started in the wrong directory;
   - a silently failed `str.replace` that dropped a JSON field;
   - a closed-form comparison that would have ignored a scalar factor (fixed before any run).

---

# Post-Phase-2 reflection checks (exploratory, after network allowlist update)

**Network.** The new allowlist covers package registries and code hosts: PyPI, conda, GitHub,
GitLab and similar. Academic hosts (arxiv.org, link.springer.com, zbmath.org, lmfdb.org,
university sites) are still refused, so literature checks remain at search-result level.

**R1. j-table independently verified.** The 13 class-number-one j-invariants in `src/cm.py`
agree with:
- PARI/GP `polclass(D)` (cypari2 from PyPI): x − j(D), and `qfbclassno(D)` = 1;
- Wolfram `1728·KleinInvariantJ(τ_D)` at 40 digits.

This upgrades a load-bearing input from snippet level to independently computed.

**R2. C3 holds in characteristic 0 (numerical).** Script: `experiments/explore_c3_char0.py`.
- Setup: x_k = ℘(k/q; [1, τ]) at 130-digit precision for two generic τ and q = 19, 31, 37,
  43 (primes not in the fit or held-out sets). Relation existence was decided by the
  singular-value ratio σ_min/σ_max < 10⁻⁶⁰.
- Result: **C3 matches 63/63 multipliers per lattice**, 17 of them below threshold.
- Separation is large. For q = 37, m = −1/3 the ratio is 1.6·10⁻¹¹ at d = 5 and 4.6·10⁻¹²³ at
  d = 6, the C3 value.
- Bug caught: the first version reported "no relation" at d = d_int. When 2d+2 > M the matrix
  is wider than tall, so it is rank-deficient automatically, but SVD returns no zero singular
  value. Fixed.
- Status: conditional on the numerical tolerance. Not a recorded evidence run.

**R3. The nodal limit follows a different law (exact, F_ℓ with ℓ ≡ 1 mod q).** Script:
`experiments/explore_c3_torus.py`.
- Setup: x_k = ζ^k + ζ^{−k}, the Chebyshev analogue.
- Result: integer multipliers have degree m (not m²), and C3 fails for most fractions.
- So C3 is genuinely elliptic: the a² term is deg[a] = a².
- In every difference observed, the nodal degree is ≤ C3, consistent with relations
  persisting under degeneration (rank can only drop at a special fibre).

**Interpretation.** R2 together with the curve independence of P2.4 strongly supports the
hypothesis that C3 is an identity among level-q Fricke functions ℘(k/q; τ) on X₁(q). This
remains unproved.

---

# Phase 3 — direction A: proving C3 for fixed q

Approved by the project lead: "proceed with direction 1, starting with q = 13", together with a
request to create a GitHub repository.
- **Repository:** the session's GitHub API token is invalid, so a new repository cannot be
  created from here. The project was put under local git (branch `main`), ready to push once the
  user creates an empty repository and it is attached.
- **Commits:** 229edff (Phases 1–2), 6dbc732 (q = 13 certificate), 69f43e6 (general prover).

## P3.1 Sanity check of the q = 13 predictions

C3 at q = 13 (M = 6, d_int = 3) predicts:
- d(±4) = 2, since 4 ≡ −1/3;
- d(±1) = 1;
- all others 3.

Observed on 5 generic F_p curves and 2 complex lattices: exact agreement.

## P3.2 Universal model and T13 (run-001, symbolic)

- **Model:** Tate normal form over ℚ(b, c). The numerator of x(6P) − x(7P) is one irreducible F
  of bidegree (7, 10), monic in b. Independent check: 165 F_p-points of F give exact order 13
  with `src/ec.py`.
- **Main identity:** det A_{4,2} ≡ 0 mod F. The other determinants (m = 2, 3, 5, 6) are non-zero.
- **Supporting items:** rank 5 (uniqueness); a non-zero 4×4 minor at d = 1; denominators only
  c, b−c, b−c²−c, meeting F only at the singular cubic; F | det over ℤ.
- **Exceptional loci (exact):**
  - m = 2: cusps only;
  - m = 3: j = 0 (4 points);
  - m = 5, 6: j = 1728 (6 points);
  - m = 4 at degree 1: exactly the 4 j = 0 points, which was **predicted before computing**.
- **Mechanism probe:** R for m = 4 has no zeros or poles at 2- or 3-division roots (negative).
- **Bugs caught:**
  - PARI's `gcd` over POLMOD coefficients returned 1 wrongly; replaced by `nffactor`.
  - PARI's predefined `y` had the wrong variable priority; replaced by `varlower`.
  - Evaluating j at a cusp divided by zero; cusps are now classified via Δ.
  - The integrality check misread a c-denominator as non-integral content, and compared PARI
    type strings incorrectly; both rewritten.

## P3.3 q = 19 and a scalable method

- **Symbolic method:** the 8×8 determinants over ℚ(b, c) did not finish (> 17 CPU-minutes,
  1.1 GB). Recorded as resource-bounded and inconclusive; it was not used as evidence.
- **Slice-Bézout certificates** (`experiments/phase3/slices.py`): non-vanishing from one slice;
  vanishing from enough slices that the slice points exceed deg F · deg G.
  - Regression: it reproduces T13.
  - Negative control: m = 2 is correctly refused.
  - Bug caught by the control: a `varlower` variable name collided with PARI's `t`, which made
    the total degree 0 and the vanishing certificate trivially true. Fixed and re-verified.
- **T19 (run-003):** see `findings_phase3.md`.

## P3.4 Exceptional points and the combined law

**F_p enumeration of X₁(q)(F_p)** (run-004):
- q = 13: p = 1009 reproduces the exact locus; p = 2003 has none, as expected, since there
  j = 0 and j = 1728 are supersingular.
- q = 19: p = 2017 gives j = 0 and j = 54000 (D = −12) exceptions only; p = 1009 gives none.

**Combined law** d_min = min(C3, N(α)), tested on CM curves for q ∈ {13, 19} and 8
discriminants each (run-005): **16/16 cases match on every multiplier**.
