# Findings — Phase 2, direction 1 (blind discovery of small-degree endomorphisms)

**Labels used in this document:**

- **[Theorem]** means established mathematics; its literature status is given in
  `research_log.md` P2.0, and was checked at search-result level only.
- **[Finite, exact]** means a statement about the listed curves and degrees, decided exactly by
  computation over F_p with no sampling error.
- **[Conjecture]** means our own claim, with the evidence stated.
- **[Speculation]** means a suggested explanation that has not been tested.

## 1. Summary

On every curve tested, a blind search over all multipliers m recovers the
**endomorphism ring** as long as the degree d is small relative to the subgroup order q. The
search asks, for each m, whether x(mP) is a rational function of x(P) of degree ≤ d. The
multipliers it finds at minimal degree d are exactly {α mod q : α ∈ End(E), N(α) = d}, and the
search distinguishes End = O from End = O_K and from Z[π].

When d is comparable to q, a **second, previously unanticipated family** of low-degree relations
appears. These relations do not come from endomorphisms and do not depend on the curve, and
their degrees follow an exact formula in (q, m), conjecture C3. That formula predicted 1454 of
1454 held-out multipliers exactly.

## 2. Endomorphisms are recovered exactly (claim C2, run-001)

[Finite, exact] The setting:
- 13 curves, one for each imaginary quadratic order of class number 1 (j = j(O) mod p);
- 4 generic controls;
- p ∈ [2^13, 2^14), subgroups of prime order q ∈ [1097, 11831];
- every m ∈ [1, q−1] and every d ≤ 20 (≤ 42 for D = −163).

| check | result |
|---|---|
| observed multipliers of minimal degree d = S_d(O) for all d | **17/17 curves** |
| End = O_K rejected for the non-maximal orders D = −12, −16, −27, −28 | 4/4 |
| End = Z[π] rejected on CM curves | 13/13 |
| order discriminant inferred blindly (trace T from m² − Tm + d ≡ 0, conductor gcd) | equals the table order, 13/13 |
| CM field of t² − 4p equals the table field (checks the j-table) | 13/13 |
| found maps have deg num = d, deg den = d − 1, den = s·h² with s \| x³+ax+b (Vélu shape) | all |
| m = 2 reproduces the closed-form doubling formula | 17/17 |
| batched screen vs independent per-m rank test | 0 disagreements / 2,404 |

Examples:
- **D = −7:** multipliers at d = 1, 2, 4, 7, 8, … with counts 2, 4, 6, 2, 8, ….
  These are r_O(d) for O = Z[(1+√−7)/2].
- **D = −163:** only ±n for n = 1, …, 6 (degrees n²) up to d = 40, then 4 multipliers at
  d = 41 = N((1 ± √−163)/2).

**Interpretation:**
- The inclusion "S_d has degree ≤ d" is a [Theorem], via the eigenvalue argument and the
  isogeny standard form.
- That *nothing else* appears (C2) is verified here as a [Finite, exact] statement for these
  17 curves and degrees.
- The identification of the observed sets with End(E) also uses Deuring's reduction theorem.
  The data supports that independently by discriminating between the orders.
- **Relation to prior work:** computing End(E) is a solved problem with better algorithms
  (Kohel 1996; Bisson–Sutherland 2011). Our method is a naive, *blind* detector of small-norm
  endomorphisms. Its interest is that it needs no knowledge of CM theory to find them.

## 3. The unexpected part: curve-independent non-endomorphism relations (C3)

### What was predicted

A relation Yψ(X) + φ(X) = 0 of bidegree (d, 1) has 2d + 2 coefficients. The subgroup gives
M = (q−1)/2 distinct pairs (x_k, x_{mk}). The naive expectation was:
- every m has a relation once 2d + 2 > M (true: [Theorem], linear algebra);
- below that, only endomorphisms have relations.

### What was observed (runs 005, 006, 007; first produced as 002–004, byte-identical)

The second part fails for small q. Non-endomorphism multipliers appear well below the
threshold d_int = ⌊M/2⌋ (smallest d with 2d+2 > M). Examples:
- **m ≡ ±1/3 (mod q):** these appear at d = M/3 whenever q ≡ 1 (mod 6). For q = 61 this is
  d = 10 against d_int = 15.
- **Other fractions:** degrees fall below the threshold for many m ≡ a/b with small a and odd b.

Several features stand out:
- **Curve independence:** for each of 10 primes q, two different generic curves over
  *different* fields F_p have **identical** d_min profiles for all m.
- **Not explained by endomorphisms:** the CM discriminants are far too large (|D_K| > 4q).
- **Real relations:** independent rank computations agree, and the maps are in lowest terms.

### Conjecture C3 [Conjecture]

For E/F_p with no non-integer endomorphism of norm < q, and G′ of prime order q (q ∤ N/q),
let M = (q−1)/2 and d_int = ⌊M/2⌋. Then the minimal degree of a rational R with
x(mP) = R(x(P)) on ⟨G′⟩ is

  d_min(m) = min( d_int, min over (a, b) of ((b−1)q + 4a² − ε_b)/(4b) ),

where:
- (a, b) runs over coprime pairs with b odd, bm ≡ ±a (mod q), and q ≡ 4a² (mod b), which is
  exactly when the numerator is divisible by 4b;
- ε_b = 0 if b ≡ 1 (mod 4) and ε_b = 2b if b ≡ 3 (mod 4).

Equivalently, with c = (4a² − q)/b, the degree is d = (q + c)/4, minus ½ when b ≡ 3 (mod 4).
The saving below the trivial threshold, about (q − 4a²)/(4b), is large only for small a and b.

Special cases:
- b = 1 gives d = m², the integer multipliers [m].
- Even denominators never beat the threshold.
- For b = 3, a = ±1, q ≡ 1 (mod 6): d = M/3.

**Evidence:**

| data set | multipliers | exact matches | relations below threshold |
|---|---|---|---|
| fit: 10 primes × 2 curves (run-006) | 908 | 908 | 194 |
| **held-out, frozen formula**: 20 new primes, new seeds (run-007) | **1454** | **1454** | **253** |

The formula was written down after inspecting about 15 cases from the fit set, then frozen in
code (`conjectured_dmin`) before the held-out run.

**CM side check (not part of C3):** on CM curves with small q (D = −3, −4, −7), the observed
d_min is never higher than the C3 value but is often lower (15/21, 15/26 and 26/39
multipliers). This is consistent with C3's relations existing on *every* curve, with
endomorphisms adding further ones.

**Consistency with Section 2:** for every b ≥ 3, the C3 degree is at least (q − 3)/6. That is
≥ 182 for every curve in run-001, which had d ≤ 42. So C3 predicts no extras there, which is
what was observed. C2 is therefore valid when d is well below q/6, and C3 describes what happens
beyond.

### Possible explanation [Speculation]

Curve independence suggests the relations are identities on the **universal** elliptic curve
over the modular curve X₁(q). On that family, x(kG′) are Fricke/Weierstrass functions of level
q, and the coefficients of the relation would be modular functions on X₁(q).

The data hints at a mechanism through cusps: cusps of X₁(q) are labelled by fractions, and orders
of Fricke functions at cusps are given by B₂-type quadratic expressions such as a²/b. If this
is right, C3 should hold in characteristic 0 and should be provable by counting zeros and poles
on X₁(q).

This explanation has **not** been tested. It is the most promising route to a proof.

### Novelty [bounded]

C3 is **apparently not stated in the literature within our search scope**, but that scope was
narrow:
- a few web searches with phrasing variants;
- the closest candidates (Lange–Winterhof, AAECC 2003; "On degrees of polynomial interpolations
  related to elliptic curve cryptography", WCC 2005; gonality bounds for X₁(q)) were seen only
  as bibliographic records, because the network policy blocks primary texts;
- the same agent performed the discovery and the verification passes.

It should be treated as a candidate observation until those papers, and literature on Fricke
functions and modular units, are checked.

## 4. What this does and does not mean

- **It is not a cryptographic shortcut.** Both C2 and C3 concern rational functions that move
  *along* a subgroup: they compute mP from P, which is already easy. They say nothing about
  recovering m from P and mP.
- C3's degrees are about q/6 or more, so on cryptographic-size curves they are astronomically
  large. The structure is mathematically interesting, not computationally useful.
- C3 is a statement about small finite sets: rational interpolation of a permutation of
  M = (q−1)/2 points. Its regime is q ≤ 199, p < 2^14.

## 5. Limitations

- **C2:**
  - verified for 17 curves at one field size, with d ≤ 20 (42);
  - the "nothing else" direction is not proved in general;
  - the parameter-count heuristic behind it turned out to be wrong in the small-q regime.
- **C3:**
  - verified on 30 generic curves (20 held out), q ∈ [59, 199], p ∈ [2^13, 2^14);
  - only q exactly dividing N with |D_K| > 4q;
  - only bidegree (d, 1) relations;
  - not proved, and no mechanism has been confirmed.
- **Literature:** all checks were at search-result level (network policy). Deuring's theorem,
  the j-table and the Vélu shape were also supported by internal consistency checks; C3's
  novelty was not.
- **Provenance:** there is no git repository, so manifests carry input hashes but no commit
  (`"commit": "unavailable"`).
