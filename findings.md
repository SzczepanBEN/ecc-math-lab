# Findings — Phase 1

**Labels used in this document:**

- **[Theorem]** means established mathematics, with a reference.
- **[Verified]** means we checked it exactly on every relevant point of the toy curves. It is
  not a proof for other curves.
- **[Conjecture]** means our own unproven claim.
- **[Artifact]** means an apparent signal that we traced to our own methods.

## 1. Summary

The sequence G, 2G, 3G, … on an elliptic curve over F_p has **two very different faces**:

1. **Statistically, it is indistinguishable from a random ordering of the same points.** No
   tested statistic detects any dependence of the coordinates on the index k: autocorrelation,
   Ljung–Box, consecutive differences, 2-D serial test, rank correlation with k, initial-segment
   discrepancy, periodogram, and linear complexity. This holds on generic, j = 0 (ordinary and
   supersingular) and j = 1728 curves alike: 1 of 49 tests had p < 0.05, against 2.5 expected by
   chance.
2. **Algebraically, it is rigidly structured.** Every pair (x_k, x_{k+ℓ}), at every lag ℓ, lies
   on one fixed curve of bidegree (2, 2): Semaev's summation polynomial S3(X, Y, x(ℓG)). Points
   are paired by negation (x_{n−k} = x_k), and on j = 0 / j = 1728 curves they are additionally
   linked across distant indices by the automorphisms (x_{λk} = βx_k).

Both faces are expected from known theory. The useful scientific content of this phase is the
precise boundary between them, and a list of the ways naive analysis produced false positives.

## 2. Confirmed relationships (all known mathematics)

| # | relationship | status | explanation |
|---|---|---|---|
| R1 | x_{n−k} = x_k, y_{n−k} = −y_k | [Verified] 100% of k, 7 curves | −P = (x, −y) **[Theorem]** |
| R2 | S3(x_k, x_{k+ℓ}, x(ℓG)) = 0, minimal bidegree (2, 2) | [Verified] all pairs, 5–7 lags × 7 curves; found blindly by linear algebra and equal to S3 coefficient by coefficient | group law; Semaev 2004 **[Theorem]** |
| R3 | at ℓ = n/2: (x_{k+ℓ} − t)(x_k − t) = 3t² + a, bidegree (1, 1) | [Verified] GEN-C6, J0-SS, J1728 | translation by 2-torsion T = (t, 0) is a Möbius map on x **[Theorem]** (elementary) |
| R4 | the multipliers m with x_{mk} = Möbius(x_k) are exactly the automorphisms: 2 / 6 / 4 of them | [Verified] 31 curves; existence **[Theorem]**; "no others" **[Conjecture]** | Aut(E) for j ≠ 0, 1728 is {±1}; j = 0: μ6; j = 1728: μ4 (Silverman, AEC III.10); GLV 2001 |
| R5 | y-value occurrence counts 0/1/2/3 in proportions 1/3 : 1/2 : 1/6 (generic) | [Verified] χ² p = 0.39, 0.88, 0.96 | number of roots of the cubic f(x) − y²; Chebotarev density for Galois group S3 **[Theorem]**, heuristic for a single p |
| R6 | j = 0, p ≡ 2 (mod 3): every y ∈ F_p occurs exactly once; #E = p + 1 | [Verified] | x ↦ x³ is a bijection of F_p **[Theorem]**; hence supersingular, embedding degree 2 (MOV / Frey–Rück) |
| R7 | j = 0 ordinary: y occurs 0 or 3 times, x³ occurs 0 or 6 times | [Verified] | orbits of the order-6 automorphism group **[Theorem]** |
| R8 | the point with y = 0 sits at k = n/2; zeros of ψ3 at k = n/3, 2n/3 | [Verified] | torsion points; division polynomials **[Theorem]** |
| R9 | no low-degree polynomial F(k, x_k) = 0 at bidegrees ≤ (60, 1) or (10, 10); linear complexity ⌈N/2⌉ | [Verified] negative result, with positive controls detected | consistent with Lange–Winterhof interpolation bounds and Hess–Shparlinski EC-LCG bounds **[Theorem]** (our tested degrees are far weaker than those bounds) |

### Generic vs j = 0 curves

The **only** differences we observed between families are in the value multiset (R5–R7) and the
extra exact symmetries (R4). The order statistics of j = 0 curves are as random-looking as
those of generic curves. All structure specific to j = 0 is captured by the automorphism group
of order 6.

## 3. Why none of this is a computational shortcut

This section is interpretation, based on established theory.

- **R2 and R3 are "relative" relations.** They tell you which x-coordinate comes ℓ steps after a
  known point, *given* x(ℓG). That is exactly what computing P + ℓG does, which is already cheap.
  They do not map a point to its index. The proven route from summation polynomials to discrete
  logarithms is Semaev's index calculus, which requires solving multivariate polynomial systems.
  For prime-field curves it is not known to beat Pollard rho.
- **R1 and R4 reduce the problem only by the size of the automorphism group.** Working on
  equivalence classes {±P}, or {ζP} for j = 0, gives a theoretical Pollard-rho speedup of
  √|Aut|: √2 generally, √6 for j = 0 (Duursma–Gaudry–Morain 1999; Wiener–Zuccherato 1998).
  In practice the gain is smaller because of fruitless cycles (Bos–Kleinjung–Lenstra 2010;
  Bernstein–Lange–Schwabe 2011). This is a constant factor and is already accounted for in
  standard security estimates, including for secp256k1.
- **Statistical regularity of value multisets (R5–R7) concerns *which* points exist, not
  *where* in the sequence they are.** Over a full period the sequence is just an enumeration of
  the point set.
- **The structure that does give real attacks is detected by invariants, not statistics.**
  Supersingular and low-embedding-degree curves (MOV / Frey–Rück) and anomalous curves (Smart;
  Satoh–Araki; Semaev) are flagged by `src/curvegen.classify` from (p, #E) alone. J0-SS is such
  a curve, and its x-sequence looked exactly as random as the others in Experiment A. **The
  coordinate sequence did not reveal the one curve in our set that is cryptographically weak
  by known theory.** This is direct evidence that coordinate-pattern analysis is the wrong lens
  for finding weaknesses.

## 4. A methodological finding: statistics are blind to algebra over F_p

Figure A2 shows (x_k, x_{k+1}) next to a shuffled version. The two are visually and
statistically indistinguishable: 2-D serial χ² p-values range from 0.08 to 0.75. Yet every
point in the left panel lies on one degree-(2,2) curve, which a 9-column linear system finds
immediately. The reason is that a curve of bounded degree over F_p has about p points, which
(by Weil-type bounds) are equidistributed in [0, p)². Embedding F_p into the reals for plotting
or correlation destroys the field structure the relation lives in.

**Implication for AI-assisted discovery** (our interpretation): pattern-finding on real-valued
embeddings of coordinates (plots, correlations, generic ML features) is the wrong inductive bias
for this domain. Hypothesis spaces native to F_p are both more powerful and more honest,
because they come with exact verification. Examples: polynomial relations via nullspaces,
cross-ratio invariants, linear complexity over F_p. In Phase 1 the native methods re-discovered
the group law (S3), 2-torsion translation and the automorphism groups *blindly*, and they found
nothing between k and x(kG).

## 5. Unexpected observations, all resolved

1. **[Artifact] y-spread looked too extreme on several curves.** The cause was a null model that
   ignored the y ↔ −y symmetry of the occupancy; it underestimated the null variance by about 2×.
   Fixed with a negation-commuting relabelling null.
2. **[Artifact] curves seemed to spread x more evenly than random sets** (KS p = 0.0006). The
   cause was an atypical 300-curve null sample, roughly a 3σ outlier: none of 120 replicate
   seeds was as extreme. A 36,000-curve pooled null matches χ²₆₃ (mean 62.90, variance 124.3).
   After this, the curve-level and relabelling nulls agree on 3 primes × 2 bin counts
   (|z| ≤ 1.14).
3. **[Artifact] the plug-in conditional entropy H(x_{k+1} | x_k) ≈ 1 bit** suggested near-
   determinism. It is identical for shuffled sequences, because each x occurs exactly twice.
4. **[Artifact] "hollow points aligned on round numbers"** in a scatter plot came from gridlines
   drawn over the markers.
5. **The single nominal ordering hit (J0-B, p = 0.017)** did not reproduce under 30 alternative
   generators (1/30 below 0.05). It is consistent with multiple testing.
6. **Mildly interesting, not anomalous:** the linear-complexity profile is *exactly* ⌈N/2⌉ at
   every N for every sequence, real or shuffled. This is expected over a large field: a profile
   deviation needs a coincidence of probability about 1/p. But it means this classical test
   cannot discriminate among random-looking sequences at this field size.

**Nothing in this phase qualifies as a new mathematical observation.** We checked the
literature and do not claim novelty for any result.

## 6. Limitations

- **Small fields:** p ≈ 2^16–2^17 (2^13–2^14 for the H2 generalisation set). Asymptotic
  statements rest on the cited theory, not on our data.
- **Few curves:** 7 main curves and 24 auxiliary ones. All are cyclic, and only one has a
  composite order with small factors (GEN-C6). Non-cyclic groups were tested only in unit tests.
- **Monte-Carlo resolution:** most nulls use 200 permutations (p-value resolution about 0.005).
  The tests on one curve are correlated, so the QQ plot is indicative rather than exact.
- **Weak index-relation test:** H3 tests polynomial degrees ≤ 60, far below the proven
  interpolation lower bounds. A negative result there is a sanity check, not strong evidence.
- **R4's "only automorphisms" part is not proved.** It holds for all 31 tested curves.
- **Real-valued embedding:** statistical tests use x/p ∈ [0, 1). Other embeddings (bits,
  characters of F_p) were not explored, except implicitly through the F_p-native tests.
- **Implementation:** affine arithmetic in pure Python. It is verified against independent
  algorithms and public test vectors, but it is not constant-time and not meant for real use.

## References

- I. Semaev, *Summation polynomials and the discrete logarithm problem on elliptic curves*, IACR
  ePrint 2004/031. https://eprint.iacr.org/2004/031
- F. Hess, I. Shparlinski, *On the linear complexity and multidimensional distribution of
  congruential generators over elliptic curves*, Des. Codes Cryptogr. 35 (2005) 111–117.
  https://research-information.bris.ac.uk/en/publications/on-the-linear-complexity-and-multidimensional-distribution-of-con/
- D. Kohel, I. Shparlinski, *On exponential sums and group generators for elliptic curves over
  finite fields*, ANTS-IV (2000). https://www.i2m.univ-amu.fr/perso/david.kohel/pub/character.pdf
- T. Lange, A. Winterhof, *Polynomial interpolation of the elliptic curve and XTR discrete
  logarithm*, COCOON 2002, LNCS 2387, 137–143. https://www.hyperelliptic.org/tanja/publications.html
- R. Gallant, R. Lambert, S. Vanstone, *Faster point multiplication on elliptic curves with
  efficient endomorphisms*, CRYPTO 2001, 190–200 (tech. report CORR 2000-53).
  https://uwaterloo.ca/combinatorics-and-optimization/sites/default/files/uploads/documents/2000_53.pdf
- I. Duursma, P. Gaudry, F. Morain, *Speeding up the discrete log computation on curves with
  automorphisms*, ASIACRYPT 1999, LNCS 1716, 103–121. https://link.springer.com/doi/10.1007/3-540-45539-6_2
- J. Bos, T. Kleinjung, A. Lenstra, *On the use of the negation map in the Pollard rho method*,
  ANTS IX (2010), LNCS 6197, 66–82. https://infoscience.epfl.ch/record/164553
- D. Bernstein, T. Lange, P. Schwabe, *On the correct use of the negation map in the Pollard rho
  method*, PKC 2011 (ePrint 2011/003). https://eprint.iacr.org/2011/003
- J. Silverman, *The Arithmetic of Elliptic Curves*, 2nd ed., Springer GTM 106: automorphism
  groups (III.10), division polynomials (Ex. 3.7), Hasse bound (V.1). Standard reference, not
  re-checked online.
- M. Wiener, R. Zuccherato, *Faster attacks on elliptic curve cryptosystems*, SAC 1998. Cited
  from memory; not re-checked online in this session.
