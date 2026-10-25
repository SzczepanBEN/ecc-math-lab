# Findings — Phase 3, direction A: making conjecture C3 rigorous (q = 13, q = 19)

**Labels used in this document:**

- **[Theorem]**: proved here by an exact computation whose certificate is in `results/phase3/`.
  The labels (E), (U) and (L) refer to the three parts listed under "Setting".
- **[Standard]**: a textbook fact used without proof.
- **[Evidence]**: finite tests that are not proofs.
- **[Conjecture]**: our own claim.

## Setting

Every elliptic curve with a point P of order ≥ 4 has a Tate normal form [Standard]:
y² + (1−c)xy − by = x³ − bx², with P = (0,0).

For a prime q, the pairs (E, P) with P of order q are the points of a plane model of X₁(q),
given by F_q(b, c) = 0, where F_q is the numerator of x(hP) − x((h+1)P) and h = (q−1)/2. We
work in its function field K. For multipliers m and degrees d, a relation
Yψ(X) + φ(X) = 0 (deg ≤ d) on the M = h pairs (x(kP), x(mkP)) exists iff the M × (2d+2)
matrix A_{m,d} has a non-zero kernel vector.

Each theorem below has three parts:
- **(E) existence:** the relation exists at the C3 degree;
- **(U) uniqueness:** it is unique up to scaling;
- **(L) lower bound:** no relation exists one degree lower, generically.

## T13 — C3 at q = 13 [Theorem] (runs 001 and 002, two independent methods)

The model: F₁₃ is irreducible over ℚ, of bidegree (7, 10) and monic in b. 165 random
F_p-points were independently checked to give exact order 13.

1. **(E, U) for m = ±4 ≡ −1/3:** det A_{4,2} = 0 in K, and the rank is exactly 5. So there is a
   unique rational map R of degree ≤ 2 with x(4Q) = R(x(Q)) for all Q ∈ ⟨P⟩.
2. **(L) for m = ±4:** generically R has degree exactly 2; some 4×4 minor of A_{4,1} is non-zero.
3. **(L) for m ∈ {±2, ±3, ±5, ±6}:** generically no relation of degree ≤ 2 exists, so the minimal
   degree is 3, where linear algebra forces one.
4. **Every field:**
   - the x(kP) have denominators c, b−c and b−c²−c only, with integral numerators;
   - these meet F₁₃ only at c = 0, where F₁₃(b, 0) = b⁷ forces the singular cubic;
   - the resultants are ±c¹¹ and ±c¹⁴, with unit constants;
   - F₁₃ divides det A_{4,2} over ℤ.

   Hence **(1) holds for every elliptic curve with a point of order 13 over every field**.
5. **The exceptional curves, determined exactly:**

   | statement | fails exactly at | explanation |
   |---|---|---|
   | m = ±2 needs degree 3 | cusps only (6 cusps) | holds for *every* elliptic curve |
   | m = ±3 needs degree 3 | the 4 points with j = 0 on ζ-eigenlines | 3² + 3 + 1 = 13: 3 is an eigenvalue of the automorphism ζ |
   | m = ±4 needs degree 2 | the same 4 points, j = 0 | 4 ≡ −ζ²: degree 1 there |
   | m = ±5 needs degree 3 | the 6 points with j = 1728 on i-eigenlines | 5² ≡ −1: the automorphism i |
   | m = ±6 needs degree 3 | the same 6 points, j = 1728 | 6 = 1 + 5: the endomorphism 1 + i, norm 2 |

   The counts 4 and 6 equal the numbers of eigenline points of X₁(13) above j = 0 and j = 1728.
   There are 2 eigenlines × 12 primitive points = 24 points, divided by |Aut(E)|: 24/6 = 4 and
   24/4 = 6.

The two methods were symbolic determinants over ℚ(b, c)/(F) (run-001) and Bézout slice
certificates (run-002). They agree.

## T19 — C3 at q = 19 [Theorem, generic statements] (run-003)

The model: F₁₉ is irreducible over ℚ, of bidegree (15, 22) and total degree 24, monic in b.
76 random F_p-points were independently checked to give exact order 19. C3 predicts d(±6) = 3
(6 ≡ −1/3) and d = 4 = d_int for every other m ≠ ±1.

1. **(E) for m = ±6:** all nine 8×8 minors of A_{6,3} vanish in K. Each is certified by Bézout:
   at least 391 slices × 15 points against bounds of at most 6816.
2. **(U):** the rank is exactly 7, so the relation is unique.
3. **(L):** for m = ±6, no relation of degree ≤ 2 generically. For m ∈ {±2, ±3, ±4, ±5, ±7,
   ±8, ±9}, no relation of degree ≤ 3 generically, so these need degree 4, where linear algebra
   forces one.
4. **Every field:**
   - the denominators of x(kP) are c, b−c, b−c²−c, b²−cb−c³, 2b²−(c²+3c)b+c² and
     b³−3cb²+(c³+3c²)b−(c⁵+c⁴+c³);
   - every resultant with F is ±c^k with a unit constant;
   - F(b, 0) = b¹⁵;
   - the x(kP) are integral over these denominators (exploratory check, not in the run);
   - F | minor over ℚ implies divisibility of the primitive part over ℤ (Gauss's lemma).

   Hence **(1) holds for every elliptic curve with a point of order 19 over every field.**
5. **Exceptional points:** the exact locus for q = 19 was not computed, because the resultants
   have degree in the thousands. Evidence on it comes from complete F_p enumerations and CM
   tests (below).

Runtime is 107 s. The direct symbolic method did not finish even its first 8×8 minor in about
19 CPU-minutes.

## The combined law on CM curves [Conjecture, with evidence]

> **C2+C3.** For E with End(E) = O and P of prime order q:
> d_min(m) = min( C3(m), min{ N(α) : α ∈ O acts as ±m on ⟨P⟩ } ).

The second term applies only when ⟨P⟩ is End-stable, i.e. when q splits in O and q ∥ #E.

**Evidence:**
- **(i) Run-005:** one CM curve per (q, D), for q ∈ {13, 19} and the 8 discriminants
  D = −3 … −19. All **16/16 cases match on every multiplier**, including:
  - j = 8000 at q = 19, where m = 6 has degree 2 (√−2) and m = 5, 7 have degree 3 (1 ± √−2);
  - inert controls, which never show extras.
- **(ii)** The exact q = 13 exceptional locus above is exactly this law.
- **(iii) Run-004:** complete enumerations of X₁(q)(F_p) for (q, p) = (13, 1009), (13, 2003),
  (19, 1009), (19, 2017). Every exceptional point has a CM j-invariant predicted by the law:
  j = 0 or 54000 at q = 19; j = 0 or 1728 at q = 13.

**Structural reason for the inert controls [Standard + elementary]:** if q is inert in the CM
field, q | N(π − 1) forces q² | #E, so no rational q-torsion line is End-stable.

## What is still open

- **The mechanism of C3.**
  - The q = 13 map R for m = −1/3 is a generic type-(2,2) rational map. It has no zeros or
    poles at 2- or 3-torsion x-coordinates.
  - The nodal (Chebyshev) degeneration obeys a different law.
  - C3 is therefore genuinely elliptic and is not explained by torsion structure.
- **C3 for general q.** Each prime now needs only a finite slice computation, but the Bézout
  bounds grow quickly with q. No uniform proof exists.
- **Novelty.** This is still unchecked against primary literature, because academic hosts
  remain blocked.

## Method lessons

- **Symbolic determinants over ℚ(b, c) do not scale.** The q = 19 attempt ran for over
  17 CPU-minutes and 1.1 GB without finishing.
- **Bézout slice certificates scale.** They prove vanishing in K from univariate computations,
  and a single slice proves non-vanishing. q = 19 takes minutes.
- **Bugs caught by negative controls and asserts:**
  - a variable-name collision made a degree bound 0, which made the certificate trivially pass;
  - PARI type checks compared string objects incorrectly;
  - PARI's predefined `y` had the wrong variable priority.

  None reached a reported result.
