"""Elliptic curve arithmetic over prime fields F_p (short Weierstrass form).

    E : y^2 = x^3 + a*x + b   over F_p, p > 3 prime, 4a^3 + 27b^2 != 0 (mod p).

Points are tuples (x, y) of Python ints in [0, p); the point at infinity
(the group identity O) is represented by ``None``.

Everything here is written for clarity and verifiability on *toy* parameters
(p up to ~2^20 for the counting routines). Arithmetic itself is exact for any
size of p, which the tests use to cross-check against published secp256k1
constants. Nothing here is constant-time or suitable for real cryptography.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

INF = None  # point at infinity


# --------------------------------------------------------------------------
# number theory helpers
# --------------------------------------------------------------------------

def is_probable_prime(n: int, rounds: int = 32) -> bool:
    """Deterministic for n < 3.3e24 (fixed Miller-Rabin bases), probabilistic beyond."""
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    for q in small:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factorize(n: int) -> dict[int, int]:
    """Trial division (fine for the toy group orders used here, n < 2^40)."""
    f: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def legendre(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def sqrt_mod(a: int, p: int) -> int | None:
    """Tonelli-Shanks square root modulo an odd prime. Returns None if a is a non-residue."""
    a %= p
    if a == 0:
        return 0
    if legendre(a, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while legendre(z, p) != -1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def cube_roots_of_unity(p: int) -> list[int]:
    """Non-trivial cube roots of 1 mod p (empty unless p = 1 mod 3)."""
    if p % 3 != 1:
        return []
    s = sqrt_mod(-3 % p, p)
    inv2 = pow(2, -1, p)
    return sorted({(-1 + s) * inv2 % p, (-1 - s) * inv2 % p})


# --------------------------------------------------------------------------
# curve
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Curve:
    p: int
    a: int
    b: int
    name: str = field(default="", compare=False)

    def __post_init__(self):
        if self.p <= 3 or not is_probable_prime(self.p):
            raise ValueError(f"p={self.p} must be a prime > 3")
        object.__setattr__(self, "a", self.a % self.p)
        object.__setattr__(self, "b", self.b % self.p)
        if self.discriminant == 0:
            raise ValueError("singular curve (4a^3 + 27b^2 = 0 mod p)")

    # ----- invariants -----
    @property
    def discriminant(self) -> int:
        p = self.p
        return (-16 * (4 * pow(self.a, 3, p) + 27 * self.b * self.b)) % p

    @property
    def j_invariant(self) -> int:
        p = self.p
        num = 1728 * 4 * pow(self.a, 3, p)
        den = (4 * pow(self.a, 3, p) + 27 * self.b * self.b) % p
        return num * pow(den, -1, p) % p

    def rhs(self, x: int) -> int:
        return (x * x * x + self.a * x + self.b) % self.p

    # ----- point predicates -----
    def is_on_curve(self, P) -> bool:
        if P is INF:
            return True
        x, y = P
        if not (0 <= x < self.p and 0 <= y < self.p):
            return False
        return (y * y - self.rhs(x)) % self.p == 0

    def lift_x(self, x: int):
        """Return a point with x-coordinate x (smaller of the two y's) or None."""
        y = sqrt_mod(self.rhs(x), self.p)
        if y is None:
            return None
        return (x % self.p, min(y, (-y) % self.p))

    def random_point(self, rng: random.Random):
        while True:
            P = self.lift_x(rng.randrange(self.p))
            if P is not None:
                if rng.random() < 0.5 and P[1] != 0:
                    P = (P[0], self.p - P[1])
                return P

    # ----- group law (affine chord-and-tangent) -----
    def neg(self, P):
        if P is INF:
            return INF
        return (P[0], (-P[1]) % self.p)

    def add(self, P, Q):
        if P is INF:
            return Q
        if Q is INF:
            return P
        p = self.p
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return INF  # Q = -P (covers 2-torsion: y1 = y2 = 0)
            lam = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    def double(self, P):
        return self.add(P, P)

    def sub(self, P, Q):
        return self.add(P, self.neg(Q))

    def mul(self, k: int, P):
        """Scalar multiplication by left-to-right double-and-add. Handles k <= 0."""
        if k < 0:
            return self.mul(-k, self.neg(P))
        R = INF
        for bit in bin(k)[2:] if k else "":
            R = self.add(R, R)
            if bit == "1":
                R = self.add(R, P)
        return R

    def mul_ladder(self, k: int, P):
        """Montgomery ladder; an independent implementation used to cross-check mul()."""
        if k < 0:
            return self.mul_ladder(-k, self.neg(P))
        R0, R1 = INF, P
        for bit in bin(k)[2:] if k else "":
            if bit == "1":
                R0, R1 = self.add(R0, R1), self.add(R1, R1)
            else:
                R0, R1 = self.add(R0, R0), self.add(R0, R1)
        return R0

    def multiples(self, G, count: int) -> list:
        """[G, 2G, ..., count*G] by repeated addition (the EC 'linear congruential' sequence)."""
        out, P = [], INF
        for _ in range(count):
            P = self.add(P, G)
            out.append(P)
        return out

    # ----- counting & structure (toy sizes only) -----
    @cached_property
    def _square_table(self) -> np.ndarray:
        p = self.p
        if p > 1 << 24:
            raise ValueError("table-based counting is for toy primes only")
        ys = np.arange(p, dtype=np.int64)
        cnt = np.zeros(p, dtype=np.int64)
        np.add.at(cnt, (ys * ys) % p, 1)  # cnt[v] = #{y : y^2 = v}
        return cnt

    def _rhs_array(self) -> np.ndarray:
        p = self.p
        x = np.arange(p, dtype=np.int64)
        return ((x * x % p) * x % p + self.a * x + self.b) % p

    def order(self) -> int:
        """#E(F_p) = 1 + sum_x #{y : y^2 = x^3+ax+b} (exact; O(p) with numpy)."""
        return 1 + int(self._square_table[self._rhs_array()].sum())

    def order_legendre(self) -> int:
        """Independent count: p + 1 + sum_x legendre(x^3+ax+b, p) (pure Python, slow)."""
        return self.p + 1 + sum(legendre(self.rhs(x), self.p) for x in range(self.p))

    def trace(self) -> int:
        return self.p + 1 - self.order()

    def all_points(self) -> list:
        pts = [INF]
        for x in range(self.p):
            P = self.lift_x(x)
            if P is None:
                continue
            pts.append(P)
            if P[1] != 0:
                pts.append(self.neg(P))
        return pts

    def point_order(self, P, N: int | None = None) -> int:
        """Exact order of P given a multiple N of it (default: the group order)."""
        if N is None:
            N = self.order()
        if self.mul(N, P) is not INF:
            raise ArithmeticError("N is not a multiple of ord(P)")
        m = N
        for q in factorize(N):
            while m % q == 0 and self.mul(m // q, P) is INF:
                m //= q
        return m

    def group_structure(self, rng: random.Random | None = None, tries: int = 40) -> tuple[int, int]:
        """Return (n1, n2) with E(F_p) = Z/n1 x Z/n2, n2 | n1, n2 | p-1.

        n1 (the group exponent) is estimated as the lcm of orders of random points; the
        result is certified when some sampled point has order n1 and N/n1 passes the
        divisibility constraints. For tiny p the tests check this against brute force.
        """
        rng = rng or random.Random(0)
        N = self.order()
        e = 1
        for _ in range(tries):
            e = math.lcm(e, self.point_order(self.random_point(rng), N))
            if e == N:
                break
        n2 = N // e
        assert e % n2 == 0 and (self.p - 1) % n2 == 0, "inconsistent group structure"
        return e, n2


def find_generator(E: Curve, rng: random.Random):
    """A point of order #E (requires E(F_p) cyclic); certified via point_order."""
    N = E.order()
    for _ in range(10_000):
        G = E.random_point(rng)
        if E.point_order(G, N) == N:
            return G
    raise ValueError("no generator found; group is probably not cyclic")
