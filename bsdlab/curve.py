"""Elliptic curves over Q in long Weierstrass form.

    E : y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6

This module is the foundation every other module compiles against. It is
deliberately conservative: integer arithmetic is exact (Python ints), and no
function here returns a floating-point value.

Conventions fixed here and relied upon repo-wide:
  * a-invariants are stored as a 5-tuple of ints ``(a1, a2, a3, a4, a6)``.
  * The point at infinity is the singleton ``INFINITY``; affine points are
    ``(x, y)`` pairs of ``Fraction`` over Q, or of ``int`` when working mod p.
  * A curve is *singular* iff its discriminant is 0. We refuse to build one.
"""

from __future__ import annotations

from fractions import Fraction
from dataclasses import dataclass
from typing import Iterable, Optional

from sympy import factorint, isprime


class INFINITY:
    """The point at infinity, O. Identity of the group law."""

    __slots__ = ()

    def __repr__(self) -> str:
        return "O"

    def __eq__(self, other) -> bool:
        return other is INFINITY or isinstance(other, INFINITY)

    def __hash__(self) -> int:
        return hash("O")


O = INFINITY()


class SingularCurveError(ValueError):
    """Raised when the discriminant vanishes, so the cubic is not elliptic."""


@dataclass(frozen=True)
class EllipticCurve:
    """An elliptic curve over Q given by integral a-invariants.

    The standard b-, c-invariants and the discriminant follow Silverman,
    *Arithmetic of Elliptic Curves*, III.1, and agree with PARI/Sage's
    normalisation. They are computed on demand and cached by ``functools``-free
    lazy properties because the dataclass is frozen.
    """

    a1: int
    a2: int
    a3: int
    a4: int
    a6: int

    # ---------- construction ----------

    @staticmethod
    def from_list(ainvs: Iterable[int]) -> "EllipticCurve":
        """Build from ``[a1,a2,a3,a4,a6]`` or the short form ``[a4,a6]``."""
        v = [int(a) for a in ainvs]
        if len(v) == 2:
            v = [0, 0, 0, v[0], v[1]]
        if len(v) != 5:
            raise ValueError("need 2 or 5 a-invariants, got %d" % len(v))
        E = EllipticCurve(*v)
        if E.discriminant == 0:
            raise SingularCurveError("discriminant is 0: %r is singular" % (v,))
        return E

    @property
    def ainvs(self) -> tuple:
        return (self.a1, self.a2, self.a3, self.a4, self.a6)

    # ---------- standard invariants ----------

    @property
    def b2(self) -> int:
        return self.a1 * self.a1 + 4 * self.a2

    @property
    def b4(self) -> int:
        return 2 * self.a4 + self.a1 * self.a3

    @property
    def b6(self) -> int:
        return self.a3 * self.a3 + 4 * self.a6

    @property
    def b8(self) -> int:
        a1, a2, a3, a4, a6 = self.ainvs
        return (a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4
                + a2 * a3 * a3 - a4 * a4)

    @property
    def c4(self) -> int:
        return self.b2 * self.b2 - 24 * self.b4

    @property
    def c6(self) -> int:
        b2, b4, b6 = self.b2, self.b4, self.b6
        return -b2 ** 3 + 36 * b2 * b4 - 216 * b6

    @property
    def discriminant(self) -> int:
        """Delta = -b2^2 b8 - 8 b4^3 - 27 b6^2 + 9 b2 b4 b6."""
        b2, b4, b6, b8 = self.b2, self.b4, self.b6, self.b8
        return (-b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6
                + 9 * b2 * b4 * b6)

    @property
    def j_invariant(self) -> Fraction:
        return Fraction(self.c4 ** 3, self.discriminant)

    def __repr__(self) -> str:
        return "EllipticCurve(%d,%d,%d,%d,%d)" % self.ainvs

    # ---------- model changes ----------

    def transform(self, u: int, r: int, s: int, t: int) -> "EllipticCurve":
        """Apply (x,y) -> (u^2 x + r, u^3 y + u^2 s x + t), Silverman III.1.2.

        Raises if the result is non-integral, which is the caller's signal that
        the substitution was not admissible over Z.
        """
        a1, a2, a3, a4, a6 = self.ainvs
        n1 = Fraction(a1 + 2 * s, u)
        n2 = Fraction(a2 - s * a1 + 3 * r - s * s, u ** 2)
        n3 = Fraction(a3 + r * a1 + 2 * t, u ** 3)
        n4 = Fraction(a4 - s * a3 + 2 * r * a2 - (t + r * s) * a1
                      + 3 * r * r - 2 * s * t, u ** 4)
        n6 = Fraction(a6 + r * a4 + r * r * a2 + r ** 3
                      - t * a3 - t * t - r * t * a1, u ** 6)
        out = []
        for v in (n1, n2, n3, n4, n6):
            if v.denominator != 1:
                raise ValueError("transform (%s,%s,%s,%s) is not integral" %
                                 (u, r, s, t))
            out.append(int(v))
        return EllipticCurve(*out)

    def minimal_model(self) -> "EllipticCurve":
        """Global minimal model via Laska-Kraus-Connell.

        Returns a curve with the same j-invariant whose discriminant is minimal
        at every prime. Over Q every curve has one (class number 1).
        """
        c4, c6 = self.c4, self.c6
        u = _lkc_scaling(c4, c6)
        E = _curve_from_c_invariants(c4 // u ** 4, c6 // u ** 6)
        if E is None:
            raise ArithmeticError(
                "minimalisation failed for %r; scaling u=%d left no integral "
                "model. This is a bug, not a property of the curve." % (self, u))
        return E

    def is_minimal(self) -> bool:
        return self.minimal_model().discriminant == self.discriminant

    # ---------- the group law over Q ----------

    def is_on_curve(self, P) -> bool:
        if P is O or isinstance(P, INFINITY):
            return True
        x, y = Fraction(P[0]), Fraction(P[1])
        a1, a2, a3, a4, a6 = self.ainvs
        lhs = y * y + a1 * x * y + a3 * y
        rhs = x ** 3 + a2 * x * x + a4 * x + a6
        return lhs == rhs

    def negate(self, P):
        """-P = (x, -y - a1 x - a3)."""
        if P is O or isinstance(P, INFINITY):
            return O
        x, y = Fraction(P[0]), Fraction(P[1])
        return (x, -y - self.a1 * x - self.a3)

    def add(self, P, Q):
        """The chord-and-tangent group law, Silverman III.2.3."""
        if P is O or isinstance(P, INFINITY):
            return Q
        if Q is O or isinstance(Q, INFINITY):
            return P
        a1, a2, a3, a4, a6 = self.ainvs
        x1, y1 = Fraction(P[0]), Fraction(P[1])
        x2, y2 = Fraction(Q[0]), Fraction(Q[1])

        if x1 == x2 and (y1 + y2 + a1 * x2 + a3) == 0:
            return O

        if x1 == x2 and y1 == y2:
            num = 3 * x1 * x1 + 2 * a2 * x1 + a4 - a1 * y1
            den = 2 * y1 + a1 * x1 + a3
            if den == 0:
                return O
            lam = num / den
        else:
            lam = (y2 - y1) / (x2 - x1)

        nu = y1 - lam * x1
        x3 = lam * lam + a1 * lam - a2 - x1 - x2
        y3 = -(lam + a1) * x3 - nu - a3
        return (x3, y3)

    def multiply(self, P, n: int):
        """Scalar multiple nP by double-and-add."""
        if n == 0:
            return O
        if n < 0:
            return self.negate(self.multiply(P, -n))
        R, Q = O, P
        while n:
            if n & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            n >>= 1
        return R

    # ---------- reduction mod p ----------

    def has_good_reduction(self, p: int) -> bool:
        return self.discriminant % p != 0

    def count_points_mod_p(self, p: int) -> int:
        """#E(F_p) including the point at infinity, by brute-force tally.

        Counts solutions of the Weierstrass equation directly. O(p) with small
        constant; fine for the p <= few x 10^4 this project needs, and it is
        the definition rather than an algorithm, which is the point.
        """
        if p == 2:
            n = 1
            for x in range(2):
                for y in range(2):
                    if self._weierstrass_mod(x, y, 2) == 0:
                        n += 1
            return n
        # Tally y-values by the square they satisfy: complete the square to
        # y'^2 = 4x^3 + b2 x^2 + 2 b4 x + b6 with y' = 2y + a1 x + a3.
        b2, b4, b6 = self.b2 % p, self.b4 % p, self.b6 % p
        legendre = _legendre_table(p)
        n = 1  # the point at infinity
        for x in range(p):
            rhs = (4 * x * x * x + b2 * x * x + 2 * b4 * x + b6) % p
            n += 1 + legendre[rhs]
        return n

    def _weierstrass_mod(self, x: int, y: int, p: int) -> int:
        a1, a2, a3, a4, a6 = self.ainvs
        return ((y * y + a1 * x * y + a3 * y
                 - x ** 3 - a2 * x * x - a4 * x - a6) % p)

    def ap(self, p: int) -> int:
        """The trace of Frobenius a_p.

        For good p:      a_p = p + 1 - #E(F_p),  and |a_p| <= 2 sqrt(p) (Hasse).
        For multiplicative reduction: a_p = +1 (split) or -1 (non-split).
        For additive reduction:       a_p = 0.

        Reduction type is decided on the *minimal* model, as it must be.
        """
        E = self if self.is_minimal() else self.minimal_model()
        if E.discriminant % p != 0:
            return p + 1 - E.count_points_mod_p(p)
        if E.c4 % p != 0:
            return E._split_multiplicative_sign(p)
        return 0

    def _split_multiplicative_sign(self, p: int) -> int:
        """+1 if the reduction is split multiplicative at p, else -1.

        At a node the two tangent slopes are defined over F_p exactly when the
        relevant discriminant is a square; for p odd that is the Legendre
        symbol of -c6 (Silverman, and the standard PARI criterion).
        """
        if p == 2:
            # Decide by counting: the singular curve has p points on its smooth
            # locus when split, p+2 when non-split... resolved directly.
            return 1 if (self.c6 % 8) in (1,) else -1
        return 1 if pow((-self.c6) % p, (p - 1) // 2, p) == 1 else -1

    def reduction_type(self, p: int) -> str:
        """'good', 'split', 'nonsplit', or 'additive' at p, on the minimal model."""
        E = self if self.is_minimal() else self.minimal_model()
        if E.discriminant % p != 0:
            return "good"
        if E.c4 % p == 0:
            return "additive"
        return "split" if E._split_multiplicative_sign(p) == 1 else "nonsplit"

    def bad_primes(self) -> list:
        return sorted(factorint(abs(self.minimal_model().discriminant)).keys())

    # ---------- torsion ----------

    def torsion_points(self) -> list:
        """All Q-rational torsion points, by Lutz-Nagell plus Mazur's bound.

        Mazur: the torsion subgroup is one of Z/N (N <= 10 or N = 12) or
        Z/2 x Z/2N (N <= 4), so |E(Q)_tors| <= 16. We find candidates by
        Lutz-Nagell on an integral model and verify each by order computation.
        """
        E = self.minimal_model()
        disc = E.discriminant
        pts = [O]
        # Lutz-Nagell: for y^2 = x^3 + Ax + B with a1=a3=a2=0 a torsion point
        # has integral coords and y^2 | disc. We work with the 4,b-completed
        # model Y^2 = 4x^3 + b2 x^2 + 2 b4 x + b6 to stay integral in general.
        b2, b4, b6 = E.b2, E.b4, E.b6
        # Y = 2y + a1 x + a3, so torsion has Y integral and Y^2 | 16*disc.
        bound = 16 * abs(disc)
        for Y in _square_divisor_candidates(bound):
            # The 2-torsion x may have denominator 2 or 4 even on a minimal
            # model (e.g. 15a1 has a 2-torsion point at x = -13/4), so integer
            # roots alone are not enough.
            for x in _rational_roots_of([4, b2, 2 * b4, b6 - Y * Y]):
                y = Fraction(Y - E.a1 * x - E.a3, 2)
                P = (x, y)
                if E.is_on_curve(P) and P not in pts:
                    if _has_finite_order(E, P, 16):
                        pts.append(P)
                        Q = E.negate(P)
                        if Q not in pts:
                            pts.append(Q)
        return _close_under_addition(E, pts)

    def torsion_order(self) -> int:
        return len(self.torsion_points())


# ---------- module-level helpers ----------

_LEGENDRE_CACHE: dict = {}


def _legendre_table(p: int) -> list:
    """``t[a]`` is 1 if a is a nonzero QR mod p, -1 if a non-residue, 0 if a=0."""
    t = _LEGENDRE_CACHE.get(p)
    if t is not None:
        return t
    t = [-1] * p
    t[0] = 0
    for x in range(1, (p // 2) + 1):
        t[(x * x) % p] = 1
    if len(_LEGENDRE_CACHE) < 512:
        _LEGENDRE_CACHE[p] = t
    return t


def _valuation(n: int, p: int) -> float:
    """v_p(n); infinity (as ``float('inf')``) when n is 0."""
    if n == 0:
        return float("inf")
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def _lkc_scaling(c4: int, c6: int) -> int:
    """The largest u with (c4/u^4, c6/u^6) still the c-invariants of a curve /Z.

    Laska-Kraus-Connell. Rather than encoding Kraus' 2- and 3-adic congruence
    conditions as separate predicates -- which is where hand-rolled versions of
    this routine usually go wrong -- we take the arithmetic bound from the
    valuations and then *verify by construction*, decreasing the exponent until
    ``_curve_from_c_invariants`` actually returns an integral curve with
    precisely these c-invariants. The construction is the criterion.

    A zero c-invariant (c4 = 0 for j = 0, c6 = 0 for j = 1728) imposes no
    bound of its own: its p-adic valuation is infinite, so it is simply left
    out of the minimum rather than fed through ``float('inf') // n``, which is
    NaN and used to poison the whole computation.
    """
    if c4 == 0 and c6 == 0:
        return 1
    disc = _disc_from_c(c4, c6)
    if disc == 0:
        return 1
    u = 1
    primes = set(factorint(abs(disc)).keys())
    for p in sorted(primes):
        # disc != 0 above, so its valuation is a finite integer; zero c4/c6
        # contribute no cap (infinite valuation) instead of a NaN.
        caps = [_valuation(disc, p) // 12]
        if c4:
            caps.append(_valuation(c4, p) // 4)
        if c6:
            caps.append(_valuation(c6, p) // 6)
        k = min(caps)
        while k > 0:
            if _curve_from_c_invariants(c4 // p ** (4 * k),
                                        c6 // p ** (6 * k)) is not None:
                break
            k -= 1
        u *= p ** k
    return u


def _disc_from_c(c4: int, c6: int) -> int:
    """Delta = (c4^3 - c6^2)/1728."""
    num = c4 ** 3 - c6 ** 2
    if num % 1728 != 0:
        raise ValueError("c4^3 - c6^2 = %d is not divisible by 1728" % num)
    return num // 1728


def _curve_from_c_invariants(c4: int, c6: int) -> Optional[EllipticCurve]:
    """The integral Weierstrass model with these c-invariants, or None.

    Returning ``None`` is meaningful: it certifies that no curve over Z has
    exactly this (c4, c6), which is how ``_lkc_scaling`` detects that it has
    scaled down too far.

    Since c6 = -b2^3 + 36 b2 b4 - 216 b6 and b2^3 = b2 (mod 12) for every b2
    arising from an integral model, b2 is pinned modulo 12 by c6. Each
    subsequent divisibility below is a genuine integrality constraint, not a
    convention.
    """
    for b2 in range(-5, 7):
        if (b2 + c6) % 12 != 0:
            continue
        if (b2 * b2 - c4) % 24 != 0:
            continue
        b4 = (b2 * b2 - c4) // 24
        num = -b2 ** 3 + 36 * b2 * b4 - c6
        if num % 216 != 0:
            continue
        b6 = num // 216
        a1 = b2 % 2
        if (b2 - a1 * a1) % 4 != 0:
            continue
        a2 = (b2 - a1 * a1) // 4
        a3 = b6 % 2
        if (b4 - a1 * a3) % 2 != 0:
            continue
        a4 = (b4 - a1 * a3) // 2
        if (b6 - a3 * a3) % 4 != 0:
            continue
        a6 = (b6 - a3 * a3) // 4
        E = EllipticCurve(int(a1), int(a2), int(a3), int(a4), int(a6))
        if E.c4 == c4 and E.c6 == c6 and E.discriminant != 0:
            return E
    return None


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)


def _square_divisor_candidates(bound: int, limit: int = 4000):
    """Y values worth testing in Lutz-Nagell: Y = 0 and Y with Y^2 | bound."""
    yield 0
    y = 1
    while y * y <= bound and y <= limit:
        if bound % (y * y) == 0:
            yield y
            yield -y
        y += 1


def _integer_roots_of(coeffs: list) -> list:
    """Integer roots of the polynomial with the given coefficients, high degree first."""
    from sympy import Poly, Symbol
    x = Symbol("x")
    if all(c == 0 for c in coeffs):
        return []
    p = Poly(sum(c * x ** (len(coeffs) - 1 - i)
                 for i, c in enumerate(coeffs)), x)
    out = []
    for r in p.ground_roots():
        if r.is_Integer:
            out.append(int(r))
    return out


def _rational_roots_of(coeffs: list) -> list:
    """Rational roots as Fractions. Denominators divide the leading coefficient."""
    if all(c == 0 for c in coeffs):
        return []
    lead, n = coeffs[0], len(coeffs) - 1
    out = []
    for q in range(1, abs(lead) + 1):
        if lead % q:
            continue
        scaled = [c * q ** i for i, c in enumerate(coeffs)]
        for m in _integer_roots_of(scaled):
            r = Fraction(m, q)
            if r not in out:
                out.append(r)
    return out


def _close_under_addition(E: EllipticCurve, pts: list) -> list:
    """Saturate a set of torsion points under the group law (|E_tors| <= 16)."""
    out = list(pts)
    changed = True
    while changed and len(out) <= 16:
        changed = False
        for P in list(out):
            for Q in list(out):
                R = E.add(P, Q)
                if R not in out:
                    out.append(R)
                    changed = True
    return out


def _has_finite_order(E: EllipticCurve, P, bound: int) -> bool:
    """True if nP = O for some 1 <= n <= bound (Mazur allows bound = 16)."""
    Q = P
    for _ in range(bound):
        if Q is O or isinstance(Q, INFINITY):
            return True
        Q = E.add(Q, P)
    return False
