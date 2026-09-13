"""Pure-int p-adic arithmetic (no Sage/PARI/sympy p-adics).

A ``Padic`` represents p^val * unit with the unit known modulo p^(K - val)
("absolute precision K": the value is known mod p^K). The zero element has
val = K, meaning v_p >= K, all precision exhausted.

Precision bookkeeping: for a op b the result precision is the *minimum of
the absolute precisions* involved. For add/sub of unequal valuations the
summand with the larger valuation is known less well than the other, and we
charge honestly for the cancellation.
"""

from __future__ import annotations


class Padic:
    """p^val * unit, unit known mod p^(K - val); zero has val = K."""

    __slots__ = ("p", "K", "val", "unit")

    def __init__(self, p: int, K: int, val: int, unit: int):
        if not (val >= 0 and val <= K):
            raise ValueError("need 0 <= val <= K (got val=%d K=%d)" % (val, K))
        self.p = p
        self.K = K
        self.val = val
        m = p ** (K - val)
        u = unit % m
        if val < K and u % p == 0:
            raise ValueError("unit %d is 0 mod %d: not in normal form" % (unit, p))
        self.unit = u

    # -- constructors ------------------------------------------------------

    @staticmethod
    def zero(p: int, K: int) -> "Padic":
        return Padic(p, K, K, 0)

    @staticmethod
    def from_int(p: int, K: int, n: int) -> "Padic":
        """Lift the integer n to precision K, computing its valuation."""
        if n == 0:
            return Padic.zero(p, K)
        val = 0
        m = n if n > 0 else -n
        while m % p == 0:
            m //= p
            val += 1
        if val >= K:
            return Padic.zero(p, K)
        return Padic(p, K, val, (n // p ** val) % p ** (K - val))

    @staticmethod
    def from_val_unit(p: int, K: int, val: int, unit: int) -> "Padic":
        return Padic(p, K, val, unit)

    # -- queries -----------------------------------------------------------

    def is_zero(self) -> bool:
        return self.val >= self.K

    def is_unit(self) -> bool:
        return self.val == 0

    def valuation(self):
        """v_p, or None if the element is zero (v_p >= K)."""
        return None if self.is_zero() else self.val

    def lift(self) -> int:
        """The known representative p^val * unit mod p^K."""
        return self.unit % self.p ** (self.K - self.val) * self.p ** self.val

    # -- arithmetic --------------------------------------------------------

    def _pad_to(self, K2: int) -> "Padic":
        """The same value at (weakly) lower absolute precision K2 <= K."""
        if K2 >= self.K:
            return self
        if self.is_zero() or self.val >= K2:
            return Padic.zero(self.p, K2)
        return Padic(self.p, K2, self.val, self.unit % self.p ** (K2 - self.val))

    def __neg__(self) -> "Padic":
        m = self.p ** (self.K - self.val)
        return Padic(self.p, self.K, self.val, (-self.unit) % m)

    def __add__(self, other: "Padic") -> "Padic":
        self._check(other)
        K = min(self.K, other.K)
        a, b = self._pad_to(K), other._pad_to(K)
        if a.is_zero():
            return b
        if b.is_zero():
            return a
        p = a.p
        if a.val == b.val:
            # Add the units in their common window. If they cancel mod p the
            # valuation rises; renormalising via from_int on the representative
            # is honest because the bracket is known mod p^(K - a.val).
            return Padic.from_int(p, K, (a.unit + b.unit) % p ** (K - a.val)
                                  * p ** a.val)
        # a.val != b.val: the deeper term cannot cancel the unit of the
        # shallower one, so no precision is lost.
        lo, hi = (a, b) if a.val < b.val else (b, a)
        return Padic(p, K, lo.val, lo.unit % p ** (K - lo.val))

    def __sub__(self, other: "Padic") -> "Padic":
        return self.__add__(other.__neg__())

    def _check(self, other: "Padic") -> None:
        if not isinstance(other, Padic) or other.p != self.p:
            raise TypeError("mixing p-adics: %r vs %r" % (self, other))

    def __mul__(self, other: "Padic") -> "Padic":
        self._check(other)
        if self.is_zero() or other.is_zero():
            return Padic.zero(self.p, min(self.K, other.K))
        K = min(self.K, other.K)
        a, b = self._pad_to(K), other._pad_to(K)
        if a.val + b.val >= K:
            return Padic.zero(a.p, K)
        return Padic(a.p, K, a.val + b.val,
                     a.unit * b.unit % a.p ** (K - a.val - b.val))

    def __pow__(self, n: int) -> "Padic":
        if n < 0:
            return Padic.__truediv__(Padic(self.p, self.K, 0, 1), self ** (-n))
        out = Padic(self.p, self.K, 0, 1)
        b, e = self, n
        while e:
            if e & 1:
                out = out * b
            e >>= 1
            if e:
                b = b * b
        return out

    def __truediv__(self, other: "Padic") -> "Padic":
        """Division, with honest precision: dividing by b (valuation vb) costs
        vb absolute digits, since only K - vb digits of 1/b are known.

        A result of negative valuation (va < vb) is refused: this class only
        represents elements of Z_p. Callers wanting Q_p with negative
        valuations should rescale.
        """
        self._check(other)
        if other.is_zero() or other.val >= other.K:
            raise ZeroDivisionError("division by zero / total precision loss")
        K = min(self.K, other.K)
        a, b = self._pad_to(K), other._pad_to(K)
        if a.is_zero():
            return Padic.zero(a.p, K)
        p, va, vb = a.p, a.val, b.val
        if va < vb:
            raise ValueError("quotient has negative valuation %d; not in Z_p"
                             % (va - vb))
        # Result: valuation va - vb, unit digits w_r = K - va; represent at
        # absolute precision K - vb so that val + w = K - vb.
        if va - vb >= K - vb:
            return Padic.zero(p, K - vb)
        inv = pow(b.unit, -1, p ** (K - vb))
        unit = a.unit * inv % p ** (K - va)
        return Padic(p, K - vb, va - vb, unit)

    def __eq__(self, other) -> bool:
        """Equality mod p^K at the common precision."""
        if not isinstance(other, Padic) or other.p != self.p:
            return NotImplemented
        K = min(self.K, other.K)
        a, b = self._pad_to(K), other._pad_to(K)
        if a.is_zero() or b.is_zero():
            return a.is_zero() and b.is_zero()
        if a.val != b.val:
            return False
        return (a.unit - b.unit) % a.p ** (K - a.val) == 0

    def __repr__(self) -> str:
        return "Padic(%d^%d * %d, K=%d)" % (self.p, self.val, self.unit, self.K)


# ---------- module-level p-adic functions ----------

def psqrt(x: Padic) -> Padic:
    """Square root by Hensel lifting on t^2 - unit.

    For odd p, a unit is a square in Z_p iff it is a square mod p; the lift
    doubles precision each Newton step (t <- (t + unit/t)/2), so from a mod-p
    root we reach K digits in O(log K) steps. Non-units and non-residues
    raise: the caller decides how to handle them (e.g. sqrt(p) is not in Z_p).
    """
    if not x.is_unit():
        raise ValueError("psqrt needs a unit (v=%s)" % x.valuation())
    p, K = x.p, x.K
    r = _sqrt_mod_p(x.unit % p, p)
    if r is None:
        raise ValueError("unit is a non-residue mod %d" % p)
    t = r
    prec = 1
    # Hensel: after each step t is correct mod p^prec, then p^2*prec.
    while prec < K:
        m = p ** min(2 * prec, K)
        t = (t + x.unit * pow(t, -1, m)) % m
        t = t * pow(2, -1, m) % m  # Newton for t^2 = unit: (t + u/t)/2
        prec = 2 * prec
    return Padic(p, K, 0, t)


def _sqrt_mod_p(a: int, p: int):
    """Square root of a mod odd prime p, or None (Tonelli-Shanks)."""
    if a % p == 0:
        return None
    if p == 2:
        return a % 2
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
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


def _teichmuller(x: Padic) -> int:
    """w(x) mod p^K: x^(p^(K-1)) equals w * g^(p^(K-1)) and
    v_p(g^(p^n) - 1) >= 1 + n, so this is w to full precision."""
    p, K = x.p, x.K
    return pow(x.unit % p ** K, p ** (K - 1), p ** K)


def _vp_int(n: int, p: int) -> int:
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def plog(x: Padic) -> Padic:
    """The plain p-adic logarithm (NOT divided by log_p(1+p)).

    Teichmuller split x = w(x) * <x>, <x> in 1 + pZ_p (p odd), then
    log <x> = sum_{n>=1} (-1)^(n+1) y^n / n with y = <x> - 1.

    Precision: y is known mod p^K, y has v_p >= 1, and v_p(y^n/n) >=
    n - log_p(n) >= K + 1 once n >= K + 2, so the series is summed to
    n = K + 2 with each term reduced mod p^(K+1); the unknown tail is
    0 mod p^K and the result is honest at absolute precision K.
    """
    if not x.is_unit():
        raise ValueError("plog needs a unit")
    p, K = x.p, x.K
    M = p ** (K + 1)
    w = _teichmuller(x)
    g = x.unit * pow(w, -1, p ** K) % p ** K  # in 1 + p Z_p, mod p^K
    y = (g - 1) % p ** K
    vy = _vp_int(y if y else p, p) or 1  # y = p^vy * u, vy >= 1 for y != 0
    total, ru = 0, 1 % M
    u = y // p ** vy if y else 1
    for n in range(1, K + 3):
        if y == 0:
            break
        # term n = (+/-) y^n / n = p^(n*vy - v_p(n)) * (u^n * (n/p^v)^-1)
        ru = ru * u % M
        v_n = n * vy - _vp_int(n, p)
        tu = ru * pow(n // p ** _vp_int(n, p), -1, M) % M
        piece = tu * p ** v_n % M
        total = (total - piece if n % 2 == 0 else total + piece) % M
    return Padic.from_int(p, K, total % p ** K)


def pexp(z: Padic) -> Padic:
    """exp(z) for v_p(z) >= 1 (p odd): sum z^n / n!.

    Convergence: v_p(z^n/n!) >= n*v_p(z) - (n - s_p(n))/(p - 1) >= n/(p-1)
    in the worst case, so summing to n = (p-1)(K+1) makes every neglected
    term 0 mod p^(K+1). Input with v_p(z) = 0 does not converge: raise.
    """
    if z.is_zero():
        return Padic(z.p, z.K, 0, 1)
    if z.val < 1:
        raise ValueError("pexp needs v_p(z) >= 1 for convergence")
    p, K = z.p, z.K
    M = p ** (K + 1)
    # Track term = p^tv * tu: multiplying by z/n adds z.val - v_p(n) to the
    # valuation, and the unit part is inverted only for n/p^v_p(n), a unit.
    tv, tu, total = 0, 1 % M, 1 % M
    for n in range(1, (p - 1) * (K + 1) + 1):
        tv += z.val - _vp_int(n, p)
        tu = tu * z.unit % M * pow(n // p ** _vp_int(n, p), -1, M) % M
        if tv >= K + 1:
            break  # this and every later term is 0 mod p^(K+1)
        total = (total + (tu * p ** tv % M if tv >= 0 else 0)) % M
    return Padic.from_int(p, K, total % p ** K)
