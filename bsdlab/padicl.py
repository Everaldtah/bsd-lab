"""Mazur-Tate-Teitelbaum p-adic L-series L_p(E,T) (E4 machinery, part 3).

Plus modular symbols [r]^+ in period units (exact Fractions, STEP 1), the
MTT p-adic measure mu on Z_p^* built from the p-stabilised symbol (STEP 2,
Stein-Wuthrich Math. Comp. 82 (2013) sec. 3), and L_p(E,T) = sum_a
mu(a+p^n Z_p)*(1+T)^{log_gamma<a>} by Riemann sums (STEP 3), coefficients
c_0..c_3 with stable-digit precision reported by level comparison.

The first draft (GLM-5.2 worker, quota-killed mid-task) was completed by the
orchestrator.  Three defects were fixed: PNum addition silently dropped the
higher-valuation summand; the measure's second term used [a/p^(n+1)]
instead of [a/p^(n-1)]; and the symbols were in units of Omega^+ = w1 while
p-adic BSD (MTT) interpolates L(E,1)/Omega_E with the Neron real period
Omega_E = c_inf * w1.  Measure values are p-integral (p >= 5), so plain
integers mod p^K are used (the buggy PNum class was deleted).
"""

from __future__ import annotations

from fractions import Fraction

from bsdlab.modsym import (ManinSpace, eigenfunctional, period_scales,
                           raw_twisted_sum, kron)
from bsdlab.padic import Padic, psqrt

__all__ = ["PlusSymbols", "unit_root", "Measure", "lp_series"]


def _vp(n: int, p: int) -> int:
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def unit_root(E, p: int, K: int) -> tuple:
    """(alpha, K2): the unit root of x^2 - a_p x + p in Z_p, as an integer
    mod p^K2, alpha = (a_p + sqrt(a_p^2-4p))/2 with the sqrt branch chosen
    so alpha = a_p mod p (Hensel via padic.psqrt).  K2 = K (psqrt's
    precision, from_int of the discriminant keeps all K digits)."""
    ap = E.ap(p)
    if ap % p == 0:
        raise ValueError("p=%d is not ordinary (a_p = %d)" % (p, ap))
    D = psqrt(Padic.from_int(p, K, ap * ap - 4 * p)).lift()
    if (D - ap) % p != 0:
        D = (-D) % p ** K  # the other square root (was a_p - D: wrong root)
    alpha = (ap + D) * pow(2, -1, p ** K) % p ** K
    if (alpha * alpha - ap * alpha + p) % p ** K or alpha % p == 0:
        raise ArithmeticError("unit root check failed at p=%d" % p)
    return alpha, K


class PlusSymbols:
    """[r]^+ for r = b/m in Q, in units of the Neron period Omega_E.

    [b/m]^+ = s^+ * ( S(1) + lam^+( coords(path_combo(b, m)) ) ) with
    s^+ = period_scales(+1), S(1) = raw_twisted_sum(E,1,N) = lam^+({0,oo}).
    Memoised on b mod m (STEP 1b makes this well-defined).
    """

    def __init__(self, E, N: int):
        self.E, self.N = E, N
        self.M = ManinSpace(N)
        self.lams = eigenfunctional(self.M, E)
        from bsdlab.periods import number_of_real_components
        # Omega^+ = w1 -> Neron period Omega_E = c_inf * w1 (MTT normalisation)
        self.s = (period_scales(E, N, self.M, self.lams)[1]
                  / number_of_real_components(E))
        self.S1 = raw_twisted_sum(E, 1, N, self.M, self.lams)
        self.lam = self.lams[1]
        self.memo = {(0, 1): self._value(0, 1)}

    def _value(self, b: int, m: int) -> Fraction:
        vec = self.M.coords(self.M.path_combo(b, m))
        inner = self.S1 + sum(l * v for l, v in zip(self.lam, vec) if v)
        return self.s * inner

    def value(self, b: int, m: int) -> Fraction:
        """[b/m]^+ (only b mod m matters; reduced before memo lookup)."""
        b %= m
        key = (b, m)
        if key not in self.memo:
            self.memo[key] = self._value(b, m)
        return self.memo[key]


def _fr_mod(q: Fraction, R: int, p: int) -> int:
    if q.denominator % p == 0:
        raise ArithmeticError("modular symbol %s not p-integral at p=%d" % (q, p))
    return q.numerator * pow(q.denominator, -1, R) % R


class Measure:
    """The MTT measure (Stein-Wuthrich 2013, sec. 3; Sage padic_lseries):
        mu(a + p^n Z_p) = alpha^-n [a/p^n]^+  -  alpha^-(n+1) [a/p^(n-1)]^+ ,
    as an integer mod p^K.  The distribution relation is the Hecke relation
    for T_p; scratch/t_padicl.py checks it."""

    def __init__(self, P: PlusSymbols, E, p: int, K: int = 30):
        self.P, self.p, self.K = P, p, K
        self.R = p ** K
        self.alpha, _ = unit_root(E, p, K)
        self.ainv = pow(self.alpha, -1, self.R)
        self.memo = {}

    def mu(self, a: int, n: int) -> int:
        p, R = self.p, self.R
        key = (a % p ** n, n)
        if key not in self.memo:
            s1 = self.P.value(a % p ** n, p ** n)
            s0 = self.P.value(a % p ** (n - 1), p ** (n - 1)) if n > 1                 else self.P.value(0, 1)
            self.memo[key] = (pow(self.ainv, n, R) * _fr_mod(s1, R, p)
                              - pow(self.ainv, n + 1, R) * _fr_mod(s0, R, p)) % R
        return self.memo[key]


def lp_series(meas: Measure, n: int, jmax: int = 4) -> list:
    """Level-n Riemann sum of L_p(E, T), T = gamma^s - 1, gamma = 1 + p:
        L_p^(n)(T) = sum_{j < p^(n-1)} sum_{a0 = 1..p-1}
                     mu(omega(a0) gamma^j + p^n Z_p) (1 + T)^j ,
    sampling each ball at its Teichmuller-times-gamma-power point exactly
    as Sage's padic_lseries.series does, so log_gamma<a> = j is an exact
    integer.  Returns [c_0..c_jmax] mod p^K.  Digits are only meaningful
    where levels n and n+1 agree (the caller compares)."""
    from math import comb
    p, R = meas.p, meas.R
    pn = p ** n
    teich = [pow(a0, p ** (meas.K - 1), R) for a0 in range(1, p)]
    cs = [0] * (jmax + 1)
    gj = 1
    for j in range(p ** (n - 1)):
        tot = sum(meas.mu(t * gj % pn, n) for t in teich) % R
        for i in range(jmax + 1):
            cs[i] = (cs[i] + comb(j, i) * tot) % R
        gj = gj * (1 + p) % pn
    return cs
