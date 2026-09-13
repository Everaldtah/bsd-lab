"""Cyclotomic p-adic heights and regulators (Mazur-Stein-Tate), p >= 5 good
ordinary.  W4p1b; written by the orchestrator (GLM-5.2 lane out of quota).

Normalisation = the one that makes MTT p-adic BSD hold in the form
    [T^r] L_p(E, T) = (1 - 1/alpha)^2 * Reg_p * prod c_v * #Sha
                      / (|tors|^2 * log_p(1+p)^r)
(Mazur-Stein-Tate, Doc. Math. Extra Vol. (2006) 577-614; Harvey, "Efficient
computation of p-adic heights", LMS JCM 2008; matches Sage's
padic_height / padic_regulator / Sha.an_padic, which serve as external
anchors in tests_research/t_pheights.py).

p-adic sigma.  On the minimal model with formal parameter t = -x/y,
invariant differential omega = f(t) dt and c = (b2 - E2(E, omega))/12,
sigma is the solution, with sigma = t + O(t^2), of
    (log sigma)' = f * ( -a1/2 - integral (x + c) f dt ).
It lies in Z_p[[t]] exactly because E2 is the p-adic (Katz) value.  Here the
series is computed in exact rational arithmetic with c built from an integer
representative of E2 mod p^M, so coefficient k carries the E2 error times
at most the p-part of k!-type denominators; M is chosen with that margin.

Height.  n = lcm(#E(F_p), lcm_v c_v); Q = nP lies in E^1(Q_p) and in the
identity component at every bad prime.  With x(Q) = a/d^2 (d > 0),
    h_p(P) = -(2/n^2) * log_p( sigma(t(Q)) / d ).
Pairing <P,Q> = (h(P+Q) - h(P) - h(Q))/2, so <P,P> = h(P); Reg_p = det.
Values are returned as (u, e): the p-adic number u / p^e, u mod p^K.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd

from bsdlab.curve import EllipticCurve


def _vp(n: int, p: int) -> int:
    if n == 0:
        return 10 ** 9
    n = abs(n)
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def _smul(a, b, L):
    out = [Fraction(0)] * L
    for i, x in enumerate(a[:L]):
        if x:
            for j in range(min(len(b), L - i)):
                if b[j]:
                    out[i + j] += x * b[j]
    return out


def _sinv(a, L):
    out = [Fraction(0)] * L
    out[0] = 1 / a[0]
    for k in range(1, L):
        s = sum(a[i] * out[k - i] for i in range(1, min(k, len(a) - 1) + 1))
        out[k] = -s / a[0]
    return out


def sigma_series(E: EllipticCurve, e2: int, N: int) -> list:
    """[sigma_0 .. sigma_N] as Fractions (sigma_0 = 0, sigma_1 = 1), E minimal,
    e2 an integer representative of E2(E, omega)."""
    a1, a2, a3, a4, a6 = (Fraction(x) for x in E.ainvs)
    L = N + 6
    # w(t) = t^3 u(t): fixed point of w = t^3 + a1 t w + a2 t^2 w + a3 w^2
    #                                     + a4 t w^2 + a6 w^3
    u = [Fraction(1)] + [Fraction(0)] * (L - 1)
    for _ in range(L):
        u2 = _smul(u, u, L)
        u3 = _smul(u2, u, L)
        new = [Fraction(0)] * L
        new[0] += 1
        for k in range(L):
            if k + 1 < L:
                new[k + 1] += a1 * u[k]
            if k + 2 < L:
                new[k + 2] += a2 * u[k]
            if k + 3 < L:
                new[k + 3] += a3 * u2[k]
            if k + 4 < L:
                new[k + 4] += a4 * u2[k]
            if k + 6 < L:
                new[k + 6] += a6 * u3[k]
        if new == u:
            break
        u = new
    iu = _sinv(u, L)                      # x = t^-2 iu,  y = -t^-3 iu
    # dx/dt = t^-3 * X,  X_k = (k - 2) iu_k ;  2y + a1 x + a3 = t^-3 * Dn
    X = [(k - 2) * iu[k] for k in range(L)]
    Dn = [Fraction(0)] * L
    for k in range(L):
        Dn[k] += -2 * iu[k]
        if k + 1 < L:
            Dn[k + 1] += a1 * iu[k]
    Dn[3] += a3
    f = _smul(X, _sinv(Dn, L), L)          # omega = f dt, f = 1 + a1 t + ...
    assert f[0] == 1
    b2 = a1 * a1 + 4 * a2
    c = (b2 - e2) / 12
    # g = (x + c) f  = t^-2 * G,  G = iu*f + c t^2 f
    G = _smul(iu, f, L)
    for k in range(L - 2):
        G[k + 2] += c * f[k]
    if G[1] != 0:
        raise ArithmeticError("residue of (x + c) omega is nonzero")
    # integral: t^-2 G -> sum_k G_k t^(k-1)/(k-1), k != 1  (offset -1)
    Ai = [Fraction(0)] * L
    for k in range(L):
        if k != 1:
            Ai[k] = G[k] / (k - 1)
    # A = (-a1/2 - Ai t^-1) f  = t^-1 * ( -a1/2 t - Ai ) f
    inner = [-x for x in Ai]
    inner[1] += -a1 / 2
    A = _smul(inner, f, L)                 # offset -1
    if A[0] != 1:
        raise ArithmeticError("leading 1/t coefficient %s != 1" % A[0])
    B = A[1:]                              # (log theta)' = B, theta = sigma/t
    if B[0] != a1 / 2:
        raise ArithmeticError("constant term %s != a1/2" % B[0])
    theta = [Fraction(1)] + [Fraction(0)] * N
    for k in range(N):
        theta[k + 1] = sum(B[i] * theta[k - i] for i in range(k + 1)) / (k + 1)
    return [Fraction(0)] + theta[:N]


def _plog_unit(u: int, p: int, M: int) -> int:
    """log_p(u) mod p^M for an integer p-unit u (log of Teichmuller part 0)."""
    R = p ** (M + 3)
    w = pow(u, p - 1, R)
    y = (w - 1) % R
    total = 0
    yk = 1
    k = 1
    while True:
        yk = yk * y % R
        if k - _vp(k, p) >= M + 3:
            break
        term = yk // p ** _vp(k, p) * pow(k // p ** _vp(k, p), -1, R)
        total = (total + (term if k % 2 else -term)) % R
        k += 1
    return total * pow(p - 1, -1, R) % p ** M


def _tamagawa_multiple(E: EllipticCurve) -> int:
    from bsdlab.reduction import tamagawa_numbers
    m = 1
    for c in tamagawa_numbers(E).values():
        m = m * c // gcd(m, c)
    return m


class PadicHeight:
    """h_p on E(Q) for the minimal model E, p >= 5 good ordinary."""

    def __init__(self, E: EllipticCurve, p: int, K: int = 10):
        from bsdlab.frobenius import e2_padic
        if not E.is_minimal():
            raise ValueError("pass the minimal model")
        if E.discriminant % p == 0 or E.ap(p) % p == 0 or p < 5:
            raise ValueError("need p >= 5 good ordinary")
        self.E, self.p, self.K = E, p, K
        n1 = p + 1 - E.ap(p)
        n2 = _tamagawa_multiple(E)
        self.n = n1 * n2 // gcd(n1, n2)
        self.e = 2 * _vp(self.n, p)
        self.N = K + self.e + 8                  # sigma terms
        self.M = self.N + 2 * (self.N // (p - 1)) + 8   # E2 digits
        self.sigma = sigma_series(E, e2_padic(E, p, self.M), self.N)
        self.S = max(_vp(s.denominator, p) for s in self.sigma[1:])

    def __call__(self, P) -> tuple:
        """(u, e) with h_p(P) = u / p^e, u known mod p^K."""
        p, K = self.p, self.K
        Q = self.E.multiply(P, self.n)
        if Q is None or Q == "O" or not isinstance(Q, tuple):
            return (0, self.e)                   # torsion
        x, y = Q
        d2 = x.denominator
        d = int(round(d2 ** 0.5)) if d2 < 2 ** 52 else _isqrt(d2)
        if d * d != d2 or d % p:
            raise ArithmeticError("nP not in the formal group at p")
        t = -x / y
        work = K + self.e + self.S + 6
        R = p ** work
        tot = Fraction(0)
        tk = Fraction(1)
        for k in range(1, self.N + 1):
            tk *= t
            tot += self.sigma[k] * tk
        val = tot / d * p ** self.S
        if val.denominator % p:
            num = val.numerator * pow(val.denominator, -1, R) % R
        else:
            raise ArithmeticError("unexpected p in denominator")
        if num % p ** self.S:
            raise ArithmeticError("sigma(t)/d is not a unit")
        unit = num // p ** self.S
        lg = _plog_unit(unit, p, work - self.S)
        nn = self.n * self.n
        nu = nn // p ** _vp(nn, p)
        h = (-2 * lg * pow(nu, -1, p ** (work - self.S))) % p ** K
        return (h, self.e)


def _isqrt(n: int) -> int:
    from math import isqrt
    return isqrt(n)


def padic_regulator(E: EllipticCurve, gens: list, p: int, K: int = 10):
    """(u, e): Reg_p = u / p^e mod p^K for a Z-basis gens of E(Q)/tors."""
    H = PadicHeight(E, p, K + 4)
    r = len(gens)
    R = p ** (K + 4)
    Mx = [[0] * r for _ in range(r)]
    hs = [H(P)[0] for P in gens]
    inv2 = pow(2, -1, R)
    for i in range(r):
        Mx[i][i] = hs[i]
        for j in range(i + 1, r):
            hij = H(E.add(gens[i], gens[j]))[0]
            Mx[i][j] = Mx[j][i] = (hij - hs[i] - hs[j]) * inv2 % R
    det = _det(Mx, R)
    return det % p ** K, r * H.e


def _det(M: list, R: int) -> int:
    r = len(M)
    if r == 1:
        return M[0][0] % R
    if r == 2:
        return (M[0][0] * M[1][1] - M[0][1] * M[1][0]) % R
    tot = 0
    for j in range(r):
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        tot += (-1) ** j * M[0][j] * _det(minor, R)
    return tot % R
