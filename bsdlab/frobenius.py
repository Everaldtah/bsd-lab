"""p-adic Frobenius on H^1 of an elliptic curve (Kedlaya) and the p-adic E2.

Written by the orchestrator after two GLM-5.2 attempts died on the output
cap (W4p1a).  Pure Python integers; no Sage/PARI.

1. MODEL.  E is replaced by its minimal model, then by the short model
       y^2 = Q(x) = x^3 + A x + B,   A = -27 c4,  B = -54 c6,
   via x_s = 36 x + 3 b2, y_s = 108 (2y + a1 x + a3).  Then
       omega_s := dx_s/(2 y_s) = omega/6,   omega = dx/(2y + a1 x + a3).
   For p >= 5 of good reduction, disc(Q) = 2^4 3^12 Delta... is a p-unit.

2. FROBENIUS (Kedlaya 2001, genus 1; Monsky-Washnitzer).  On H^1_MW(-)
   with basis w_i = x^i dx/y (i = 0, 1) the p-power lift x -> x^p,
   y -> y^p (1 + Er/Q^p)^(1/2), Er := Q(x^p) - Q(x)^p (divisible by p), gives
       Frob(w_i) = p x^(p i + p - 1) sum_k binom(-1/2, k) Er^k
                   / Q^((p-1)/2 + p k)  dx/y .
   Truncation after k = N; Er^k is divisible by p^k.

3. REDUCTION in cohomology.
   Finite poles: with S Q + T Q' = 1 (S, T p-integral since disc(Q) is a
   p-unit), P = P1 Q + P0 and P0 T = q Q + r give
       P / Q^m = (P1 + P0 S + q Q') / Q^(m-1) + r Q' / Q^m,
   and d(r / y^(2m-1)) yields  r Q' dx/y^(2m+1) == 2 r'/(2m-1) dx/y^(2m-1).
   Pole at infinity (level 0): d(x^j y) gives
       x^(j+2) dx/y == -[(2j+1) A x^j + 2 j B x^(j-1)] / (2j+3) dx/y.
   Divisions by p in (2m-1), (2j+3) are absorbed by scaling every numerator
   by p^S up front (S = total p-valuation of all divisors that can occur),
   so all arithmetic is exact in Z/p^M; the result is known mod p^(M-S).

4. E2 (Katz; Mazur-Stein-Tate 2006).  Let U be the unit-root eigenline of
   Frobenius.  With omega_s = w_0/2, eta_s = w_1/2 (same matrix as w_0, w_1),
       eta_s - (E2(E, omega_s)/12) omega_s  in  U.
   If v0 w_0 + v1 w_1 spans U, then E2(E, omega_s) = -12 v0/v1, and by
   weight 2, E2(E, omega) = E2(E, omega_s)/36.
   Gates (scratch/t_frob.py): charpoly(F) = x^2 - a_p x + p, and the Serre
   congruence E2(E, omega) == E_{p+1}(E, omega) mod p with
   E4(E, omega) = c4, E6(E, omega) = -c6 (the Tate-curve normalisation).
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

from bsdlab.curve import EllipticCurve


def _vp(n: int, p: int) -> int:
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def _pmul(a: list, b: list, R: int) -> list:
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return [c % R for c in out]


def _padd(a: list, b: list, R: int) -> list:
    if len(a) < len(b):
        a, b = b, a
    out = list(a)
    for i, y in enumerate(b):
        out[i] = (out[i] + y) % R
    return out


def _divmod_monic(P: list, Q: list, R: int):
    """P = q Q + r with Q monic (list, lowest degree first), deg r < deg Q."""
    P = list(P)
    dq = len(Q) - 1
    if len(P) <= dq:
        return [], P
    q = [0] * (len(P) - dq)
    for i in range(len(P) - 1, dq - 1, -1):
        c = P[i] % R
        if c:
            q[i - dq] = c
            for j in range(dq + 1):
                P[i - dq + j] = (P[i - dq + j] - c * Q[j]) % R
    return q, P[:dq]


def short_model(E: EllipticCurve):
    """(A, B, Emin): y^2 = x^3 + A x + B for the minimal model Emin."""
    Em = E if E.is_minimal() else E.minimal_model()
    return -27 * Em.c4, -54 * Em.c6, Em


def _bezout(A: int, B: int, p: int, R: int):
    """S (deg 1), T (deg 2) mod R with S Q + T Q' = 1, Q = x^3 + A x + B."""
    # unknowns s0 s1 t0 t1 t2; equations for x^0..x^4
    M = [[Fraction(B), 0, Fraction(A), 0, 0],
         [Fraction(A), Fraction(B), 0, Fraction(A), 0],
         [0, Fraction(A), Fraction(3), 0, Fraction(A)],
         [Fraction(1), 0, 0, Fraction(3), 0],
         [0, Fraction(1), 0, 0, Fraction(3)]]
    rhs = [Fraction(1), 0, 0, 0, 0]
    n = 5
    M = [[Fraction(x) for x in row] + [Fraction(r)] for row, r in zip(M, rhs)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        inv = 1 / M[c][c]
        M[c] = [x * inv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    sol = [M[i][n] for i in range(n)]
    for x in sol:
        if x.denominator % p == 0:
            raise ArithmeticError("Bezout not p-integral: bad reduction at p?")
    vals = [x.numerator * pow(x.denominator, -1, R) % R for x in sol]
    return vals[0:2], vals[2:5]


def frobenius_matrix(E: EllipticCurve, p: int, K: int, N: int = None):
    """2x2 integer matrix F mod p^K: F[r][c] = coefficient of w_r in
    Frob(w_c), basis w_i = x^i dx/y on the short model of E.  p >= 5 good."""
    if p < 5:
        raise ValueError("p >= 5 required")
    A, B, Em = short_model(E)
    if Em.discriminant % p == 0:
        raise ValueError("bad reduction at p=%d" % p)
    if N is None:
        N = K + 3
    mmax = (p - 1) // 2 + p * N
    dmax = 3 * mmax + p + 8
    S = (sum(_vp(2 * m - 1, p) for m in range(1, mmax + 1))
         + sum(_vp(2 * j + 3, p) for j in range(dmax + 1)))
    M = K + S + 2
    R = p ** M
    Q = [B % R, A % R, 0, 1]
    Qd = [A % R, 0, 3]
    Qxp = [0] * (3 * p + 1)
    Qxp[0], Qxp[p], Qxp[3 * p] = B % R, A % R, 1
    Qpow = [1]
    for _ in range(p):
        Qpow = _pmul(Qpow, Q, R)
    Er = _padd(Qxp, [(-c) % R for c in Qpow], R)
    Sb, Tb = _bezout(A, B, p, R)
    inv4 = pow(4, -1, R)
    scale = p * p ** S % R
    cols = []
    for i in (0, 1):
        levels = {}
        Ek = [1]
        for k in range(N + 1):
            coef = (-1) ** k * comb(2 * k, k) * pow(inv4, k, R) * scale % R
            num = [0] * (p * i + p - 1) + [c * coef % R for c in Ek]
            m = (p - 1) // 2 + p * k
            levels[m] = _padd(levels.get(m, []), num, R)
            Ek = _pmul(Ek, Er, R)
        spent = 0
        P = []
        for m in range(mmax, 0, -1):
            P = _padd(P, levels.pop(m, []), R)
            if not any(P):
                P = []
                continue
            P1, P0 = _divmod_monic(P, Q, R)
            q2, r = _divmod_monic(_pmul(P0, Tb, R), Q, R)
            Anew = _padd(_padd(P1, _pmul(P0, Sb, R), R), _pmul(q2, Qd, R), R)
            v = _vp(2 * m - 1, p)
            u_inv = pow((2 * m - 1) // p ** v, -1, R)
            rd = [(j * r[j]) % R for j in range(1, len(r))]
            term = []
            for c in rd:
                x = 2 * c * u_inv % R
                if x % p ** v:
                    raise ArithmeticError("scaling budget exhausted (m=%d)" % m)
                term.append(x // p ** v)
            spent += v
            P = _padd(Anew, term, R)
        P = list(P) + [0, 0]
        for d in range(len(P) - 1, 1, -1):
            c = P[d] % R
            if not c:
                continue
            j = d - 2
            v = _vp(2 * j + 3, p)
            u_inv = pow((2 * j + 3) // p ** v, -1, R)
            for idx, mult in ((j, (2 * j + 1) * A), (j - 1, 2 * j * B)):
                if idx < 0 or mult == 0:
                    continue
                x = c * mult * u_inv % R
                if x % p ** v:
                    raise ArithmeticError("scaling budget exhausted (j=%d)" % j)
                P[idx] = (P[idx] - x // p ** v) % R
            P[d] = 0
            spent += v
        # Every stored numerator is (true value) * p^S throughout: dividing a
        # term by p^v consumes divisibility headroom, not the global scale.
        # (First version unscaled by p^(S - spent): off by p^spent.)
        if spent > S:
            raise ArithmeticError("spent more than budget")
        rest = S
        col = []
        for c in P[:2]:
            if c % p ** rest:
                raise ArithmeticError("Frobenius entry not p-integral")
            col.append((c // p ** rest) % p ** K)
        cols.append(col)
    return [[cols[0][0], cols[1][0]], [cols[0][1], cols[1][1]]]


def e2_padic(E: EllipticCurve, p: int, K: int, N: int = None) -> int:
    """E2(E, omega) mod p^K (omega the minimal model's invariant differential),
    as an integer in [0, p^K).  p >= 5 good ORDINARY."""
    Kw = K + 2
    F = frobenius_matrix(E, p, Kw, N)
    R = p ** Kw
    if (F[0][0] + F[1][1]) % p == 0:
        raise ValueError("p=%d is supersingular" % p)
    for v in ((0, 1), (1, 0)):
        for _ in range(Kw + 2):
            v = ((F[0][0] * v[0] + F[0][1] * v[1]) % R,
                 (F[1][0] * v[0] + F[1][1] * v[1]) % R)
        if v[1] % p:
            break
    else:
        raise ArithmeticError("unit-root line has v1 = 0 mod p")
    e2_short = -12 * v[0] * pow(v[1], -1, R) % R
    return e2_short * pow(36, -1, R) % p ** K


def eisenstein_poly(k: int):
    """E_k (normalised, constant term 1) as {(i, j): Fraction} with
    E_k = sum c_ij E4^i E6^j, 4i + 6j = k; solved from q-expansions."""
    from sympy import bernoulli
    mons = [(i, (k - 4 * i) // 6) for i in range(k // 4 + 1)
            if (k - 4 * i) >= 0 and (k - 4 * i) % 6 == 0]
    n = len(mons) + 1

    def ek(w):
        c = Fraction(-2 * w) / Fraction(str(bernoulli(w)))
        return [Fraction(1)] + [c * sum(d ** (w - 1) for d in range(1, m + 1)
                                        if m % d == 0) for m in range(1, n)]

    def mul(a, b):
        return [sum(a[i] * b[m - i] for i in range(m + 1)) for m in range(n)]

    e4, e6 = ek(4), ek(6)
    rows = []
    for i, j in mons:
        s = [Fraction(1)] + [Fraction(0)] * (n - 1)
        for _ in range(i):
            s = mul(s, e4)
        for _ in range(j):
            s = mul(s, e6)
        rows.append(s)
    target = ek(k)
    m = len(mons)
    Mx = [[rows[c][r] for c in range(m)] + [target[r]] for r in range(m)]
    for c in range(m):
        piv = next(r for r in range(c, m) if Mx[r][c] != 0)
        Mx[c], Mx[piv] = Mx[piv], Mx[c]
        inv = 1 / Mx[c][c]
        Mx[c] = [x * inv for x in Mx[c]]
        for r in range(m):
            if r != c and Mx[r][c] != 0:
                f = Mx[r][c]
                Mx[r] = [x - f * y for x, y in zip(Mx[r], Mx[c])]
    coeffs = {mons[c]: Mx[c][m] for c in range(m)}
    for r in range(m, n):     # one extra q-coefficient as a consistency check
        assert sum(coeffs[mons[c]] * rows[c][r] for c in range(m)) == target[r]
    return coeffs


def eisenstein_value_mod_p(E: EllipticCurve, k: int, p: int) -> int:
    """E_k(E, omega) mod p with E4 = c4, E6 = -c6 of the minimal model."""
    Em = E if E.is_minimal() else E.minimal_model()
    e4, e6 = Em.c4, -Em.c6
    tot = 0
    for (i, j), c in eisenstein_poly(k).items():
        if c.denominator % p == 0:
            raise ArithmeticError("E_%d not p-integral at p=%d" % (k, p))
        tot += c.numerator * pow(c.denominator, -1, p) * pow(e4, i, p) * pow(e6, j, p)
    return tot % p
