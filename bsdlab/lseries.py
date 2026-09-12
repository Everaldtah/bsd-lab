"""The Hasse-Weil L-function L(E,s) of an elliptic curve over Q.

L(E,s) = sum_{n>=1} a_n n^{-s}, where the a_n are the Fourier coefficients
of the newform attached to E, computed from a_p data by multiplicativity.

The completed L-function is

    Lambda(s) = N^{s/2} (2pi)^{-s} Gamma(s) L(E,s),   Lambda(s) = w Lambda(2-s),

with w the root number.  Writing alpha = sqrt(N)/(2pi) and x_n = 2pi n/sqrt(N)
(so x_n = n/alpha), the exact rapidly-convergent identity used throughout is

    Lambda(s) = sum_n a_n [Phi_n(s) + w Phi_n(2-s)],
    Phi_n(s) = (alpha/n)^s Gamma(s, x_n),

with Gamma(s,x) the upper incomplete gamma function.  Differentiating r
times at the central point s=1 (where the two terms collide) gives

    Lambda^(r)(1)/r! = (1 + w(-1)^r) * sum_n a_n G_r(x_n),
    G_r(x) = (1/r!) int_1^infty (log u)^r exp(-x u) du,

the factor (1 + w(-1)^r) expressing that Lambda^(r)(1) vanishes unless r has
the parity of w (w=+1 -> r even, w=-1 -> r odd).  Finally L^(r)(1)/r! is
read off from the Taylor coefficients of L = Lambda/g, g(s) = alpha^s
Gamma(s), by the convolution L_r = (sum_j Lambda_{r-j} h_j)/alpha where
h(s) = g(1)/g(s) = alpha^{1-s} Gamma(1)/Gamma(s); log h has coefficients
d_1 = gamma_E - log(alpha) (gamma_E the Euler-Mascheroni constant) and
d_k = (-1)^{k+1} zeta(k)/k for k >= 2, exponentiated by the standard
recurrence.  For r = 0 the whole chain collapses to the familiar

    L(E,1) = (1+w) * sum_n (a_n/n) exp(-2 pi n / sqrt(N)).

Number of terms: the smoothed summand carries exp(-x_n), so truncating at
n_terms = TERMS_PER_SQRT_N*sqrt(N) + TERMS_CONSTANT = 10 sqrt(N) + 100 puts
the first omitted term near exp(-20 pi) ~ 4e-28, below the working
precision for prec <= 50 or so; raise n_terms for higher precision.

The root number is computed as a product of local root numbers over the bad
primes, w = -prod w_p (the leading -1 is the archimedean place), following
the convention: split multiplicative w_p = -1, non-split w_p = +1.  Every
value so computed is cross-checked against the sign read off numerically
from the functional equation itself via the Fricke involution: for
f(iy) = sum a_n exp(-2pi n y) one has exactly f(-1/(Nz)) = eps N z^2 f(z)
with w = -eps, so at y = 2/sqrt(N),

    w = f(i/(2 sqrt(N))) / (4 f(2i/sqrt(N))).

Additive primes are NOT implemented locally and raise NotImplementedError;
callers fall back to the numeric sign above.
"""

from __future__ import annotations

import math
from typing import Optional

from mpmath import mp, mpf
from sympy import factorint

from bsdlab.curve import EllipticCurve

try:
    from bsdlab.reduction import conductor as _conductor
except ImportError:  # bsdlab.reduction is written in parallel; optional
    _conductor = None

DEFAULT_PREC = 60
#: terms in smoothed sums: n_terms = TERMS_PER_SQRT_N*sqrt(N) + TERMS_CONSTANT
TERMS_PER_SQRT_N = 10
TERMS_CONSTANT = 100
#: |L^(r)(1)/r!| below 10^-(prec/2) counts as zero to working precision.
ZERO_THRESHOLD = mpf(10) ** (-DEFAULT_PREC // 2)


def _zero_threshold(prec: int) -> mpf:
    """10^-(prec/2): the ZERO_THRESHOLD rule scaled to `prec` digits."""
    return mpf(10) ** (-(prec // 2))


def _resolve_N(E: EllipticCurve, N: Optional[int]) -> int:
    """The conductor, taken from the explicit argument when given."""
    if N is not None:
        return int(N)
    if _conductor is None:
        raise ImportError(
            "N was not given and bsdlab.reduction is unavailable; pass the "
            "conductor N explicitly")
    return int(_conductor(E))


def _default_n_terms(N: int) -> int:
    """TERMS_PER_SQRT_N*sqrt(N) + TERMS_CONSTANT, as documented above."""
    return TERMS_PER_SQRT_N * math.isqrt(N) + TERMS_CONSTANT


def an_coefficients(E: EllipticCurve, n_max: int, N: Optional[int] = None) -> list:
    """a_1..a_{n_max} by sieve: 1-indexed, index 0 unused (0).

    For primes p: good reduction gives a_p = p + 1 - #E(F_p) (via
    ``E.ap``), split multiplicative a_p = 1, non-split a_p = -1, additive
    a_p = 0.  Prime powers by the recurrences

        good p: a_{p^{k+1}} = a_p a_{p^k} - p a_{p^{k-1}}
        bad  p: a_{p^k} = (a_p)^k

    and a_{mn} = a_m a_n for coprime m, n, since the local Euler factors
    (1 - a_p p^{-s} + p^{1-2s})^{-1} (good) or (1 - a_p p^{-s})^{-1} (bad)
    determine every coefficient.  One smallest-prime-factor sieve fills
    everything; no n is factored on its own.
    """
    if n_max < 1:
        raise ValueError("n_max must be >= 1")
    a = [0] * (n_max + 1)
    a[1] = 1
    # Smallest-prime-factor table; spf[p] == p marks p prime.
    spf = list(range(n_max + 1))
    for i in range(2, math.isqrt(n_max) + 1):
        if spf[i] == i:
            for j in range(i * i, n_max + 1, i):
                if spf[j] == j:
                    spf[j] = i
    bad = set(E.bad_primes())
    for p in range(2, n_max + 1):
        if spf[p] != p:
            continue
        if p not in bad:
            ap = E.ap(p)
            # a_{p^{k+1}} = a_p a_{p^k} - p a_{p^{k-1}}, starting a_{p^0}=1.
            cur, prev = ap, 1
            pk = p
            while pk <= n_max:
                a[pk] = cur
                cur, prev = ap * cur - p * prev, cur
                pk *= p
        else:
            ap = {"split": 1, "nonsplit": -1}.get(E.reduction_type(p), 0)
            pk = p
            while pk <= n_max:
                a[pk] = a[pk // p] * ap  # a_{p^k} = (a_p)^k
                pk *= p
    for n in range(2, n_max + 1):
        if spf[n] == n:
            continue  # primes and their powers were filled above
        p = spf[n]
        m = n
        while m % p == 0:
            m //= p
        # m == 1 means n is a prime power (already set); otherwise split.
        if m > 1:
            a[n] = a[n // m] * a[m]
    return a


def _local_root_number(E: EllipticCurve, p: int, vN: int) -> int:
    """w_p at a bad prime p from the reduction type.

    split multiplicative  -> -1 ; non-split multiplicative -> +1 (exact,
    classical).  Additive reduction needs the Halberstadt/Kraus local
    tables: they are not implemented here, by design rather than by guess.
    """
    kind = E.reduction_type(p)
    if kind == "split":
        return -1
    if kind == "nonsplit":
        return 1
    if p >= 5:
        raise NotImplementedError(
            "local root number at additive p=%d (potentially good, "
            "v_p(N)=%d): the p>=5 additive table (Halberstadt) is not "
            "implemented; use root_number_numeric" % (p, vN))
    raise NotImplementedError(
        "local root number at additive p=%d (v_p(N)=%d): additive "
        "reduction at 2 or 3 is not implemented; use root_number_numeric"
        % (p, vN))


def root_number(E: EllipticCurve, N: Optional[int] = None) -> int:
    """The sign w of the functional equation Lambda(s) = w Lambda(2-s).

    Computed as w = -prod_{p | N} w_p (the leading -1 is the archimedean
    place), then cross-checked against the numerical sign implied by the
    functional equation (``root_number_numeric``); a mismatch raises.
    """
    N = _resolve_N(E, N)
    w = -1
    for p, vN in factorint(N).items():
        w *= _local_root_number(E, p, vN)
    if root_number_numeric(E, N) != w:
        raise ArithmeticError(
            "local root number product %d disagrees with the functional-"
            "equation sign for %r (N=%d)" % (w, E, N))
    return w


def root_number_numeric(E: EllipticCurve, N: Optional[int] = None,
                        prec: int = DEFAULT_PREC) -> int:
    """w read off the functional equation via the Fricke involution.

    For f(iy) = sum a_n exp(-2 pi n y) the exact Fricke relation is
    f(i/(Ny)) = w N y^2 f(iy) (see module docstring).  Evaluating both
    sides at y = 2/sqrt(N) -- off the self-dual point y = 1/sqrt(N), with
    both q-series decaying like exp(-pi n/sqrt(N)) -- gives
    w = f(i/(2 sqrt(N))) / (4 f(2i/sqrt(N))).  The result is +1 or -1; we
    round and insist the ratio sits within 0.01 of an integer sign.
    """
    N = _resolve_N(E, N)
    with mp.workdps(prec):
        n_max = int(12 * mp.sqrt(N)) + 50
        a = an_coefficients(E, n_max, N)
        sN = mp.sqrt(N)
        hi = sum(a[n] * mp.exp(-4 * mp.pi * n / sN) for n in range(1, n_max + 1))
        lo = sum(a[n] * mp.exp(-mp.pi * n / sN) for n in range(1, n_max + 1))
        ratio = lo / (4 * hi)
        w = 1 if ratio > 0 else -1
        if abs(ratio - w) > mpf("0.01"):
            raise ArithmeticError(
                "Fricke ratio %s is not close to +-1 for %r (N=%d)"
                % (mp.nstr(ratio, 6), E, N))
        return w


def sign(E: EllipticCurve, N: Optional[int], prec: int) -> int:
    """root_number where computable, else the numeric Fricke sign."""
    try:
        return root_number(E, N)
    except NotImplementedError:
        return root_number_numeric(E, N, prec)


def _G(r: int, x) -> mpf:
    """G_r(x) = (1/r!) int_1^infty (log u)^r e^{-x u} du.

    This is Phi_n^(r)(1)/r! for x = x_n, the building block of the
    central-point derivatives of Lambda (module docstring).  For r = 0 it
    collapses to e^{-x}/x.  The integrand decays exponentially, so a
    plain quadrature over [1, inf) is fast and stable.
    """
    if r == 0:
        return mp.exp(-x) / x
    c = mpf(1) / mp.factorial(r)

    def integrand(u):
        return mp.power(mp.log(u), r) * mp.exp(-x * u) * c
    return mp.quad(integrand, [1, mp.inf])


def _h_coeffs(alpha, rmax: int) -> list:
    """Taylor coefficients h_0..h_rmax of h(s) = g(1)/g(s) at s=1.

    log h(s) = (1-s) log(alpha) - log Gamma(s) + log Gamma(1), whose
    Taylor coefficients at s=1 are d_1 = euler - log(alpha) (from
    psi(1) = -euler) and d_k = (-1)^{k+1} zeta(k)/k for k >= 2 (from
    psi^{(k-1)}(1) = (-1)^k (k-1)! zeta(k)).  The h_j follow from the
    standard exp-series recurrence h_j = (1/j) sum_{k<=j} k d_k h_{j-k}.
    """
    log_alpha = mp.log(alpha)
    d = [None, mp.euler - log_alpha]
    for k in range(2, rmax + 1):
        d.append((-1) ** (k + 1) * mp.zeta(k) / k)
    h = [mpf(1)]
    for j in range(1, rmax + 1):
        h.append(sum(k * d[k] * h[j - k] for k in range(1, j + 1)) / j)
    return h


def _completed_moments(E: EllipticCurve, rmax: int, N: Optional[int] = None,
                       prec: int = DEFAULT_PREC) -> tuple:
    """(w, alpha, [L_0..L_rmax]) with L_r = L^(r)(1)/r!.

    Chain (module docstring): Lambda_r = (1 + w(-1)^r) sum_n a_n G_r(x_n),
    then L_r = (sum_j Lambda_{r-j} h_j) / alpha.
    """
    N = _resolve_N(E, N)
    with mp.workdps(prec):
        w = sign(E, N, prec)
        alpha = mp.sqrt(N) / (2 * mp.pi)
        n_terms = _default_n_terms(N)
        a = an_coefficients(E, n_terms, N)
        xs = [2 * mp.pi * n / mp.sqrt(N) for n in range(1, n_terms + 1)]
        lam = []
        for r in range(rmax + 1):
            s = mp.mpf(0)
            for n in range(1, n_terms + 1):
                if a[n]:
                    s += a[n] * _G(r, xs[n - 1])
            lam.append((1 + w * (-1) ** r) * s)
        h = _h_coeffs(alpha, rmax)
        L = [sum(lam[r - j] * h[j] for j in range(r + 1)) / alpha
             for r in range(rmax + 1)]
        return w, alpha, L


def lambda_value(E: EllipticCurve, s, n_terms: Optional[int] = None,
                 N: Optional[int] = None, prec: int = DEFAULT_PREC) -> mpf:
    """Lambda(E,s) = N^{s/2}(2pi)^{-s}Gamma(s)L(E,s) by the exact afe.

    Uses the identity Lambda(s) = sum_n a_n[Phi_n(s) + w Phi_n(2-s)] with
    Phi_n(s) = (alpha/n)^s Gamma(s, x_n) (module docstring), exact up to
    the truncation at n_terms (default 10 sqrt(N) + 100 terms; the tail is
    damped by exp(-2 pi n_terms/sqrt(N)) ~ exp(-20 pi)).
    """
    N = _resolve_N(E, N)
    with mp.workdps(prec):
        s = mp.mpf(s) if not isinstance(s, complex) else mp.mpc(s)
        w = sign(E, N, prec)
        alpha = mp.sqrt(N) / (2 * mp.pi)
        if n_terms is None:
            n_terms = _default_n_terms(N)
        a = an_coefficients(E, n_terms, N)
        total = mp.mpf(0)
        for n in range(1, n_terms + 1):
            if not a[n]:
                continue
            x = 2 * mp.pi * n / mp.sqrt(N)
            phi = (alpha / n) ** s * mp.gammainc(s, x, mp.inf)
            phi_tilde = (alpha / n) ** (2 - s) * mp.gammainc(2 - s, x, mp.inf)
            total += a[n] * (phi + w * phi_tilde)
        return +total


def l_value(E: EllipticCurve, N: Optional[int] = None,
            prec: int = DEFAULT_PREC) -> mpf:
    """L(E,1) = (1+w) sum_n (a_n/n) exp(-2 pi n/sqrt(N)).

    For w = -1 the functional equation forces L(E,1) = 0 and the formula
    returns exactly that; for w = +1 it is the r=0 case of the general
    central-point machinery (see ``_completed_moments``).
    """
    return _completed_moments(E, 0, N=N, prec=prec)[2][0]


def l_derivative(E: EllipticCurve, order: int, N: Optional[int] = None,
                 prec: int = DEFAULT_PREC) -> mpf:
    """L^(order)(E,1)/order!, via the G_r-smoothed approximate FE.

    L^(r)(1)/r! = (sum_j Lambda_{r-j} h_j)/alpha with
    Lambda_r = (1 + w(-1)^r) sum_n a_n G_r(2 pi n/sqrt(N)); see the module
    docstring for the derivation of G_r and h.
    """
    if order < 0:
        raise ValueError("order must be >= 0")
    return _completed_moments(E, order, N=N, prec=prec)[2][order]


def analytic_rank(E: EllipticCurve, N: Optional[int] = None,
                  prec: int = DEFAULT_PREC, max_rank: int = 4) -> int:
    """The analytic rank: least r with L^(r)(1)/r! != 0.

    The functional equation forces Lambda^(r)(1) = 0 whenever
    (-1)^r != w, and that vanishing is built into the evaluation (the
    factor (1 + w(-1)^r)), so only r of the parity implied by the root
    number are ever tested: even r for w = +1, odd r for w = -1.  A
    wrong-parity value can therefore never be reported nonzero; the parity
    cross-check is structural.  "Nonzero" means exceeding the threshold
    10^-(prec/2) (ZERO_THRESHOLD at default precision).
    """
    w, _, L = _completed_moments(E, max_rank, N=N, prec=prec)
    threshold = _zero_threshold(prec)
    for r in range(0 if w == 1 else 1, max_rank + 1, 2):
        if abs(L[r]) > threshold:
            if r % 2 != (0 if w == 1 else 1):
                raise ArithmeticError(
                    "first nonzero derivative r=%d has the wrong parity "
                    "for w=%d on %r" % (r, w, E))
            return r
    raise ArithmeticError(
        "no nonzero L^(r)(1)/r! for r <= %d at prec=%d for %r; raise "
        "max_rank or prec" % (max_rank, prec, E))


def leading_coefficient(E: EllipticCurve, N: Optional[int] = None,
                        prec: int = DEFAULT_PREC) -> tuple:
    """(r, L^(r)(1)/r!) at the analytic rank r: the BSD leading term."""
    r = analytic_rank(E, N=N, prec=prec)
    return r, l_derivative(E, r, N=N, prec=prec)
