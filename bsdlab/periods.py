"""Periods and elliptic logarithms.

Everything here is computed from the arithmetic-geometric mean, following
Cremona, *Algorithms for Modular Elliptic Curves*, section 3.7 (equivalently
Cohen, *A Course in Computational Algebraic Number Theory*, algorithm 7.4.7).
No CAS elliptic machinery is used: the AGM is implemented by hand in mpmath
arithmetic, which is the independence this project is about.

Conventions fixed here and relied upon by ``bsd``:
  * We work on the completed model Y^2 = 4x^3 + b2 x^2 + 2 b4 x + b6 with
    Y = 2y + a1 x + a3, whose roots e1, e2, e3 are the branch points.
  * ``period_lattice`` returns generators (w1, w2) with w1 real > 0 and
    Im(w2/w1) > 0.
  * ``real_period`` returns the BSD / LMFDB Omega: the least positive real
    period TIMES the number of real components -- Omega = w1 when the
    discriminant is negative (one real component), Omega = 2*w1 when it is
    positive (two). A factor of 2 here propagates straight into #Sha, so
    this convention is asserted in the tests, not just documented.
"""

from __future__ import annotations

from fractions import Fraction

from mpmath import mp, mpf, mpc, sqrt as msqrt, pi as mpi, polyroots, quad

from bsdlab.curve import EllipticCurve, INFINITY, O


def _agm(a, b, tol):
    """Complex arithmetic-geometric mean, iterated to within ``tol``.

    For complex arguments the iteration (a,b) -> ((a+b)/2, sqrt(ab)) only
    converges when Re(a * conj(b)) > 0, so that a and b never straddle the
    negative real axis of the sqrt branch; callers pick the signs, we just
    iterate. About 5-6 sweeps reach 70-digit agreement.
    """
    while abs(a - b) > tol * abs(a):
        a, b = (a + b) / 2, msqrt(a * b)
    return a


def _cubic_roots(E: EllipticCurve):
    """Roots e1, e2, e3 of 4x^3 + b2 x^2 + 2 b4 x + b6, ordered by case.

    disc > 0: three real roots, returned descending e1 > e2 > e3.
    disc < 0: the real root e1 first, then e2 = alpha + i beta with beta > 0
    and e3 = conj(e2).  Which of the conjugate pair is "e2" fixes the sign of
    Im(w2); both orderings give the same lattice.
    """
    if E.discriminant == 0:
        raise ValueError("singular curve: the cubic has a repeated root")
    rs = polyroots([mpf(4), mpf(E.b2), mpf(2 * E.b4), mpf(E.b6)],
                   extraprec=mp.dps)
    if E.discriminant > 0:
        return tuple(sorted((r.real for r in rs), reverse=True))
    # The real root is the one with negligible imaginary part; of the
    # conjugate pair, "e2" is the one with positive imaginary part.
    e1 = min(rs, key=lambda r: abs(mpc(r).imag)).real
    e2 = max(rs, key=lambda r: mpc(r).imag)
    e3 = min(rs, key=lambda r: mpc(r).imag)
    return e1, e2, e3


def _half_periods(E: EllipticCurve):
    """The Weierstrass half-periods (w_re, w_im); 2*w_re and 2*w_im generate.

    With a = sqrt(e1-e3), b = sqrt(e1-e2), c = sqrt(e2-e3), the half-periods
    of Y^2 = 4(x-e1)(x-e2)(x-e3) are

        w_re = pi / (2 AGM(a, b))          (real in both cases)
        w_im = pi*i / (2 AGM(a, c))        (imaginary iff disc > 0)

    Derivation for the real case: w_re = K(k)/sqrt(e1-e3) with
    k^2 = (e2-e3)/(e1-e3), and K(k) = pi / (2 AGM(1, k')) with
    k' = sqrt(1-k^2) = sqrt((e1-e2)/(e1-e3)); rescaling the AGM by a
    turns 1, k' into a, b.  The same formulas hold for one real root by
    analytic continuation.  Branches: with principal square roots, a and b
    are real and positive when all roots are real; when e2, e3 are a
    conjugate pair, b = conj(a) by construction (e1 - e2 = conj(e1 - e3)),
    and then AGM(a, b) collapses to the real AGM(Re a, |a|) after one
    sweep.  c = sqrt(e2 - e3) is positive imaginary in the real case and
    has argument pi/4 in the complex case; in both AGM(a, c) converges.
    A naive "flip b if Re(conj(a) b) < 0" rule is wrong: it fires exactly
    when e1 < Re(e2) and wrecks the real period.
    """
    e1, e2, e3 = _cubic_roots(E)
    a = msqrt(e1 - e3)
    b = msqrt(e1 - e2) if E.discriminant > 0 else a.conjugate()
    c = msqrt(e2 - e3)
    tol = mpf(10) ** (-(mp.dps - 5))
    w_re = mpi / (2 * _agm(a, b, tol))
    w_im = mpi * mpc(0, 1) / (2 * _agm(a, c, tol))
    return w_re, w_im


def period_lattice(E: EllipticCurve, prec: int = 60):
    """Generators (w1, w2) of the period lattice, w1 real > 0, Im(w2/w1) > 0.

    w1 is the least positive real period of the lattice; w2 is the second
    generator, purely imaginary when disc > 0 (rectangular lattice).
    """
    with mp.workdps(prec + 10):
        w_re, w_im = _half_periods(E)
        # polyroots returns mpc even for real roots, so scrub the type:
        # w1 is mathematically real, and callers rely on that.
        w1 = mpf((2 * w_re).real)
        w2 = 2 * w_im
        if (w2 / w1).imag < 0:
            w2 = -w2
        with mp.workdps(prec):
            return +mpf(w1), +mpc(w2)


def real_period(E: EllipticCurve, prec: int = 60) -> mpf:
    """The BSD real period Omega = w1 (1 component) or 2*w1 (2 components)."""
    w1, _ = period_lattice(E, prec)
    return w1 * number_of_real_components(E)


def number_of_real_components(E: EllipticCurve) -> int:
    """2 if the discriminant is positive (two real components), else 1."""
    return 2 if E.discriminant > 0 else 1


def elliptic_logarithm(E: EllipticCurve, P, prec: int = 60):
    """z with (wp(z), wp'(z)) = (x(P), Y(P)): the elliptic logarithm of P.

    P is an affine point (Fraction, Fraction) on E; write Y = 2y + a1 x + a3
    for its coordinate on the completed model Y^2 = 4x^3 + b2 x^2 + 2 b4 x + b6.
    Then z = integral of d x / Y from x to +infinity (Cremona 3.7).

    On the unbounded component (x >= e1; all of E(R) when disc < 0) the
    result is a real mpf in [0, w1).  On the bounded component (disc > 0,
    e3 <= x <= e2) it is z = w_im + t with 0 <= t < 2*w_re, returned as an
    mpc; the canonical-height code downstream only ever feeds the real
    return through, but both are provided since both occur (e.g. the
    generator of 37a1 is on the bounded component).

    Sign: the parametrisation z -> (wp(z), wp'(z)) has wp'(z) ~ -2/z^3
    near 0, so a point with Y > 0 (asymptotic to +2x^(3/2)) sits near
    z = w1, i.e. z = 2*w_re - u, while Y < 0 gives z = u.  Equivalently
    Y = wp'(z); the test suite asserts this directly.
    """
    if P is O or isinstance(P, INFINITY):
        return mpf(0)
    with mp.workdps(prec + 10):
        e1, e2, e3 = _cubic_roots(E)
        xf = Fraction(P[0])
        yf = Fraction(P[1])
        x = mpf(xf.numerator) / xf.denominator
        Y = (2 * yf + E.a1 * xf + E.a3)
        Ysign = 1 if Y >= 0 else -1
        w_re, w_im = _half_periods(E)
        f = lambda t: 1 / msqrt(4 * t ** 3 + E.b2 * t * t
                                + 2 * E.b4 * t + E.b6)
        if Y == 0:
            # 2-torsion: z is a half-period, wp(w_re)=e1, wp(w_im)=e3,
            # wp(w_re+w_im)=e2.  (The branch below would integrate into a
            # branch point and return noise.)
            for root, w in ((e1, w_re), (e2, w_re + w_im), (e3, w_im)):
                if abs(x - root) < mpf(10) ** (-(mp.dps - 10)):
                    z = w
                    break
        elif E.discriminant > 0 and x < e2:
            # Bounded component: z = w_im + t, t = int from e3 to x of
            # dx / (2 sqrt((e1-x)(e2-x)(x-e3))), sqrt positive on (e3, e2).
            # Note the orientation is opposite to the unbounded oval: here
            # Y > 0 pairs with t itself, Y < 0 with 2*w_re - t.
            g = lambda t: 1 / (2 * msqrt((e1 - t) * (e2 - t) * (t - e3)))
            t = quad(g, [e3, x])
            if Ysign < 0:
                t = 2 * w_re - t
            z = w_im + t
        else:
            u = quad(f, [x, mp.inf])
            z = u if Ysign < 0 else 2 * w_re - u
        with mp.workdps(prec):
            if isinstance(z, mpc) and abs(z.imag) > mpf(10) ** -prec:
                return +mpc(z)
            return +mpf(z.real if isinstance(z, mpc) else z)
