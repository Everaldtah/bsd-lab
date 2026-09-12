"""The Mordell-Weil group: point search, canonical heights, regulator.

Canonical heights are computed as a SUM OF LOCAL HEIGHTS following Silverman,
"Computing canonical heights with little (or no) factorization", Math. Comp. 51
(1988) -- Section 4 for the archimedean place (Tate's rapidly convergent
theta-series, with the translation trick when |x| < 1/2) and Section 5 for the
finite places.  This is the same algorithm Sage uses, and it converges to full
working precision in ~100 terms, so no 4^-n naive limit is ever taken.

Normalisation: hhat(P) is the height with hhat(nP) = n^2 hhat(P), calibrated
against 37a1 where hhat((0,0)) = 0.0511114082... (LMFDB).  The regulator is
det(<P_i,P_j>) for the height pairing <P,Q> = (hhat(P+Q)-hhat(P)-hhat(Q))/2.

All point arithmetic is exact (Fractions); mpmath enters only when logarithms
must be taken.
"""

from __future__ import annotations

from fractions import Fraction
from math import isqrt

from mpmath import mp, mpf

from bsdlab.curve import EllipticCurve, INFINITY, O

# Tolerance used by independent_points: a pivot (or determinant) of the height
# pairing matrix whose magnitude is below this fraction of the largest matrix
# entry is treated as zero, i.e. the corresponding point as dependent.  Heights
# of searched points are O(1), and pairing entries are computed to ~60 digits,
# so 1e-25 leaves an enormous margin in both directions.
INDEPENDENCE_TOLERANCE = mpf("1e-25")


# ---------- naive height ----------

def naive_height(P) -> mpf:
    """h(P) = log max(|num(x)|, |den(x)|) for x(P) in lowest terms; h(O) = 0."""
    if P is O or isinstance(P, INFINITY):
        return mpf(0)
    x = Fraction(P[0])
    return mp.log(max(abs(x.numerator), x.denominator))


# ---------- local heights (Silverman 1988) ----------

def _to_minimal(E: EllipticCurve, P):
    """Map P onto E's global minimal model, exactly, via (u, r, s, t).

    curve.py's minimal_model() does not record the transformation it applied,
    so we recover it from the invariants: u from Delta/Delta_min (a 12th power),
    then r, s, t from b2, a1, a3.  The result is verified to lie on the minimal
    model -- the construction is the criterion, as in curve._lkc_scaling.
    """
    Emin = E.minimal_model()
    if Emin.ainvs == E.ainvs:
        return Emin, P
    from sympy import integer_nthroot
    u = integer_nthroot(abs(E.discriminant) // abs(Emin.discriminant), 12)[0]
    r = (u * u * Emin.b2 - E.b2) // 12
    s = (u * Emin.a1 - E.a1) // 2
    t = (u ** 3 * Emin.a3 - E.a3 - r * E.a1) // 2
    if P is O or isinstance(P, INFINITY):
        return Emin, O
    x, y = Fraction(P[0]), Fraction(P[1])
    Q = (u * u * x + r, u ** 3 * y + u * u * s * x + t)
    if not Emin.is_on_curve(Q):
        raise ArithmeticError("point %r does not map onto the minimal model "
                              "of %r" % (P, E))
    return Emin, Q


def _archimedean_local_height(E: EllipticCurve, P, prec: int) -> mpf:
    """lambda_infinity(P), Silverman (1988) section 4, on the given model.

    Tate's series for the sigma function, evaluated in the local parameter
    t = 1/x when |x| >= 1/2 and t = 1/(x+1) on the r = -1 translate otherwise;
    the series may hop between the two models, which is what `beta` tracks.
    Each term gains a factor 4, so ~0.51*dps bits of precision need a like
    number of terms.
    """
    if P is O or isinstance(P, INFINITY):
        return mpf(0)
    x = Fraction(P[0])
    b2, b4, b6, b8 = E.b2, E.b4, E.b6, E.b8
    with mp.workdps(prec + 15):
        xm = mpf(x.numerator) / mpf(x.denominator)
        # translate the model by r = -1 when |x| is small
        b2p, b4p = b2 - 12, b4 - b2 + 6
        b6p, b8p = b6 - 2 * b4 + b2 - 4, b8 - 3 * b6 + 3 * b4 - b2 + 3

        def fz(T, b4_, b6_, b8_):
            return 1 - T * T * (b4_ + T * (2 * b6_ + T * b8_))

        def fw(T, b2_, b4_, b6_):
            return T * (4 + T * (b2_ + T * (2 * b4_ + T * b6_)))

        beta = abs(xm) >= mpf("0.5")
        t = 1 / xm if beta else 1 / (xm + 1)
        lam = -mp.log(abs(t))
        mu = mpf(0)
        four_to_n = mpf(1)
        H = max(4, abs(b2), 2 * abs(b4), 2 * abs(b6), abs(b8))
        nterms = int(mp.ceil(mp.mpf("0.51") * mp.prec + mpf("0.5")
                              + mpf("0.75") * mp.log(7 + 4 * mp.log(H) / 3)))
        for _ in range(nterms):
            if beta:
                w = fw(t, b2, b4, b6)
                z = fz(t, b4, b6, b8)
                if abs(w) <= 2 * abs(z):
                    mu += four_to_n * mp.log(abs(z))
                    t = w / z
                else:
                    mu += four_to_n * mp.log(abs(z + w))
                    t = w / (z + w)
                    beta = not beta
            else:
                w = fw(t, b2p, b4p, b6p)
                z = fz(t, b4p, b6p, b8p)
                if abs(w) <= 2 * abs(z):
                    mu += four_to_n * mp.log(abs(z))
                    t = w / z
                else:
                    mu += four_to_n * mp.log(abs(z - w))
                    t = w / (z - w)
                    beta = not beta
            four_to_n /= 4
        return mp.mpf(lam + mu / 4)


def _vp(q: Fraction, p: int) -> float:
    """v_p of a nonzero rational; +infinity (as a float) for q = 0."""
    if q == 0:
        return float("inf")
    n, d = q.numerator, q.denominator
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    while d % p == 0:
        d //= p
        v -= 1
    return v


def _non_archimedean_local_height(E: EllipticCurve, P, prec: int) -> mpf:
    """Sum of the finite-place local heights, on a minimal model.

    On a minimal model the total finite contribution is
        log den(x)  +  sum over bad p not dividing den(x) of r_p log p,
    where r_p is Silverman (1988) section 5's valuation expression.  (Bad
    primes dividing den(x) contribute r_p = -v_p(den(x)) there, cancelling
    exactly against the log den(x) term, so they may simply be skipped.)
    """
    if P is O or isinstance(P, INFINITY):
        return mpf(0)
    a1, a2, a3, a4, _ = E.ainvs
    b2, b4, b6, b8 = E.b2, E.b4, E.b6, E.b8
    x, y = Fraction(P[0]), Fraction(P[1])
    den = x.denominator
    total = mpf(0)
    for p in E.bad_primes():
        if den % p == 0:
            continue
        N = _vp(Fraction(E.discriminant), p)
        A = _vp(3 * x * x + 2 * a2 * x + a4 - a1 * y, p)
        B = _vp(2 * y + a1 * x + a3, p)
        if A <= 0 or B <= 0:
            vx = _vp(x, p)
            r = Fraction(-int(vx)) if vx < 0 else Fraction(0)
        elif _vp(Fraction(E.c4), p) == 0:
            n = min(int(B), Fraction(N, 2))
            r = Fraction(-n * (N - n), N)
        else:
            C = _vp(3 * x ** 4 + b2 * x ** 3 + 3 * b4 * x * x
                    + 3 * b6 * x + b8, p)
            r = Fraction(-2 * B, 3) if C >= 3 * B else Fraction(-C, 4)
        if r:
            total += mpf(r.numerator) / mpf(r.denominator) * mp.log(p)
    total += mp.log(den)
    return total


def canonical_height(E: EllipticCurve, P, prec: int = 60) -> mpf:
    """The Neron-Tate height hhat(P), as the sum of Silverman's local heights.

    Torsion points have hhat = 0 exactly, which we decide by the group law
    (Mazur bounds the order by 16) rather than by a numerical threshold.
    """
    if P is O or isinstance(P, INFINITY):
        return mpf(0)
    for n in range(1, 17):
        if E.multiply(P, n) is O or isinstance(E.multiply(P, n), INFINITY):
            return mpf(0)
    Emin, Q = _to_minimal(E, P)
    with mp.workdps(prec + 10):
        h = (_archimedean_local_height(Emin, Q, prec)
             + _non_archimedean_local_height(Emin, Q, prec))
        return mp.mpf(h)


def height_pairing(E: EllipticCurve, P, Q, prec: int = 60) -> mpf:
    """<P,Q> = (hhat(P+Q) - hhat(P) - hhat(Q)) / 2.

    The subtraction must happen at full working precision: hhat values are
    O(1) numbers agreeing to many digits, and differencing them at mpmath's
    default 15 digits would leave the pairing accurate to only ~1e-17.
    """
    with mp.workdps(prec + 10):
        return mp.mpf((canonical_height(E, E.add(P, Q), prec)
                       - canonical_height(E, P, prec)
                       - canonical_height(E, Q, prec)) / 2)


# ---------- point search ----------

def search_points(E: EllipticCurve, height_bound: int = 100) -> list:
    """Every rational point with x = a/b^2, |a| <= height_bound*b^2, b <= sqrt(height_bound).

    For each candidate x the curve has a rational y exactly when the quadratic
    discriminant 4x^3 + b2 x^2 + 2 b4 x + b6 is a square in Q.  Writing
    x = a/b^2, that discriminant is N/b^6 with N an integer, so the test is
    `isqrt(N)**2 == N` -- exact, never a float judgement.  O is not returned.
    """
    a1, a2, a3, a4, a6 = E.ainvs
    b2, b4, b6 = E.b2, E.b4, E.b6
    points = []
    seen = set()
    bmax = isqrt(height_bound)
    for b in range(bmax + 1):
        b2v = b * b
        for a in range(-height_bound * (b2v or 1), height_bound * (b2v or 1) + 1):
            N = 4 * a ** 3 + b2 * a * a * b2v + 2 * b4 * a * b2v * b2v \
                + b6 * b2v * b2v * b2v
            if N < 0:
                continue
            s = isqrt(N)
            if s * s != N:
                continue
            x = Fraction(a, b2v) if b else Fraction(a)
            # y = (+/- sqrt(N)/b^3 - a1 x - a3) / 2, cleared of denominators
            for sign in ((1,) if s == 0 else (1, -1)):
                if b:
                    y = Fraction(sign * s - a1 * a * b - a3 * b ** 3,
                                 2 * b ** 3)
                else:
                    y = Fraction(sign * s - a1 * a - a3, 2)
                P = (x, y)
                if E.is_on_curve(P) and P not in seen:
                    seen.add(P)
                    points.append(P)
    return points


# ---------- independence, regulator, rank ----------

def _gram_determinant(entries: list, tol: mpf) -> mpf:
    """det of the symmetric matrix given as a flat row-major list of mpf.

    Fraction-free Gaussian elimination is not available over the reals, so we
    pivot on magnitude; `tol` is the absolute pivot size below which the matrix
    is declared singular, and then det() is the exact 0 it computed.  The
    caller scales `tol` by the largest entry (see INDEPENDENCE_TOLERANCE).
    """
    n = isqrt(len(entries))
    m = [list(row) for row in (entries[i * n:(i + 1) * n] for i in range(n))]
    det = mpf(1)
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) <= tol:
            return mpf(0)
        if piv != col:
            m[col], m[piv] = m[piv], m[col]
            det = -det
        det *= m[col][col]
        inv = 1 / m[col][col]
        for r in range(col + 1, n):
            f = m[r][col] * inv
            if f:
                for c in range(col, n):
                    m[r][c] -= f * m[col][c]
    return det


def independent_points(E: EllipticCurve, points, prec: int = 60) -> list:
    """A maximal subset of `points` independent modulo torsion.

    Points are tried in order of increasing canonical height; each is kept
    when the height pairing matrix of the kept set together with it is
    nonsingular at INDEPENDENCE_TOLERANCE.  Torsion points (and O) are
    dropped first -- decided by the group law, not by a height threshold.
    """
    live = []
    for P in points:
        if P is O or isinstance(P, INFINITY):
            continue
        if any(E.multiply(P, n) is O or isinstance(E.multiply(P, n), INFINITY)
               for n in range(1, 17)):
            continue
        live.append(P)
    live.sort(key=lambda P: canonical_height(E, P, prec))
    kept = []
    with mp.workdps(prec + 10):
        for P in live:
            trial = kept + [P]
            n = len(trial)
            entries = [height_pairing(E, trial[i], trial[j], prec)
                       for i in range(n) for j in range(n)]
            scale = max(abs(e) for e in entries) if entries else mpf(1)
            if _gram_determinant(entries, INDEPENDENCE_TOLERANCE * scale) != 0:
                kept = trial
    return kept


def regulator(E: EllipticCurve, generators, prec: int = 60) -> mpf:
    """det of the height pairing matrix <P_i, P_j>; 1 for the empty list (rank 0)."""
    gens = [P for P in generators
            if P is not O and not isinstance(P, INFINITY)]
    n = len(gens)
    if n == 0:
        return mpf(1)
    with mp.workdps(prec + 10):
        entries = [height_pairing(E, gens[i], gens[j], prec)
                   for i in range(n) for j in range(n)]
        return _gram_determinant(entries, mpf(0))


def rank_lower_bound(E: EllipticCurve, height_bound: int = 100,
                     prec: int = 60) -> tuple:
    """A LOWER BOUND only: naive search can never prove the rank is not higher.

    Returns (rank_found, generators): a maximal independent subset of the
    points found by ``search_points`` at this bound.  For curves whose
    generators are outside the searched box the true rank is strictly larger.
    """
    gens = independent_points(E, search_points(E, height_bound), prec)
    return len(gens), gens
