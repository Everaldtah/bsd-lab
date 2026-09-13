"""Descent via 2-isogeny for elliptic curves with rational 2-torsion.

Moves a rational point of order 2 to (0,0), giving
    E : y^2 = x^3 + a x^2 + b x,
with 2-isogenous curve
    E' : Y^2 = X^3 - 2a X^2 + b' X,   b' = a^2 - 4b.
The descent map sends E'(Q)/phi(E(Q)) into Q*/(Q*)^2, represented by
squarefree divisors d of b (resp. b'); the class d is in the image iff
    H_d : N^2 = d M^4 + a M^2 e^2 + (b/d) e^4
has a nontrivial rational point.
"""
from fractions import Fraction
from math import gcd, isqrt

from bsdlab.curve import INFINITY


def _apply_transform(ainvs, u, r, s, t):
    """Silverman's change x=u^2 X+r, y=u^3 Y+u^2 s X+t on (a1..a6), Fractions."""
    a1, a2, a3, a4, a6 = (Fraction(v) for v in ainvs)
    return (
        (a1 + 2 * s) / u,
        (a2 - s * a1 + 3 * r - s * s) / (u * u),
        (a3 + r * a1 + 2 * t) / (u ** 3),
        (a4 - s * a3 + 2 * r * a2 - (r * s + t) * a1 + 3 * r * r - 2 * s * t) / (u ** 4),
        (a6 + r * a4 + r * r * a2 + r ** 3 - t * a3 - t * t - r * t * a1) / (u ** 6),
    )


def two_torsion_model(E):
    """Return (a, b, to_E) for E ~ y^2 = x^3 + a x^2 + b x with integers a, b.

    to_E is a callable mapping a point (X, Y) on the reduced model back to a
    point on E.  Returns None when E has no rational point of order 2 -- the
    only failure mode, so the reason is always "no rational 2-torsion".
    """
    two_tors = [P for P in E.torsion_points()
                if P != INFINITY and E.multiply(P, 2) == INFINITY]
    if not two_tors:
        return None
    x0, y0 = (Fraction(v) for v in two_tors[0])
    # Move P to (0,0): shift x by x0, y by y0.  Then a6 = 0 (P on curve) and
    # a3 = 0 (tangent at P horizontal since 2P = O).
    ainvs = _apply_transform(E.ainvs, Fraction(1), x0, Fraction(0), y0)
    # Kill a1 via y -> Y - (a1/2) X, leaving y^2 = x^3 + a2 x^2 + a4 x.
    a1 = ainvs[0]
    ainvs = _apply_transform(ainvs, Fraction(1), Fraction(0), -a1 / 2, Fraction(0))
    a2, a4 = ainvs[1], ainvs[3]
    if ainvs[2] != 0 or ainvs[4] != 0:
        raise NotImplementedError("2-torsion move failed: a3 or a6 nonzero")
    # Clear denominators: x = X / L^2, y = Y / L^3 scales a2 -> a2*L^2, a4 -> a4*L^4.
    L = a2.denominator * a4.denominator  # not minimal, but integral suffices
    a, b = a2 * L * L, a4 * L ** 4
    if b == 0:
        raise NotImplementedError("degenerate b = 0 (singular reduced curve)")

    def to_E(P):
        if P == INFINITY:
            return INFINITY
        X, Y = P
        x = Fraction(X) / (L * L) + x0
        y = Fraction(Y) / (L ** 3) - a1 / 2 * (x - x0) + y0
        return (x, y)

    return (int(a), int(b), to_E)


def selmer_phi(a, b, search_bound=1000):
    """Local-solvability Selmer bound and global-image data for H_d, d | b.

    Applies to E : y^2 = x^3 + a x^2 + b x; the descent map
    E'(Q)/phi(E(Q)) -> Q*/(Q*)^2 has image among squarefree divisors d of b.

    Prime set tested: p = 2 and p | b*b'.  These are exactly the primes of
    bad reduction for the pair (E, E'), since disc(E) = 16 b^2 b' and
    disc(E') = 256 b b'^2; for any p >= 3 outside them the cover
    C_d : N^2 = quartic has good reduction, so #C_d(F_p) >= p + 1 - 2*sqrt(p)
    > 0 by Hasse, and a nonsingular F_p-point lifts to Q_p by Hensel.
    """
    b_prime = a * a - 4 * b
    bad = {2}
    for n in (abs(b), abs(b_prime)):
        p = 2
        while p * p <= n:
            if n % p == 0:
                bad.add(p)
                while n % p == 0:
                    n //= p
            p += 1
        if n > 1:
            bad.add(n)
    selmer, glob, undec, pts = [], [], [], {}
    for d in _squarefree_divisors(b):
        bod = b // d
        if not _solvable_real(d, a, bod):
            continue
        if not all(_hensel_solvable(d, a, bod, p, 3 + 2 * (1 if p == 2 else 0))
                   for p in bad):
            continue
        selmer.append(d)
        found = _has_global_point(d, a, bod, search_bound)
        if found:
            glob.append(d)
            pts[d] = found
        else:
            undec.append(d)
    return {'selmer': sorted(selmer), 'global': sorted(glob),
            'undecided': sorted(undec), 'points': pts}


def _pow2_at_least(g):
    """Smallest power of 2 >= max(g, 1)."""
    return 1 << (max(g, 1) - 1).bit_length()


def rank_bounds(E, search_bound=1000):
    """2-isogeny rank bounds for E; all None if no rational 2-torsion.

    rank_upper = log2 |Sel^phi| + log2 |Sel^phihat| - 2, since
    2^rank = |Sel^phi| |Sel^phihat| / (4 |Sha[phi]| |Sha[phihat]|).
    rank_lower uses the found global classes: their closure in Q*/(Q*)^2 is a
    subgroup of the true image, so has order p2(g) >= g (a power of 2), and
    2^rank = |im^phi| |im^phihat| / 4 >= p2(g1) p2(g2) / 4.
    sha_phi_order_bound = |Sel^phi| / p2(g1): an UPPER bound on |Sha[phi]|
    (equal to it iff the found classes exhaust the image).
    """
    model = two_torsion_model(E)
    if model is None:
        return {'rank_lower': None, 'rank_upper': None,
                'sha_phi_order_bound': None, 'certified': None,
                'reason': 'E has no rational point of order 2; '
                          'descent via 2-isogeny does not apply'}
    a, b, _ = model
    bp = a * a - 4 * b
    ph = selmer_phi(a, b, search_bound)
    phh = selmer_phi(-2 * a, bp, search_bound)
    s1, s2 = len(ph['selmer']), len(phh['selmer'])
    g1, g2 = len(ph['global']), len(phh['global'])
    for s in (s1, s2):
        if s < 1 or s & (s - 1):
            raise RuntimeError(f"Selmer size {s} not a power of 2: local-test bug")
    prod = s1 * s2
    if prod < 4 or prod % 4 or (q := prod // 4) & (q - 1):
        raise RuntimeError(f"Selmer product {prod} not 4 * power of 2: local-test bug")
    rank_upper = q.bit_length() - 1
    lower_prod = _pow2_at_least(g1) * _pow2_at_least(g2)
    rank_lower = max(0, (lower_prod // 4).bit_length() - 1)
    sha_bound = s1 // _pow2_at_least(g1)
    certified = rank_lower == rank_upper
    reason = (f"2-isogeny descent on y^2=x^3+({a})x^2+({b})x: "
              f"|Sel^phi|={s1}, |Sel^phihat|={s2} -> rank<={rank_upper}; "
              f"{g1}+{g2} global classes found (bound {search_bound}) -> "
              f"rank>={rank_lower}; undecided d: phi={ph['undecided']}, "
              f"phihat={phh['undecided']}; certified={certified}")
    return {'rank_lower': rank_lower, 'rank_upper': rank_upper,
            'sha_phi_order_bound': sha_bound, 'certified': certified,
            'reason': reason}


def _squarefree_divisors(b):
    """Squarefree divisors d (up to sign) of the nonzero integer b."""
    n = abs(b)
    primes = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            primes.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        primes.append(n)
    out = []
    for i in range(1 << len(primes)):
        d = 1
        for j, q in enumerate(primes):
            if (i >> j) & 1:
                d *= q
        out.append(d)
        out.append(-d)
    return sorted(out)


def _solvable_real(d, a, b_over_d):
    """Real solvability of N^2 = d M^4 + a M^2 e^2 + b_over_d e^4.

    Write s = M^2/e^2 >= 0 (with e=0 giving s=+inf, value d): the quartic
    takes the values d*s^2 + a*s + b_over_d for s in [0, inf].  A solution
    over R exists iff this quadratic is >= 0 somewhere on [0, inf]:
      - d > 0        : s -> inf gives +inf, solvable;
      - b_over_d > 0 : s = 0 gives b_over_d > 0, solvable;
      - d < 0 and b_over_d < 0 (so b > 0): both ends negative, and the
        downward parabola has its vertex at s* = -a/(2d) >= 0 only when
        a >= 0, with vertex value (4b - a^2)/(4d) >= 0 iff a^2 >= 4b.
        Solvable exactly when a > 0 and a^2 >= 4b (at a = 0 the vertex
        value is b < 0), i.e. when b' = a^2 - 4b >= 0 with a > 0.
    """
    if d > 0 or b_over_d > 0:
        return True
    return a > 0 and a * a >= 4 * d * b_over_d


def _is_square_mod_pk(x, p, k):
    """Is there n with n^2 = x mod p^k?  x integer, p odd."""
    m = p ** k
    x %= m
    if x == 0:
        return True
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    if v % 2:
        return False
    return pow(x, (p - 1) // 2, p) == 1  # unit square mod p^k iff QR mod p


def _quartic_has_square(c4, c2, c0, p, k):
    """Does h(t) = c4 t^4 + c2 t^2 + c0 take a square value mod p^k, t in Z/p^k?

    A point on N^2 = dM^4 + aM^2e^2 + (b/d)e^4 with e (resp. M) a p-adic unit
    gives t = M/e (resp. s = e/M) in Z_p with h(t) a square mod p^k.  For odd
    p, unit values h(t0) mod p are decided by a Legendre symbol; only classes
    with p | h(t0) need refining to mod p^k.  For p = 2, k <= 6, brute force.
    """
    h = lambda t: c4 * t % p ** k * t % p ** k * t % p ** k * t + c2 * t % p ** k * t + c0
    if p == 2:
        m = 1 << k
        sq = {n * n % m for n in range(m)}
        return any(h(t) % m in sq for t in range(m))
    for t0 in range(p):
        h0 = h(t0)
        if h0 % p:
            if pow(h0 % p, (p - 1) // 2, p) == 1:
                return True
        else:
            p2 = p * p
            for j in range(p2):
                if _is_square_mod_pk(h(t0 + p * j), p, k):
                    return True
    return False


def _hensel_solvable(d, a, b_over_d, p, k):
    """Q_p solvability of N^2 = dM^4 + aM^2e^2 + (b/d)e^4 via mod-p^k search.

    A primitive solution has e a unit (then t = M/e and dt^4 + at^2 + b/d is a
    square in Z_p) or M a unit (then s = e/M and (b/d)s^4 + as^2 + d is a
    square in Z_p).  k = 3 + 2*v_p(2) makes mod-p^k solvability equivalent to
    Q_p solvability for these quartics (k = 5 at p = 2, k = 3 at odd p).
    """
    return (_quartic_has_square(d, a, b_over_d, p, k)
            or _quartic_has_square(b_over_d, a, d, p, k))


def _has_global_point(d, a, b_over_d, bound):
    """Search coprime M, e in a box for N^2 = d M^4 + a M^2 e^2 + (b/d) e^4.

    Returns (M, e, N) or None; None proves nothing (record as UNDECIDED).
    Only M^2, e^2 occur, so it suffices to scan M, e >= 0 with gcd(M, e) = 1.
    """
    from math import isqrt as _isqrt
    for e in range(bound + 1):
        e2 = e * e
        e4 = e2 * e2
        for M in range(bound + 1):
            if e == 0 and M == 0:
                continue
            if e and M and gcd(M, e) != 1:
                continue
            val = d * M ** 4 + a * M * M * e2 + b_over_d * e4
            if val >= 0:
                s = _isqrt(val)
                if s * s == val:
                    return (M, e, s)
    return None
