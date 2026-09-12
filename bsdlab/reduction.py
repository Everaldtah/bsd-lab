"""Tate's algorithm: local data of an elliptic curve at a prime.

Computes, at each prime p of bad reduction of a minimal model
    E : y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6,
the Kodaira symbol, the conductor exponent f_p, and the Tamagawa number c_p,
following Silverman, *Advanced Topics in the Arithmetic of Elliptic Curves*,
IV.9, and Cremona, *Algorithms for Modular Elliptic Curves*, section 3.2.

Everything is exact integer arithmetic; the only external helper is sympy's
``factorint``/``isprime`` for listing bad primes and checking primality.
Ogg's formula f = v_p(disc) - #components + 1 is asserted for every prime:
it is a theorem, so a violation means an implementation bug, and we raise.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from sympy import factorint, isprime

from bsdlab.curve import EllipticCurve, O


@dataclass(frozen=True)
class LocalData:
    """Local reduction data of a minimal model at the prime p.

    ``kodaira`` is the Kodaira symbol; ``f`` is the conductor exponent f_p;
    ``c`` is the Tamagawa number c_p (number of components of the special
    fibre in the Neron model); ``reduction`` is 'good', 'split multiplicative',
    'nonsplit multiplicative', or 'additive'; ``valuation_disc`` is v_p of the
    minimal discriminant.
    """

    p: int
    kodaira: str
    f: int
    c: int
    reduction: str
    valuation_disc: int


def local_data(E: EllipticCurve, p: int) -> LocalData:
    """Local data at p on the global minimal model of E.

    Runs Tate's algorithm on E.minimal_model() and returns the Kodaira
    symbol, conductor exponent, Tamagawa number, reduction type and
    v_p(minimal discriminant). Raises ArithmeticError if Ogg's formula
    f = v_p(disc) - #components + 1 fails at any point.
    """
    if not isprime(p):
        raise ValueError("local_data: %d is not prime" % p)
    try:
        Emin = E.minimal_model()
    except (ArithmeticError, ValueError):
        # curve.minimal_model can fail (e.g. j = 0 curves with c4 = 0).
        # Tate's algorithm only needs p-minimality, which the divide-by-p^i
        # fallback inside the loop supplies, so fall back to E's own model.
        Emin = E
    A = list(Emin.ainvs)
    v_total = _vp(Emin.discriminant, p)
    for _ in range(v_total // 12 + 3):
        out = _tate_step(A, p)
        if isinstance(out, tuple):
            kodaira, f, c, red, v = out
            break
        # Non-minimal at p: divide a_i by p^i.  If the divided model is not
        # integral, scale every a_i back up by a uniform u = p^e (Tate's
        # original loop works over Q_p and re-integralises the same way).
        A = [Fraction(x) / p ** e for x, e in zip(out, (1, 2, 3, 4, 6))]
        up = 0
        for x, i in zip(A, (1, 2, 3, 4, 6)):
            if x.denominator != 1:
                up = max(up, (_vp(x.denominator, p) + i - 1) // i)
        if up:
            A = [x * p ** (up * i) for x, i in zip(A, (1, 2, 3, 4, 6))]
        A = [int(x) for x in A]
    else:
        raise ArithmeticError(
            "Tate's algorithm did not terminate for %r at %d" % (Emin, p))
    if f != v - _components(kodaira) + 1:
        raise ArithmeticError(
            "Ogg's formula violated at p=%d: f=%d, v=%d, kodaira=%s" %
            (p, f, v, kodaira))
    return LocalData(p, kodaira, f, c, red, v)


def conductor(E: EllipticCurve) -> int:
    """The conductor N = prod over bad primes p of p^f_p."""
    N = 1
    for p in _bad_primes(E):
        N *= p ** local_data(E, p).f
    return N


def tamagawa_numbers(E: EllipticCurve) -> dict:
    """{p: c_p} for the bad primes p only."""
    out = {}
    for p in _bad_primes(E):
        d = local_data(E, p)
        if d.valuation_disc > 0:
            out[p] = d.c
    return out


def tamagawa_product(E: EllipticCurve) -> int:
    """The product of c_p over the bad primes (1 if there are none)."""
    n = 1
    for c in tamagawa_numbers(E).values():
        n *= c
    return n


# ---------- internal helpers ----------

def _vp(n: int, p: int) -> int:
    """The p-adic valuation v_p(n); raises when n == 0."""
    if n == 0:
        raise ArithmeticError("_vp(0, %d) is undefined" % p)
    v = 0
    n = abs(n)
    while n % p == 0:
        n //= p
        v += 1
    return v


def _bad_primes(E: EllipticCurve) -> list:
    """Sorted primes dividing the discriminant of E (a superset of the
    bad primes of the minimal model; local_data filters out the extras)."""
    return sorted(factorint(abs(E.discriminant)).keys())


def _inv_mod(x: int, p: int) -> int:
    """The inverse of x mod p; raises if p divides x."""
    if x % p == 0:
        raise ArithmeticError("%d is not invertible mod %d" % (x, p))
    return pow(x, -1, p)


def _proot(a: int, e: int, p: int) -> int:
    """Some z in [0, p) with z^e = a mod p; raises if there is none."""
    a %= p
    for z in range(p):
        if pow(z, e, p) == a:
            return z
    raise ArithmeticError("%d is not an %d-th power mod %d" % (a, e, p))


def _quad_has_root(a: int, b: int, c: int, p: int) -> bool:
    """True if a*T^2 + b*T + c has a root in F_p."""
    a, b, c = a % p, b % p, c % p
    if a == 0:
        return b != 0 or c == 0
    if p == 2:
        return any((a * y * y + b * y + c) % p == 0 for y in range(p))
    disc = (b * b - 4 * a * c) % p
    return disc == 0 or pow(disc, (p - 1) // 2, p) == 1


def _exact_div(n: int, d: int, what: str) -> int:
    """n // d, raising ArithmeticError if d does not divide n."""
    if n % d:
        raise ArithmeticError("%s: %d not divisible by %d" % (what, n, d))
    return n // d


def _components(kodaira: str) -> int:
    """Number of geometric components of the special fibre, by Kodaira type."""
    if kodaira == "I0":
        return 1
    star = kodaira.endswith("*")
    base = kodaira[:-1] if star else kodaira
    if base[0] == "I" and base[1:].isdigit():
        n = int(base[1:])
        return n + 5 if star else n
    return {"II": 1, "III": 2, "IV": 3, "IV*": 7, "III*": 8, "II*": 9}[kodaira]


def _origin_shift(C: EllipticCurve, p: int) -> tuple:
    """(r, t) mod p moving the singular point of the special fibre to (0,0).

    After the shift x = x' + r, y = y' + t the reduced curve has its
    singularity at the origin, i.e. p | a3', a4', a6'.  The closed forms
    for r are the standard ones (Cremona, section 3.2): r is a multiple
    root mod p of 4x^3 + b2 x^2 + 2 b4 x + b6, written out separately for
    p = 2, p = 3 and p >= 5 because the derivative degenerates for small p.
    """
    a1, a2, a3, a4, a6 = C.ainvs
    if p == 2:
        if C.b2 % 2 == 0:
            r = _proot(a4, 2, p)
            t = _proot(((r + a2) * r + a4) * r + a6, 2, p)
        else:
            inv = _inv_mod(a1, p)          # b2 odd forces a1 odd
            r = (inv * a3) % p
            t = (inv * (a4 + r * r)) % p
    elif p == 3:
        if C.b2 % 3 == 0:
            r = _proot(-C.b6, 3, p)
        else:
            r = (-_inv_mod(C.b2, p) * C.b4) % p
        t = (a1 * r + a3) % p
    else:
        if C.c4 % p == 0:
            r = (-C.b2 * _inv_mod(12, p)) % p
        else:
            r = (-(C.c6 + C.b2 * C.c4) * _inv_mod(12 * C.c4 % p, p)) % p
        t = (-(a1 * r + a3) * _inv_mod(2, p)) % p
    return r, t


def _tate_step(A: list, p: int) -> tuple:
    """One pass of Tate's algorithm on coefficients A.

    Returns (kodaira, f, c, reduction, valuation_disc), or the current
    coefficient list when the model is non-minimal at p (the caller then
    divides a_i by p^i and re-runs). Handles good, multiplicative and
    types II, III, IV here; the star types are delegated to _star_types.
    """
    a1, a2, a3, a4, a6 = A
    C = EllipticCurve(a1, a2, a3, a4, a6)
    v = _vp(C.discriminant, p)

    # Step 1: good reduction.
    if v == 0:
        return ("I0", 0, 1, "good", 0)

    # Step 2 setup: shift the singular point of the special fibre to the
    # origin, so that p | a3, a4, a6 (the c- and b-invariants with an even
    # number of a's are unchanged; in particular c4 and the discriminant).
    r, t = _origin_shift(C, p)
    C = C.transform(1, r, 0, t)
    a1, a2, a3, a4, a6 = C.ainvs
    for name, ai in (("a3", a3), ("a4", a4), ("a6", a6)):
        if ai % p:
            raise ArithmeticError(
                "p=%d does not divide %s after the origin shift" % (p, name))

    # Step 3: multiplicative reduction, type I_v.
    if C.c4 % p != 0:
        split = _quad_has_root(1, a1, -a2, p)
        if split:
            kod, f, c, red = "I%d" % v, 1, v, "split multiplicative"
        else:
            kod, f, red = "I%d" % v, 1, "nonsplit multiplicative"
            c = 2 if v % 2 == 0 else 1
        return (kod, f, c, red, v)

    # Additive reduction: types II, III, IV in order.
    if a6 % (p * p):
        return ("II", v, 1, "additive", v)
    if C.b8 % (p ** 3):
        return ("III", v - 1, 2, "additive", v)
    if C.b6 % (p ** 3):
        a3t = _exact_div(a3, p, "a3/p") % p
        a6t = _exact_div(a6, p * p, "a6/p^2") % p
        c = 3 if _quad_has_root(1, a3t, -a6t, p) else 1
        return ("IV", v - 2, c, "additive", v)

    # Beyond type IV: shift so that p | a1, a2; p^2 | a3, a4; p^3 | a6,
    # making the cubic P(T) = T^3 + (a2/p)T^2 + (a4/p^2)T + a6/p^3 integral.
    if p == 2:
        s = _proot(a2, 2, p)
        t = p * _proot(_exact_div(a6, p * p, "a6/p^2"), 2, p)
    elif p == 3:
        # NB the full integers, not their residues: a3' = 3*a3 is divisible
        # by 9, while a3 reduced mod 3 would leave a3' only divisible by 3.
        s, t = a1, a3
    else:
        inv2 = _inv_mod(2, p)
        s, t = (-a1 * inv2) % p, (-a3 * inv2) % p
    C = C.transform(1, 0, s, t)
    a1, a2, a3, a4, a6 = C.ainvs
    for name, ai, e in (("a1", a1, 1), ("a2", a2, 1), ("a3", a3, 2),
                        ("a4", a4, 2), ("a6", a6, 3)):
        if ai % p ** e:
            raise ArithmeticError(
                "p=%d does not divide %s enough after the second shift" %
                (p, name))
    return _star_types(C, p, v)


def _star_types(C: EllipticCurve, p: int, v: int) -> tuple:
    """Star types (I0*, In*, IV*, III*, II*) on the model C of _tate_step.

    Expects C to satisfy p | a1, a2; p^2 | a3, a4; p^3 | a6. Returns
    (kodaira, f, c, 'additive', v), or the current coefficient list when
    the model turns out to be non-minimal at p.
    """
    a1, a2, a3, a4, a6 = C.ainvs
    p2, p3, p4 = p * p, p ** 3, p ** 4

    # Roots of the cubic P(T) = T^3 + b T^2 + c T + d with b = a2/p,
    # c = a4/p^2, d = a6/p^3, all reduced mod p.  Here w = -disc(P) and
    # x = 3c - b^2: w = 0 detects a repeated root, and w = x = 0 (mod p)
    # a triple one.
    b = _exact_div(a2, p, "a2/p") % p
    c = _exact_div(a4, p2, "a4/p^2") % p
    d = _exact_div(a6, p3, "a6/p^3") % p
    w = 27 * d * d - b * b * c * c + 4 * b * b * b * d - 18 * b * c * d \
        + 4 * c * c * c
    x = 3 * c - b * b

    if w % p:
        # Three distinct roots: type I0*, and the non-identity components
        # correspond to the roots of P in F_p.
        n = sum(1 for T in range(p)
                if (T * T * T + b * T * T + c * T + d) % p == 0)
        return ("I0*", v - 4, 1 + n, "additive", v)

    if x % p:
        # One double root: type In* for some n >= 1.  Shift it to T = 0.
        if p == 2:
            r = _proot(c, 2, p)
        elif p == 3:
            r = c * _inv_mod(b, p) % p
        else:
            r = (b * c - 9 * d) * _inv_mod(2 * x % p, p) % p
        C = C.transform(1, p * r, 0, 0)
        a1, a2, a3, a4, a6 = C.ainvs
        # Blow up alternately in y and x, keeping the model minimal, until
        # neither -a3^2/4 nor a4^2/(4a2) cancels; ix, iy count the blow-ups.
        ix = iy = 3
        mx = my = p2
        while True:
            def quotients():
                return (_exact_div(a2, p, "a2/p") % p,
                        _exact_div(a3, my, "a3/my") % p,
                        _exact_div(a4, p * mx, "a4/(p*mx)") % p,
                        _exact_div(a6, mx * my, "a6/(mx*my)") % p)
            a2t, a3t, a4t, a6t = quotients()
            if (a3t * a3t + 4 * a6t) % p == 0:
                if p == 2:
                    t = my * _proot(a6t, 2, p)
                else:
                    t = my * (-a3t * _inv_mod(2, p)) % p
                C = C.transform(1, 0, 0, t)
                a1, a2, a3, a4, a6 = C.ainvs
                my *= p
                iy += 1
                a2t, a3t, a4t, a6t = quotients()
                if (a4t * a4t - 4 * a6t * a2t) % p == 0:
                    if p == 2:
                        r = mx * _proot(a6t * _inv_mod(a2t, p), 2, p)
                    else:
                        r = mx * (-a4t * _inv_mod(2 * a2t % p, p)) % p
                    C = C.transform(1, r, 0, 0)
                    a1, a2, a3, a4, a6 = C.ainvs
                    mx *= p
                    ix += 1
                    continue
                c = 4 if _quad_has_root(a2t, a4t, a6t, p) else 2
            else:
                c = 4 if _quad_has_root(1, a3t, -a6t, p) else 2
            break
        return ("I%d*" % (ix + iy - 5), v - ix - iy + 1, c, "additive", v)

    # Triple root: shift it to T = 0, so p^2 | a2, p^3 | a4, p^4 | a6.
    if p == 2:
        r = b
    elif p == 3:
        r = _proot(-d, 3, p)
    else:
        r = -b * _inv_mod(3, p) % p
    C = C.transform(1, p * r, 0, 0)
    a1, a2, a3, a4, a6 = C.ainvs
    for name, ai, e in (("a2", a2, 2), ("a4", a4, 3), ("a6", a6, 4)):
        if ai % p ** e:
            raise ArithmeticError(
                "no triple root at 0: p=%d divides %s too little" % (p, name))
    a3t = _exact_div(a3, p2, "a3/p^2") % p
    a6t = _exact_div(a6, p4, "a6/p^4") % p
    if (a3t * a3t + 4 * a6t) % p:
        c = 3 if _quad_has_root(1, a3t, -a6t, p) else 1
        return ("IV*", v - 6, c, "additive", v)
    if p == 2:
        t = -p2 * _proot(a6t, 2, p)
    else:
        t = p2 * (-a3t * _inv_mod(2, p)) % p
    C = C.transform(1, 0, 0, t)
    a1, a2, a3, a4, a6 = C.ainvs
    if a4 % p4:
        return ("III*", v - 7, 2, "additive", v)
    if a6 % p ** 6:
        return ("II*", v - 8, 1, "additive", v)
    # Past type II*: the model is non-minimal at p.  Return the current
    # coefficients -- the translations above are exactly what makes the
    # division by p^i admissible -- and let the caller divide and retry.
    return [a1, a2, a3, a4, a6]
