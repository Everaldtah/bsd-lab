"""Exact rational modular symbols for Gamma_0(N) (E4 machinery, part 2).

Computes L(E(x)chi_d, 1)/Omega_d as an EXACT Fraction from Manin symbols
for Gamma_0(N) (Stein, *Modular Forms, a Computational Approach*, ch. 8;
Cremona, *Algorithms for Modular Elliptic Curves*, ch. II).  No floating
point anywhere on the output path; numpy floats bootstrap the Hecke
eigenvector only, which is then rationalised and re-verified exactly.

Pipeline:
  1. Manin symbols: pairs (c:d) in P^1(Z/NZ) with gcd(c,d,N)=1.
  2. Manin relations  x + x*sigma = 0  and  x + x*tau + x*tau^2 = 0
     (sigma = [[0,-1],[1,0]], tau = [[0,-1],[1,-1]], right action
     (u:v).g = (u,v) . g on row vectors); quotient by the relation
     lattice over Q via rref.
  3. Hecke T_l (l not dividing N) on the quotient via Merel's Heilbronn
     matrices (a > b >= 0, d > c >= 0, ad - bc = l).  NOTE: the naive
     l+1 double-coset reps [[1,i],[0,l]] + [[l,0],[0,1]] act on PATH
     symbols but do NOT preserve the Manin relation lattice; using them
     yields non-rational "eigenvalues" and was the first implementation
     bug here.
  4. The star involution * : {a,b} -> {-a,-b} (matrix [[-1,0],[0,1]])
     splits the 2-dimensional newform eigenspace; the eigenfunctionals
     lambda^+/- (contragredient, row vectors) are bootstrapped in floats
     and re-verified EXACTLY (asserts lambda*T_l = a_l*lambda for
     l = 2,3,5,7,11,13 with l not dividing N, and lambda*star =
     +-lambda).  The Fricke W_N is not needed: the star parity of chi_d
     alone selects the branch (the mismatched parity sums to exact 0).
  5. Path expansion of {0, b/m} along the continued-fraction convergents
     of b/m (consecutive convergents are Farey neighbours, so each piece
     is directly a Manin symbol up to sign).
  6. S(d) = sum_{b mod m} chi_d(b) * lambda({0, b/m}) (raw_twisted_sum),
     exact but in the eigenvector's arbitrary scale.  period_scales fixes
     that scale once per (curve, sign) against Omega^+/- of E, giving the
     exact L(E,chi_d,1)*sqrt|d|/Omega^eps (algebraic_twisted_l).
     HISTORY: the first version pinned one constant C to 16 anchor values;
     on 60 unseen twists it was wrong by a factor 2 or 4 on 10 of them
     (the eigenvector scale is curve- and sign-dependent).  Anchors alone
     cannot validate a normaliser -- scratch/oos_modsym.py is the gate.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd

from sympy.functions.combinatorial.numbers import jacobi_symbol

from bsdlab.curve import EllipticCurve

__all__ = ["ManinSpace", "twisted_l_ratio", "algebraic_twisted_l",
           "raw_twisted_sum", "period_scales", "omega_pm", "kron",
           "heilbronn_matrices"]


def _units(n: int) -> list:
    """The units of Z/nZ."""
    return [u for u in range(1, n) if gcd(u, n) == 1] + ([1] if n == 1 else [])


def p1_symbols(n: int) -> list:
    """Representatives of P^1(Z/nZ): classes (c:d), gcd(c,d,n) = 1 mod units.

    Iterating pairs in lexicographic order, the first unseen pair of a unit
    orbit is the canonical representative; the whole orbit is then marked
    seen, so the total work is #classes * phi(n) plus the scan.
    """
    units = _units(n)
    seen = set()
    out = []
    for c in range(n):
        for d in range(n):
            if gcd(gcd(c, d), n) != 1 or (c, d) in seen:
                continue
            out.append((c, d))
            for u in units:
                seen.add((c * u % n, d * u % n))
    return out


def _rrel(u: int, v: int, m: tuple, n: int) -> tuple:
    """Right action of the integer matrix m = (a, b, c, d) on the Manin
    symbol (u:v): row (u,v) times m, reduced mod n but NOT canonicalised."""
    a, b, c, d = m
    return ((u * a + v * c) % n, (u * b + v * d) % n)


SIGMA = (0, -1, 1, 0)    # order 4 in SL_2(Z); relation x + x*sigma = 0
TAU = (0, -1, 1, -1)     # order 3 up to -I; relation x + x*tau + x*tau^2 = 0
STAR = (-1, 0, 0, 1)     # det -1, induces * : {a,b} -> {-a,-b}


def _manin_relations(symbols: list, n: int, index: dict) -> list:
    """Manin relation rows as sparse vectors {symbol_index: coeff}: for each
    symbol x, the rows x + x*sigma and x + x*tau + x*tau^2 (all coeffs 1)."""
    out = []
    for (u, v) in symbols:
        for ms, length in ((SIGMA, 2), (TAU, 3)):
            vec = {}
            x, y = u, v
            for _ in range(length):
                i = index[(x, y)]
                vec[i] = vec.get(i, 0) + 1
                x, y = _rrel(x, y, ms, n)
            out.append(vec)
    return [v for v in out if v]


def _rref_sparse(rows: list, ncols: int) -> tuple:
    """Sparse rref over Fraction.  Returns (pivot_rows, pivot_cols, free_cols);
    pivot_rows[i] is normalised with a 1 in pivot_cols[i] and zeros in the
    other pivot columns."""
    rows = [dict((c, Fraction(v)) for c, v in r.items()) for r in rows]
    pivot_rows, pivot_cols = [], []
    used = set()
    for col in range(ncols):
        pr = None
        for r in rows:
            if r.get(col) and id(r) not in used:
                pr = r
                break
        if pr is None:
            continue
        used.add(id(pr))
        inv = Fraction(1) / pr[col]
        for c in pr:
            pr[c] *= inv
        for r in rows:
            if r is pr or not r.get(col):
                continue
            f = r[col]
            for c, val in pr.items():
                new = r.get(c, Fraction(0)) - f * val
                if new:
                    r[c] = new
                elif c in r:
                    del r[c]
        pivot_rows.append(pr)
        pivot_cols.append(col)
    free = [c for c in range(ncols) if c not in pivot_cols]
    return pivot_rows, pivot_cols, free


class ManinSpace:
    """The space of weight-2 modular symbols for Gamma_0(N) over Q, as a
    quotient of the free space on Manin symbols by the Manin relations."""

    def __init__(self, N: int):
        self.N = N
        self.symbols = p1_symbols(N)
        self.index = {}
        units = _units(N)
        for i, (c, d) in enumerate(self.symbols):
            self.index[(c, d)] = i
            for u in units:
                self.index[(c * u % N, d * u % N)] = i
        rels = _manin_relations(self.symbols, N, self.index)
        pivot_rows, pivot_cols, free = _rref_sparse(rels, len(self.symbols))
        self.dim = len(free)
        self.basis_symbols = free
        pos = {f: j for j, f in enumerate(free)}
        # symcoord[i] = image of symbol i in the quotient basis (free cols).
        self.symcoord = []
        pivrow = {p: r for p, r in zip(pivot_cols, pivot_rows)}
        for i in range(len(self.symbols)):
            vec = [Fraction(0)] * self.dim
            if i in pos:
                vec[pos[i]] = Fraction(1)
            else:
                r = pivrow[i]
                for f, val in r.items():
                    if f in pos:
                        vec[pos[f]] = -val
            self.symcoord.append(vec)

    def coords(self, combo: dict) -> list:
        """Coordinates of a Q-linear combination {symbol_index: coeff}."""
        vec = [Fraction(0)] * self.dim
        for i, c in combo.items():
            sc = self.symcoord[i]
            for j, v in enumerate(sc):
                if v:
                    vec[j] += c * v
        return vec

    def _apply(self, m: tuple, combo: dict) -> dict:
        """Right action of an integer matrix m on a symbol combination."""
        out = {}
        for i, c in combo.items():
            u, v = self.symbols[i]
            j = self.index[_rrel(u, v, m, self.N)]
            out[j] = out.get(j, 0) + c
        return {i: c for i, c in out.items() if c}

    def free_symbols(self) -> list:
        """Indices of the quotient basis symbols (the free rref columns)."""
        return self.basis_symbols

    def op_matrix(self, mats: list) -> list:
        """Matrix of x -> sum_m x*m on quotient coords, as column vectors."""
        basis = self.free_symbols()
        cols = []
        for i in basis:
            combo = {}
            for m in mats:
                combo = {**combo}
                for j, c in self._apply(m, {i: 1}).items():
                    combo[j] = combo.get(j, 0) + c
            combo = {i: c for i, c in combo.items() if c}
            cols.append(self.coords(combo))
        return cols

    def hecke_cols(self, l: int) -> list:
        """T_l (l not dividing N) on quotient coords, as column vectors.

        Uses Merel's Heilbronn matrices (Stein, ch. `General Modular
        Symbols', Prop. 1.29: {(a b; c d) : a > b >= 0, d > c >= 0,
        ad - bc = n} satisfies condition C_n), each acting by right
        multiplication on (u:v) with summands dropped when the image
        fails gcd(u', v', N) = 1.  The naive p+1 double-coset reps do NOT
        preserve the Manin relation lattice (they act on path symbols,
        not Manin symbols); the Heilbronn set does.
        """
        if gcd(l, self.N) != 1:
            raise ValueError("l=%d divides N=%d" % (l, self.N))
        basis = self.free_symbols()
        cols = []
        for i in basis:
            combo = {}
            for m in heilbronn_matrices(l):
                u, v = self.symbols[i]
                x, y = _rrel(u, v, m, self.N)
                if gcd(gcd(x, y), self.N) != 1:
                    continue
                j = self.index[(x, y)]
                combo[j] = combo.get(j, 0) + 1
            cols.append(self.coords(combo))
        return cols

    def star_cols(self) -> list:
        """The star involution * : {a,b} -> {-a,-b} on quotient coords."""
        return self.op_matrix([STAR])

    def path_combo(self, b: int, m: int) -> dict:
        """{0, b/m} as a combination of Manin symbols, via the convergents
        of the continued fraction of b/m (consecutive convergents are Farey
        neighbours, hence each piece is exactly one Manin symbol).

        Convention: for g = [[a,b],[c,d]] in SL_2(Z) the symbol (c:d) is
        the path {g(0), g(inf)} = {b/d, a/c}.  For consecutive convergents
        c_i = p_i/q_i with det_i = p_i q_{i+1} - p_{i+1} q_i = +-1: if
        det_i = +1 the matrix [[p_i, p_{i+1}],[q_i, q_{i+1}]] lies in
        SL_2(Z) and gives {c_i, c_{i+1}} = -(q_i : q_{i+1}); if det_i = -1
        the column-swapped matrix does, giving +(q_{i+1} : q_i).
        """
        convs = _convergents(b, m)
        combo = {}
        # integer prefix: {0, b//m} = sum of {k-1, k}, each Farey with
        # symbol (1:1); needed when b >= m (only the m = 1 twist, where
        # the path is {0, 1} = {g(0), g(inf)} for g = [[0,1],[1,1]]).
        for k in range(1, b // m + 1):
            one = self.index[(1 % self.N, 1 % self.N)]
            combo[one] = combo.get(one, 0) + 1
        for (p_i, q_i), (p_j, q_j) in zip(convs, convs[1:]):
            det = p_i * q_j - p_j * q_i
            if det == 1:
                sym, sign = (q_i, q_j), -1
            else:
                sym, sign = (q_j, q_i), 1
            k = self.index[(sym[0] % self.N, sym[1] % self.N)]
            combo[k] = combo.get(k, 0) + sign
        return {i: c for i, c in combo.items() if c}


def _convergents(b: int, m: int) -> list:
    """Convergents (p, q) of b/m from 0/1 up to b/m, Euclid-style.

    Uses the subtractive form of the extended Euclidean algorithm so that
    every intermediate (p_i, q_i) appears with p_i q_{i+1} - p_{i+1} q_i
    = +-1 (Farey neighbours), which is exactly what path_combo needs.
    """
    convs = [(0, 1), (1, 0)]
    # convs[-1] is the infinity placeholder until the first digit is folded
    r, s = b, m
    while s:
        a = r // s
        r, s = s, r % s
        convs.append((a * convs[-1][0] + convs[-2][0],
                      a * convs[-1][1] + convs[-2][1]))
    # convs = [p_{-2}, p_{-1}=inf, p_0 = (b//m)/1, ..., p_r = b/m]
    assert convs[2] == (b // m, 1)
    return convs[2:]


HECKE_PRIMES = (2, 3, 5, 7, 11, 13)


def kron(a: int, n: int) -> int:
    """The full Kronecker symbol (a/n) for n != 0, including even n.

    (bsdlab.twists.kronecker bottoms out in jacobi_symbol, which rejects
    even moduli > 2; the twisted sums here need chi_d(b) for all b.)
    """
    a, n = int(a), int(n)
    if n == 0:
        raise ValueError("n == 0")
    out = 1
    if n < 0:
        n = -n
        if a < 0:
            out = -out
    while n % 2 == 0:
        n //= 2
        if a % 2 == 0:
            return 0
        out *= 1 if a % 8 in (1, 7) else -1
    if n == 1:
        return out
    return out * int(jacobi_symbol(a % n, n))


def heilbronn_matrices(n: int) -> list:
    """Matrices [[a,b],[c,d]] with a > b >= 0, d > c >= 0, ad - bc = n.

    Stein, *Modular Forms, a Computational Approach*, `General Modular
    Symbols', Prop. 1.29: this set satisfies Merel's condition C_n, so it
    may be used to compute T_n on Manin symbols.  For n = 2 it is exactly
    the four matrices (2 0;0 1), (1 0;0 2), (2 1;0 1), (1 0;1 2).
    """
    out = []
    for b in range(n + 1):
        for c in range(n + 1):
            t = n + b * c
            for a in range(b + 1, t + 1):
                if t % a:
                    continue
                d = t // a
                if d > c:
                    out.append((a, b, c, d))
    return out


def _dense(cols: list):
    """Column vectors -> numpy square array."""
    import numpy as np
    A = np.array([[float(x) for x in col] for col in cols])
    return A.T


def _rationalise(vec) -> list:
    """A float vector scaled so its dominant component is 1, as Fractions."""
    k = max(range(len(vec)), key=lambda i: abs(vec[i]))
    scaled = [x / vec[k] for x in vec]
    return [Fraction(x).limit_denominator(10 ** 9) for x in scaled]


def _row_action(v: list, cols: list) -> list:
    """The row vector v times the matrix given by column vectors:
    (v*M)_i = v . cols[i]."""
    n = len(v)
    out = []
    for i in range(n):
        col = cols[i]
        out.append(sum(v[j] * col[j] for j in range(n) if v[j] and col[j]))
    return out


def eigenfunctional(M: ManinSpace, E: EllipticCurve) -> dict:
    """Hecke eigenfunctionals lambda^+/- on the Manin space of E's level.

    Returns {+1: lam_plus, -1: lam_minus}: row vectors over Q with
    lambda*T_l = a_l(E)*lambda for the HECKE_PRIMES with l not dividing N,
    and lambda*star = sign*lambda.  numpy floats bootstrap the kernel of the
    stacked (T_l^T - a_l I, star^T - sign I); the rationalised vector is
    then re-verified by EXACT rational linear algebra (asserts).
    """
    import numpy as np
    N = M.N
    aps = {l: E.ap(l) for l in HECKE_PRIMES if l % N}
    T = {l: _dense(M.hecke_cols(l)) for l in aps}
    S = _dense(M.star_cols())
    out = {}
    for sign in (1, -1):
        rows = [T[l].T - float(a) * np.eye(M.dim) for l, a in aps.items()]
        rows.append(S.T - float(sign) * np.eye(M.dim))
        A = np.vstack(rows)
        _, sv, VT = np.linalg.svd(A)
        v = _rationalise(VT[-1])
        cols = {l: M.hecke_cols(l) for l in aps}
        scols = M.star_cols()
        for l, a in aps.items():
            lhs, rhs = _row_action(v, cols[l]), [a * x for x in v]
            assert lhs == rhs, "lambda*H_%d != a_%d*lambda (E=%r)" % (l, l, E)
        lhs = _row_action(v, scols)
        assert lhs == [sign * x for x in v], "lambda*star != %d*lambda" % sign
        out[sign] = v if sign == 1 else _primitive(v)
    return out


def _primitive(v: list) -> list:
    """Scale a rational vector to primitive integers, dominant entry > 0.

    The odd (star = -1) eigenfunctional carries the period Omega^- in its
    scale; the dominant-coordinate-1 normalisation leaves a curve-arbitrary
    rational there, while the primitive-integer normalisation makes the
    odd-character normaliser universal (C = kron(-1, N)/2 on the anchors).
    """
    L = 1
    for x in v:
        L = L * x.denominator // gcd(L, x.denominator)
    w = [int(x * L) for x in v]
    g = 0
    for x in w:
        g = gcd(g, x)
    if g:
        w = [x // g for x in w]
    return [Fraction(x) for x in w]


def raw_twisted_sum(E: EllipticCurve, d: int, N: int, M: ManinSpace = None,
                    lams: dict = None) -> Fraction:
    """The UNNORMALISED exact symbol sum for chi_d (gcd(d, N) = 1):

        S(d) = sum_{b mod m, gcd(b,m)=1} chi_d(b) * lambda^eps({0, b/m}),

    m = |d|, eps = chi_d(-1); S(1) = -lambda^+({0, oo}).  Exact, but in
    the arbitrary scale of the eigenvector (dominant coordinate 1 / primitive
    integers).  Use `algebraic_twisted_l` for a period-normalised value.
    """
    if M is None:
        M = ManinSpace(N)
    if lams is None:
        lams = eigenfunctional(M, E)
    if d == 1:
        inf_sym = M.index[(0, 1)]
        return -sum(l * v for l, v in zip(lams[1], M.symcoord[inf_sym]))
    m = abs(d)
    lam = lams[kron(d, -1)]
    tot = Fraction(0)
    for b in range(1, m + 1):
        if gcd(b, m) != 1:
            continue
        vec = M.coords(M.path_combo(b, m))
        tot += kron(d, b) * sum(l * v for l, v in zip(lam, vec))
    return tot


def omega_pm(E: EllipticCurve, prec: int = 50):
    """(Omega^+, Omega^-): the least positive real period and the least
    positive purely imaginary period (as a positive real) of E's lattice.
    Rectangular lattice (disc > 0): Omega^- = Im(w2); otherwise w2 has real
    part w1/2 and the least imaginary period is 2*w2 - w1, so 2*Im(w2)."""
    from bsdlab.periods import period_lattice
    w1, w2 = period_lattice(E, prec)
    om_minus = w2.imag if E.discriminant > 0 else 2 * w2.imag
    return w1, om_minus


def _identify(x, max_den: int = 1000, tol_exp: int = -12) -> Fraction:
    """A small-height rational equal to the mpf x to 10^tol_exp, else raise."""
    from mpmath import mpf
    q = Fraction(str(x)).limit_denominator(max_den)
    if abs(x - mpf(q.numerator) / q.denominator) > mpf(10) ** tol_exp:
        raise ArithmeticError("no rational of denominator <= %d within "
                              "1e%d of %s" % (max_den, tol_exp, x))
    return q


_SCALES: dict = {}


def period_scales(E: EllipticCurve, N: int, M: ManinSpace = None,
                  lams: dict = None, prec: int = 50, dmax: int = 400) -> dict:
    """{+1: s^+, -1: s^-} with  L(E,chi_d,1)*sqrt|d|/Omega^eps = s^eps*S(d).

    WHY: the eigenvector scale is arbitrary per curve and per sign (the
    2026-09-13 out-of-sample check found the anchor-pinned C wrong by
    factors 2 and 4 on 10/60 unseen twists).  By the theory of modular
    symbols the period-normalised values are rational, so ONE nonzero
    twist per sign fixes the scale: the smallest fundamental d of that
    parity, coprime to N, with S(d) != 0, is evaluated numerically at
    `prec` digits and the ratio identified as a small rational (residual
    < 1e-12, denominator <= 1000 enforced; the L-series delivers ~17 digits).  Every other d is then exact; scratch/oos_modsym.py
    is the out-of-sample gate that this calibration generalises.
    """
    key = (tuple(E.ainvs), N)
    if key in _SCALES:
        return _SCALES[key]
    from mpmath import mp, sqrt
    from bsdlab import twists
    from bsdlab.lseries import l_value
    from bsdlab.reduction import conductor
    if M is None:
        M = ManinSpace(N)
    if lams is None:
        lams = eigenfunctional(M, E)
    omegas = dict(zip((1, -1), omega_pm(E, prec)))
    cands = [1] + [d for d in twists.fundamental_discriminants(dmax)
                   if d != 1 and gcd(d, N) == 1]
    out = {}
    with mp.workdps(prec):
        for eps in (1, -1):
            for d in cands:
                if kron(d, -1) != eps:
                    continue
                S = raw_twisted_sum(E, d, N, M, lams)
                if S == 0:
                    continue
                Ed = twists.twist(E, d).minimal_model() if d != 1 else E
                L = l_value(Ed, conductor(Ed), prec=prec)
                A = L * sqrt(abs(d)) / omegas[eps]
                out[eps] = _identify(A / (mp.mpf(S.numerator) / S.denominator))
                break
            else:
                raise ArithmeticError("no nonzero chi_d, eps=%d, |d|<=%d"
                                      % (eps, dmax))
    _SCALES[key] = out
    return out


def algebraic_twisted_l(E: EllipticCurve, d: int, N: int,
                        M: ManinSpace = None, lams: dict = None) -> Fraction:
    """EXACT L(E,chi_d,1)*sqrt|d| / Omega^{chi_d(-1)}_E  (gcd(d, N) = 1).

    This is the normalisation the Mazur-Tate-Teitelbaum measure needs
    (periods of E itself, not of the twist)."""
    if M is None:
        M = ManinSpace(N)
    if lams is None:
        lams = eigenfunctional(M, E)
    s = period_scales(E, N, M, lams)
    return s[kron(d, -1)] * raw_twisted_sum(E, d, N, M, lams)


def twisted_l_ratio(E: EllipticCurve, d: int, N: int, M: ManinSpace = None,
                    lams: dict = None) -> Fraction:
    """L(E_d, 1)/Omega(E_d) as a Fraction (E_d the minimal twist, gcd(d,N)=1).

    = algebraic_twisted_l * Omega^eps_E / (sqrt|d| * Omega(E_d)), where the
    last factor is a lattice index (a small rational: 1/2, 1, 2, ...)
    identified numerically at 50 digits.  The symbol sum is exact; the two
    normalising rationals are identified-and-checked, not guessed.
    """
    from mpmath import mp, sqrt
    from bsdlab import twists
    from bsdlab.periods import real_period
    alg = algebraic_twisted_l(E, d, N, M, lams)
    if alg == 0:
        return Fraction(0)
    eps = kron(d, -1)
    Ed = twists.twist(E, d).minimal_model() if d != 1 else E
    with mp.workdps(50):
        om = omega_pm(E, 50)[0 if eps == 1 else 1]
        index = _identify(om / (sqrt(abs(d)) * real_period(Ed, 50)))
    return alg * index
