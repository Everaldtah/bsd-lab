"""Interval certificates for the analytic order of Sha.

`bsd.sha_analytic` evaluates

    #Sha_an = ( L^(r)(1)/r! ) * |E_tors|^2 / ( Omega * Reg * prod c_p )

in floating point and then *trusts* that the result is near an integer
(SHA_INTEGRALITY_TOLERANCE).  This module replaces that trust with a
rigorous-looking error budget assembled from three a-posteriori bounds:

  L-series   the truncation tail of the smoothed central-point sums,
             bounded via Hasse |a_n| <= 2 sqrt(n) and the closed form
             G_r(x) <= e^{-x} sum_{k<=r} x^k / (k! x^{r+1}), summed as a
             geometric series (terms decay like exp(-2 pi / sqrt(N)));
             the same bound replaces lseries' 10^(-prec/2) zero-test in
             the rank scan (`analytic_rank_bounded`), because for large N
             the tail is precision-independent and sinks below that
             threshold (571b1 at prec 80: tail 3e-39 vs threshold 1e-40,
             so the plain scan misreports rank 0 for a rank-2 curve);
  Omega      the AGM last-iterate gap: the AGM limit M lies between the
             two final iterates (a_n decreases to M, b_n increases to M),
             so the relative error is at most the final |a-b|/|a|, which
             we track through our own replay of the AGM;
  rounding   two-precision differencing: every ingredient is computed at
             dps=prec and dps=prec2>prec; |v_prec - v_prec2| * SAFETY is
             taken as the rounding component.  This is an engineering
             proxy, not a theorem -- see the honesty note in the budget.

Exact inputs (torsion, Tamagawa product, number of real components) carry
no interval.  The result is a certificate: #Sha_an in [lo, hi], and if
[lo, hi] sits strictly inside (n - 1/2, n + 1/2) for an integer n >= 1,
then #Sha = n is certified *conditional on strong BSD and on the
un-certified rank* (rank <= 1 is Gross-Zagier + Kolyvagin; rank >= 2
remains a lower bound only, as recorded by bsd.analyse).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import mpmath
from mpmath import mp, mpf

from bsdlab.curve import EllipticCurve

#: two-precision differencing safety factor (lab methodology, PLAN.md L1).
SAFETY = mpf(10)
#: geometric-ratio cap for the L tail; needs q < 1, checked at runtime.


def iv_mul(u, v):
    """Endpoint product of two positive intervals (lo, hi)."""
    return (u[0] * v[0], u[1] * v[1])


def iv_div(u, v):
    """Endpoint quotient of two positive intervals (lo, hi)."""
    return (u[0] / v[1], u[1] / v[0])


def iv_from_error(value, err):
    """The interval [value - err, value + err]; err must be >= 0."""
    if err < 0:
        raise ValueError('error bound must be nonnegative')
    return (value - err, value + err)


def iv_narrow(*ivals):
    """Intersection of positive intervals; raises if empty."""
    lo = max(i[0] for i in ivals)
    hi = min(i[1] for i in ivals)
    if not lo <= hi:
        raise ValueError('interval intersection is empty')
    return (lo, hi)


@dataclass
class ErrorBudget:
    """Every error component entering the #Sha interval, in absolute units."""
    l_tail: mpf = mpf(0)
    l_round: mpf = mpf(0)
    omega_agm: mpf = mpf(0)
    omega_round: mpf = mpf(0)
    regulator: mpf = mpf(0)

    def total_relative(self, sha_mid):
        """(sum of propagated absolute errors) / |sha_mid|."""
        if sha_mid == 0:
            return mpf('inf')
        tot = (self.l_tail + self.l_round + self.omega_agm
               + self.omega_round + self.regulator)
        return tot / abs(sha_mid)


@dataclass
class ShaCertificate:
    """A machine-checkable statement  #Sha_an in [lo, hi]."""
    label: Optional[str] = None
    rank: Optional[int] = None
    rank_certified: bool = False
    lo: Optional[mpf] = None
    hi: Optional[mpf] = None
    sha_certified: Optional[int] = None
    budget: ErrorBudget = field(default_factory=ErrorBudget)
    notes: list = field(default_factory=list)


def _Ghat(k: int, x) -> mpf:
    """The dominating envelope of G_k(x), exact and explicit.

    G_k(x) = (1/k!) int_1^inf (log u)^k e^{-x u} du  <=  (1/k!) int_1^inf
    u^k e^{-x u} du  (log u <= u on [1, inf))  =  Gamma(k+1, x)/(k! x^{k+1})
    =  e^{-x} sum_{i=0}^{k} x^{i-k-1} / i!   (closed form for integer k).
    A sum of decreasing powers of x, hence decreasing in x; equality for
    k = 0.
    """
    return mp.exp(-x) * sum(x ** (i - k - 1) / mp.factorial(i)
                            for i in range(k + 1))


def _lambda_tail(N: int, k: int, n_terms: int) -> mpf:
    """Bound on the truncation error of Lambda_k = (1 + w(-1)^k) sum a_n G_k(x_n).

    Terms beyond n_terms carry |a_n| <= 2 sqrt(n) (Hasse), so the tail of
    Lambda_k is at most 2 * sum_{n>n_terms} 2 sqrt(n) G_k(x_n), bounded by
    4 times a geometric series: term(n) = sqrt(n) Ghat_k(x_n) satisfies
    term(n+1)/term(n) <= sqrt((n+1)/n) * exp(-2 pi/sqrt(N)) =: q(n) <= q,
    and q < 1 for every conductor (checked).  The (1 + w(-1)^k) factor is
    at most 2 and already included.
    """
    q = mpf(n_terms + 2) / mpf(n_terms + 1)
    q = mp.sqrt(q) * mp.exp(-2 * mp.pi / mp.sqrt(N))
    if q >= 1:
        raise ArithmeticError(
            'geometric ratio %s >= 1 for N=%d; the tail is not controlled'
            % (mp.nstr(q, 6), N))
    first = mp.sqrt(n_terms + 1) * _Ghat(k, 2 * mp.pi * (n_terms + 1)
                                         / mp.sqrt(N))
    return 4 * first * q / (1 - q)


def _sum_bound(N: int, k: int) -> mpf:
    """4 * sum_{n>=1} sqrt(n) Ghat_k(x_n): a bound on |Lambda_k|'s full
    (untruncated) sum of |a_n G_k| terms, used to size the rounding noise
    of adding them up (each of the n_terms additions rounds once, by at
    most ulp ~ 10^{1-prec} of the partial sum).

    The step-ratio bound term(n+1)/term(n) <= sqrt((n+1)/n) * exp(-2 pi/
    sqrt(N)) is decreasing in n and equals 1 at n ~= sqrt(N)/(4 pi), so
    for N > ~329 the series is NOT decreasing from n = 1: sum the terms
    explicitly up to the peak, then close geometrically with the (now < 1
    and decreasing) ratio bound."""
    step = mpf(2) * mp.pi / mp.sqrt(N)
    head = mpf(0)
    n = 1
    while mp.sqrt(mpf(n + 1) / n) * mp.exp(-step) >= 1:
        head += mp.sqrt(mpf(n)) * _Ghat(k, step * n)
        n += 1
    q = mp.sqrt(mpf(n + 1) / n) * mp.exp(-step)
    peak = mp.sqrt(mpf(n)) * _Ghat(k, step * n)
    # sum <= head + peak + peak*q + peak*q^2 + ... = head + peak/(1 - q)
    return 4 * (head + peak / (1 - q))


def _l_truncation_error(N: int, w: int, alpha, h, r: int,
                        n_terms: int) -> mpf:
    """Rigorous bound on the truncation error of L_r = (sum_j Lam_{r-j}
    h_j)/alpha, each Lambda_k tail via `_lambda_tail` (zero when the
    (1 + w(-1)^k) factor kills Lambda_k structurally)."""
    tot = mpf(0)
    for j in range(r + 1):
        k = r - j
        if (1 + w * (-1) ** k) == 0:
            continue
        tot += abs(h[j]) * _lambda_tail(N, k, n_terms)
    return tot / alpha


def analytic_rank_bounded(E: EllipticCurve, N: int, prec: int = 80,
                          max_rank: int = 4) -> tuple:
    """(r, L_r, floor_r): the first parity-allowed r with |L_r| above its
    own noise floor = analytic truncation tail + estimated rounding.

    This replaces lseries' 10^(-prec/2) zero-test, which silently assumes
    the truncation tail sits below that threshold.  It does not for large
    N at high precision (571b1: tail e^{-89} ~ 3e-39 vs threshold 1e-40
    at prec 80, so plain analytic_rank reports 0 instead of 2); and
    dually, at LOW precision the rounding of the cancelled sum dominates
    the tail (389a1 at prec 30: cancellation noise ~1e-28 vs tail ~1e-41),
    so the tail alone would misread rank 0 there.  L_r is declared
    nonzero exactly when |L_r| exceeds both: the tail rigorously, the
    rounding by a generous estimate (n_terms additions, each ulp of the
    partial sums, times SAFETY).
    """
    from bsdlab import lseries

    with mp.workdps(prec):
        w, alpha, L = lseries._completed_moments(E, max_rank, N=N, prec=prec)
        h = lseries._h_coeffs(alpha, max_rank)
        n_terms = lseries._default_n_terms(N)
        ulp = mpf(10) ** (1 - prec) * SAFETY * 4 * n_terms
        for r in range(0 if w == 1 else 1, max_rank + 1, 2):
            tail = _l_truncation_error(N, w, alpha, h, r, n_terms)
            noise = sum(abs(h[j]) * _sum_bound(N, r - j)
                        for j in range(r + 1)) * ulp / alpha
            if abs(L[r]) > tail + noise:
                return r, L[r], tail + noise
        raise ArithmeticError(
            'no |L_r| above its noise floor for r <= %d (N=%d, prec=%d): '
            'raise max_rank or prec' % (max_rank, N, prec))


def leading_coefficient_interval(E: EllipticCurve, N: int, prec: int = 50,
                                 prec2: int = 80) -> tuple:
    """(r, (lo, hi), tail, round) for L^(r)(1)/r! at the analytic rank r.

    r comes from the tail-aware `analytic_rank_bounded` at both
    precisions; a disagreement is fatal (raises) rather than averaged
    away.  The truncation tail is the analytic `_lambda_tail` propagated
    through the convolution L_r = (sum_j Lambda_{r-j} h_j)/alpha; the
    rounding component is the two-precision diff times SAFETY, floored at
    a few ulp of the prec2 value (SAFETY cannot rescue a computation that
    agrees only because both precisions hit the same rounding pattern).
    """
    from bsdlab import lseries

    r_hi, lead_hi, tail = analytic_rank_bounded(E, N, prec=prec2)
    with mp.workdps(prec):
        r_lo, lead_lo, _ = analytic_rank_bounded(E, N, prec=prec)
    if r_hi != r_lo:
        raise ArithmeticError(
            'analytic rank unstable across prec %d/%d (%d vs %d) for %r'
            % (prec, prec2, r_lo, r_hi, E))

    with mp.workdps(prec2):
        ulp_floor = abs(lead_hi) * mpf(10) ** (2 - prec2)
        rnd = abs(lead_hi - lead_lo) * SAFETY + ulp_floor
        return r_hi, iv_from_error(lead_hi, tail + rnd), tail, rnd


def _agm_gap(a, b, tol):
    """Run the periods._agm loop, returning (M, final |a-b|/min(|a|,|b|)).

    The AGM limit M of the initial pair lies strictly between the two
    iterates at every stage (a_n decreases to M, b_n increases to M), so
    the final gap IS a rigorous a-posteriori bound on the truncation
    error of the iteration itself -- before input rounding, which the
    two-precision diff of the caller covers.
    """
    while abs(a - b) > tol * abs(a):
        a, b = (a + b) / 2, mp.sqrt(a * b)
    return a, abs(a - b) / min(abs(a), abs(b))


def real_period_interval(E: EllipticCurve, prec: int = 50,
                         prec2: int = 80) -> tuple:
    """(lo, hi) for Omega with (agm, round) components.

    Omega = 2 * w_re * ncomponents with w_re = pi / (2 AGM(a, b)), so the
    only iteration whose truncation matters is AGM(a, b); we replay it at
    prec2 recording the final gap.  polyroots' own error on the branch
    points is not bounded analytically here -- the two-precision diff
    times SAFETY is the lab's mandated proxy for it.
    """
    from bsdlab import periods

    # The calls sit in explicit contexts: periods.real_period does its last
    # multiplication at the *ambient* precision, so an unguarded call would
    # silently truncate Omega to the default 15 digits.
    with mp.workdps(prec2 + 10):
        omega_hi = periods.real_period(E, prec2)
    with mp.workdps(prec + 10):
        omega_lo = periods.real_period(E, prec)
    with mp.workdps(prec2 + 10):
        e1, e2, e3 = periods._cubic_roots(E)
        a = mp.sqrt(e1 - e3)
        b = mp.sqrt(e1 - e2) if E.discriminant > 0 else a.conjugate()
        tol = mpf(10) ** (-(mp.dps - 5))
        _, gap = _agm_gap(a, b, tol)
        agm = abs(omega_hi) * (gap + mpf(10) ** (2 - prec2))
        rnd = abs(omega_hi - omega_lo) * SAFETY
        return iv_from_error(omega_hi, agm + rnd), agm, rnd


def regulator_interval(E: EllipticCurve, generators, prec: int = 50,
                       prec2: int = 80) -> tuple:
    """(lo, hi) for the regulator, with its single error component.

    Silverman's local-height series converges to full working precision
    in ~0.51*dps terms (mordellweil.py), so there is no analytic tail to
    bound: the two-precision diff times SAFETY, floored at a few ulp,
    carries both the series truncation and the rounding.  Rank 0 has
    regulator exactly 1 (empty generator list) and zero error.
    """
    from bsdlab import mordellweil

    gens = [P for P in generators
            if P is not None and not _is_zero_point(P)]
    if not gens:
        return (mpf(1), mpf(1)), mpf(0)
    with mp.workdps(prec2):
        r_hi = mordellweil.regulator(E, gens, prec2)
    with mp.workdps(prec):
        r_lo = mordellweil.regulator(E, gens, prec)
    with mp.workdps(prec2):
        err = (abs(r_hi - r_lo) * SAFETY
               + abs(r_hi) * mpf(10) ** (2 - prec2))
        return iv_from_error(r_hi, err), err


def _is_zero_point(P) -> bool:
    """The point at infinity, however spelled (O or INFINITY)."""
    from bsdlab.curve import O, INFINITY
    return P is O or isinstance(P, INFINITY)


def certify_sha_interval(sha_iv) -> Optional[int]:
    """The unique integer n with [lo, hi] inside (n - 1/2, n + 1/2), or None.

    Strict containment: the endpoints may not sit exactly on a half-integer
    boundary, where the nearest integer would be ambiguous.  mpf
    comparisons are exact, so the strictness is decidable, not fuzzy.
    """
    lo, hi = sha_iv
    if not (mpf(0) < lo <= hi):
        return None
    n = int(mp.floor(lo + mpf('0.5')))
    if n >= 1 and lo > n - mpf('0.5') and hi < n + mpf('0.5'):
        return n
    return None


def certify_curve(E: EllipticCurve, label: Optional[str] = None,
                  prec: int = 50, prec2: int = 80,
                  height_bound: int = 200) -> ShaCertificate:
    """Assemble the full error budget and certify #Sha_an for one curve.

    The pipeline of `bsd.analyse` runs once at prec2 (single source of
    truth for conductor, Tamagawa numbers, torsion, generators, rank);
    each floating ingredient is then re-derived with its interval, and
    the exact integers enter as degenerate intervals.  The certificate
    states `#Sha_an in [lo, hi]` and, when that pins an integer, the
    integer itself -- always conditional on strong BSD, and on the
    algebraic rank being the analytic one (r <= 1: Gross--Zagier +
    Kolyvagin; r >= 2: only a search lower bound, noted).
    """
    from bsdlab import bsd

    E = E.minimal_model()
    cert = ShaCertificate(label=label)
    # analyse runs at the LOWER precision: its 10^(-prec/2) zero-test is
    # calibrated for the truncation tail there, while the certificate's
    # own rank comes from the tail-aware scan below and is cross-checked.
    d = bsd.analyse(E, label=label, prec=prec, height_bound=height_bound)
    cert.notes.extend(d.notes)
    cert.rank = d.analytic_rank
    cert.rank_certified = d.rank_certified
    if d.sha_analytic is None or d.conductor is None:
        cert.notes.append('pipeline incomplete: no #Sha interval possible')
        return cert

    try:
        r_cert, (l_lo, l_hi), l_tail, l_round = leading_coefficient_interval(
            E, d.conductor, prec, prec2)
        (o_lo, o_hi), o_agm, o_round = real_period_interval(E, prec, prec2)
        (r_lo, r_hi), r_err = regulator_interval(E, d.generators, prec, prec2)
    except Exception as exc:                        # noqa: BLE001
        cert.notes.append(f'interval computation failed: {exc!r}')
        return cert

    if r_cert != d.analytic_rank:
        cert.notes.append(f'RANK DISAGREEMENT: tail-aware scan {r_cert} vs '
                          f'analyse {d.analytic_rank}; certificate uses the '
                          'tail-aware rank')
        cert.rank = r_cert
        cert.rank_certified = r_cert in (0, 1)

    cert.budget = ErrorBudget(l_tail=l_tail, l_round=l_round,
                              omega_agm=o_agm, omega_round=o_round,
                              regulator=r_err)

    # Endpoint arithmetic and the half-integer test run at prec2: at the
    # ambient default of 15 digits the additions and comparisons would
    # round the interval onto its centre and the certificate would lie.
    with mp.workdps(prec2 + 10):
        num = iv_mul((l_lo, l_hi), (mpf(d.torsion_order ** 2),) * 2)
        den = iv_mul(iv_mul((o_lo, o_hi), (r_lo, r_hi)),
                     (mpf(d.tamagawa_product),) * 2)
        sha_iv = iv_div(num, den)
        cert.lo, cert.hi = sha_iv
        cert.sha_certified = certify_sha_interval(sha_iv)
        if (d.sha_analytic is not None and r_cert == d.analytic_rank
                and not sha_iv[0] <= d.sha_analytic <= sha_iv[1]):
            cert.notes.append('INTERNAL: heuristic #Sha_an outside the '
                              'rigorous interval')
            cert.sha_certified = None

    if r_cert != d.analytic_rank:
        # the containment checks below would compare quantities taken at
        # different ranks and so prove nothing
        return cert
    if not (l_lo <= d.leading_coefficient <= l_hi
            and o_lo <= d.real_period <= o_hi):
        cert.notes.append('INTERNAL: analyse value outside its own interval')
        cert.sha_certified = None
        return cert
    if cert.sha_certified is None and not cert.notes:
        cert.notes.append('interval too wide: #Sha_an not pinned to one '
                          'integer at this precision')
    return cert


def _fmt(x, n=6) -> str:
    return '-' if x is None else mp.nstr(x, n)


def main(argv=None) -> int:
    """Certify the reference set:  python -m bsdlab.certify [--prec P P2].

    Prints one row per curve -- the interval, its half-width relative to
    its centre, the certified integer (or '-'), and whether the rank
    behind it is Gross--Zagier + Kolyvagin certified -- and writes
    data/certificates.json for machine checking.
    """
    import argparse
    import json
    from pathlib import Path

    ap = argparse.ArgumentParser(description='interval certificates for #Sha')
    ap.add_argument('--ref', default='data/lmfdb_reference.json')
    ap.add_argument('--out', default='data/certificates.json')
    ap.add_argument('--prec', type=int, nargs=2, default=[50, 80],
                    metavar=('PREC', 'PREC2'))
    ap.add_argument('--labels', nargs='*', default=[])
    args = ap.parse_args(argv)
    prec, prec2 = args.prec

    rows = json.loads(Path(args.ref).read_text())
    if args.labels:
        rows = [r for r in rows if r['Clabel'] in args.labels]

    results = []
    print(f'{"label":<9}{"rk":>3}{"GZK":>4}{"#Sha_an interval":<44}'
          f'{"rel.err":>10}{"cert":>6}')
    for row in rows:
        from bsdlab.curve import EllipticCurve
        E = EllipticCurve.from_list(row['ainvs'])
        try:
            c = certify_curve(E, label=row['Clabel'], prec=prec, prec2=prec2)
        except Exception as exc:                    # noqa: BLE001
            print(f'{row["Clabel"]:<9}  ERROR {exc!r}')
            results.append({'label': row['Clabel'], 'error': repr(exc)})
            continue
        mid = None
        rel = None
        ival = '-'
        if c.lo is not None:
            with mp.workdps(prec2):
                mid = (c.lo + c.hi) / 2
                half = (c.hi - c.lo) / 2
                rel = half / mid
            ival = f'{mp.nstr(mid, 18)} +/- {mp.nstr(half, 3)}'
        ival = (f'[{_fmt(c.lo, 17)}, {_fmt(c.hi, 17)}]'
                if c.lo is not None else '-')
        print(f'{c.label:<9}{c.rank if c.rank is not None else -1:>3}'
              f'{"y" if c.rank_certified else "n":>4}{ival:<44}'
              f'{_fmt(rel, 3):>10}'
              f'{c.sha_certified if c.sha_certified is not None else "-":>6}')
        results.append({
            'label': c.label, 'rank': c.rank, 'rank_certified':
            c.rank_certified,
            'lo': _fmt(c.lo, 25), 'hi': _fmt(c.hi, 25),
            'sha_certified': c.sha_certified,
            'sha_ref': row.get('sha'),
            'budget': {k: _fmt(v, 6) for k, v in vars(c.budget).items()},
            'notes': c.notes,
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=1))
    n_ok = sum(1 for r in results if r.get('sha_certified') is not None)
    print(f'\n{n_ok}/{len(results)} certificates pin an integer; '
          f'wrote {args.out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
