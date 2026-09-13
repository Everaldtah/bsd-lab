"""L3 — the regulator oracle: score candidate Mordell-Weil generators
against the BSD prediction.

For a curve E/Q with analytic rank r and candidate generators P_1..P_r,

    Reg_meas = det(<P_i, P_j>)   (canonical height pairing)

while BSD predicts (with #Sha = 1, full generators)

    Reg_BSD = L^(r)(1)/r! * |E_tors|^2 / (prod c_p * Omega).

If the P_i are independent and generate a finite-index subgroup, the ratio
Reg_meas / Reg_BSD equals #Sha_an * index^2 -- in particular an integer.
The oracle computes BOTH sides as W1b-grade intervals and reports the
ratio together with a verdict:

    PASS         ratio interval contains an integer (and r matches);
    FAIL         ratio interval excludes every integer, or r mismatched;
    DEPENDENT    the height-pairing determinant vanishes;
    INSUFFICIENT fewer generators than the analytic rank.

This is the falsifier every Route-B construction is scored against.
"""

from dataclasses import dataclass, field
from typing import Optional

from mpmath import mp, mpf

from .certify import (leading_coefficient_interval, real_period_interval,
                      regulator_interval, iv_div, iv_mul)

VERDICTS = ('PASS', 'FAIL', 'DEPENDENT', 'INSUFFICIENT')


@dataclass
class OracleResult:
    label: Optional[str] = None
    analytic_rank: Optional[int] = None
    n_generators: int = 0
    reg_measured: Optional[tuple] = None    # (lo, hi)
    reg_bsd: Optional[tuple] = None         # (lo, hi)
    ratio: Optional[tuple] = None           # (lo, hi) for reg_measured/reg_bsd
    nearest_integer: Optional[int] = None
    distance: Optional[mpf] = None          # |ratio - nearest| in interval radii
    verdict: str = 'FAIL'
    notes: list = field(default_factory=list)

    def note(self, msg: str) -> None:
        self.notes.append(msg)


def _check_points(E, points) -> list:
    """Drop None/O, verify each remaining point lies on E (else raise)."""
    from .curve import O, INFINITY

    kept = []
    for P in points:
        if P is None or P is O or isinstance(P, INFINITY):
            continue
        if not E.is_on_curve(P):
            raise ValueError('point %r is not on the curve %r' % (P, E))
        kept.append(P)
    return kept


def test_generators(E, generators, label=None, prec=50, prec2=80
                    ) -> OracleResult:
    """Score `generators` against the BSD-predicted regulator."""
    from . import reduction

    gens = _check_points(E, generators)
    res = OracleResult(label=label, n_generators=len(gens))

    N = reduction.conductor(E)
    r_an, l_iv, _tail, _rnd = leading_coefficient_interval(
        E, N, prec=prec, prec2=prec2)
    o_iv, _agm, _ornd = real_period_interval(E, prec=prec, prec2=prec2)
    res.analytic_rank = r_an
    with mp.workdps(prec2):
        tors = E.torsion_order()
        tam = reduction.tamagawa_product(E)
        # Reg_BSD = L^(r)(1)/r! * tors^2 / (tam * Omega)
        num = iv_mul(l_iv, (mpf(tors ** 2),) * 2)
        den = iv_mul(o_iv, (mpf(tam),) * 2)
        res.reg_bsd = iv_div(num, den)

    m_iv, m_err = regulator_interval(E, gens, prec=prec, prec2=prec2)
    res.reg_measured = m_iv
    with mp.workdps(prec2):
        res.ratio = iv_div(m_iv, res.reg_bsd)
    if m_iv[0] <= 0 <= m_iv[1] and gens:
        res.verdict = 'DEPENDENT'
        res.note('height-pairing determinant compatible with 0')
        return res
    return _verdict(res, prec2)


def _verdict(res: OracleResult, prec: int) -> OracleResult:
    """PASS iff the ratio interval traps an integer and ranks agree."""
    r, n = res.analytic_rank, res.n_generators
    if n < r:
        res.verdict = 'INSUFFICIENT'
        res.note('%d generators for analytic rank %d' % (n, r))
        return res
    if n > r:
        res.verdict = 'FAIL'
        res.note('RANK EXCESS: %d independent generators, analytic rank %d'
                 % (n, r))
        return res
    with mp.workdps(prec):
        lo, hi = res.ratio
        k = int(mp.nint(mp.mpf((lo + hi) / 2)))
        res.nearest_integer = k
        # distance of k OUTSIDE the interval; 0 means k is trapped inside
        res.distance = max(mpf(0), lo - k, mp.mpf(k) - hi)
        if res.distance == 0:
            res.verdict = 'PASS'
            res.note('ratio = %d x (index^2 * Sha) consistent' % k
                     if k > 1 else 'ratio consistent with 1')
        elif res.distance > 0:
            res.verdict = 'FAIL'
            res.note('ratio interval [%s, %s] excludes every integer'
                     % (mp.nstr(lo, 8), mp.nstr(hi, 8)))
    return res


def main(argv=None) -> int:
    """CLI: python -m bsdlab.oracle --ainvs 0 -1 1 -10 -20 P1 [P2..]

    Each point is 'x,y' with x, y rationals spelled 'a/b' or 'a'.
    """
    import argparse
    from fractions import Fraction

    from .curve import EllipticCurve

    def frac(s):
        n, _, d = s.partition('/')
        return Fraction(int(n), int(d or 1))

    ap = argparse.ArgumentParser(description='BSD regulator oracle')
    ap.add_argument('--ainvs', type=int, nargs=5, required=True)
    ap.add_argument('--label', default=None)
    ap.add_argument('--prec', type=int, nargs=2, default=[50, 80],
                    metavar=('PREC', 'PREC2'))
    ap.add_argument('points', nargs='*', default=[])
    args = ap.parse_args(argv)

    E = EllipticCurve.from_list(args.ainvs)
    pts = [tuple(map(frac, p.split(','))) for p in args.points]
    res = test_generators(E, pts, label=args.label, prec=args.prec[0],
                          prec2=args.prec[1])
    with mp.workdps(args.prec[1]):
        print(f'curve        {args.ainvs} ({args.label or "?"})')
        print(f'analytic rank {res.analytic_rank}, generators {res.n_generators}')
        print(f'Reg_meas     [{mp.nstr(res.reg_measured[0], 17)},'
              f' {mp.nstr(res.reg_measured[1], 17)}]')
        print(f'Reg_BSD      [{mp.nstr(res.reg_bsd[0], 17)},'
              f' {mp.nstr(res.reg_bsd[1], 17)}]')
        mid = (res.ratio[0] + res.ratio[1]) / 2
        half = (res.ratio[1] - res.ratio[0]) / 2
        print(f'ratio        {mp.nstr(mid, 12)} +/- {mp.nstr(half, 3)}'
              f'  -> nearest int {res.nearest_integer}')
        print(f'VERDICT      {res.verdict}')
        for msg in res.notes:
            print(f'  note: {msg}')
    return 0 if res.verdict == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
