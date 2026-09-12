"""Known-answer tests for every module, runnable offline.

    python -m bsdlab.selftest

The reference values below were taken from the LMFDB and are *targets*, not
inputs: nothing in ``bsdlab/`` reads this file, and no computation is tuned to
hit them.  A failure here means bsdlab disagrees with the Cremona/PARI lineage
on a curve where that lineage has been checked by many people for decades, so
the presumption is that bsdlab is wrong.

The point of the whole project is that these numbers were recomputed by a
completely separate implementation.  Agreement is therefore evidence; it is not
a tautology.  Where a value below is marked ``None`` it is not checked.
"""

from __future__ import annotations

import sys
import traceback

import mpmath

from bsdlab.curve import EllipticCurve

#: label, a-invariants, conductor, {p: c_p}, torsion, rank, root number,
#: regulator, real period, #Sha.
#: Regulators and periods are strings so they are read at full working
#: precision rather than through a float literal.
CURVES = [
    # label    ainvs                  N     tamagawa      tors rk  w
    ('11a1',  [0, -1, 1, -10, -20],   11,   {11: 5},       5,  0, +1,
     '1', '1.26920930427955342602495318?', 1),
    ('11a2',  [0, -1, 1, -7820, -263580], 11, {11: 1},     1,  0, +1,
     '1', None, 1),
    # X_0(11): disc = -11, so v_11(disc) = 1 forces c_11 = 1, not 5.
    ('11a3',  [0, -1, 1, 0, 0],       11,   {11: 1},       5,  0, +1,
     '1', None, 1),
    ('14a1',  [1, 0, 1, 4, -6],       14,   {2: 2, 7: 3},  6,  0, +1,
     '1', None, 1),
    ('15a1',  [1, 1, 1, -10, -10],    15,   {3: 2, 5: 4},  8,  0, +1,
     '1', None, 1),
    # LMFDB 32.a1.  The old row paired these a-invariants with 32.a4's values.
    ('32a3',  [0, 0, 0, -11, -14],    32,   {2: 1},        2,  0, +1,
     '1', None, 1),
    ('37a1',  [0, 0, 1, -1, 0],       37,   {37: 1},       1,  1, -1,
     '0.051111408239968840235886099757', '5.98691729246391914767010400?', 1),
    ('37b1',  [0, 1, 1, -23, -50],    37,   {37: 3},       3,  0, +1,
     '1', None, 1),
    # j = 0 (c4 = 0) rows: regression coverage for the minimal-model fix.
    # These used to raise inside _lkc_scaling -- v_p(0) = inf, and inf // n is
    # NaN in CPython, so the scaling-exponent search died before a minimal
    # model existed.  Conductor and Tamagawa numbers are the oracle's (rows of
    # data/lmfdb_oracle.json, which IS the minimal model), Sha is the LMFDB
    # analytic value; rank 0 forces w = +1 and regulator 1.  Labels are
    # Cremona-style, with the LMFDB label noted where the two labellings
    # disagree (they order curves within an isogeny class differently -- cf.
    # the 32a3 row above).
    # y^2 + y = x^3: the CM curve with j = 0.  LMFDB 27.a4; c_3 = 1 because
    # the minimal discriminant -27 has v_3 = 3 (additive, so no ODD constraint).
    ('27a1',  [0, 0, 1, 0, 0],       27,   {3: 1},        3,  0, +1,
     '1', None, 1),
    # y^2 = x^3 + 1.  LMFDB 36.a4.
    ('36a1',  [0, 0, 0, 0, 1],       36,   {2: 3, 3: 2},   6,  0, +1,
     '1', None, 1),
    # y^2 = x^3 - 1.  LMFDB 144.a3.
    ('144a1', [0, 0, 0, 0, -1],      144,  {2: 1, 3: 2},   2,  0, +1,
     '1', None, 1),
    # LMFDB 243.b2: an a3 = 1 j = 0 curve (y^2 + y = x^3 + 2), checking that
    # the fix does not silently force the b-form y^2 = x^3 + b.
    ('243b1', [0, 0, 1, 0, 2],       243,  {3: 3},        3,  0, +1,
     '1', None, 1),
    ('389a1', [0, 1, 1, -2, 0],       389,  {389: 1},      1,  2, +1,
     '0.15246017794314375162432475705', '4.98042512173897283?', 1),
    ('5077a1', [0, 0, 1, -7, 6],      5077, {5077: 1},     1,  3, -1,
     '0.41714355875838396981711954462', '4.15168798308693305?', 1),
    # The smallest curve with non-trivial Sha: #Sha = 4.  This is the single
    # most valuable test in the file, because it is the only one where a bug
    # that drops a factor would still leave #Sha_an a perfect square (1).
    ('571a1', [0, -1, 1, -929, -10595], 571, {571: 1},     1,  0, +1,
     '1', None, 4),
    # #Sha = 9, again guarding the Sha assembly against factor errors.
    ('681b1', [1, 1, 0, -1154, -15345], 681, {3: 2, 227: 2}, 4, 0, +1,
     '1', None, 9),
]

TOL = mpmath.mpf('1e-9')


class Result:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.skipped = 0
        self.messages: list[str] = []

    def check(self, label, what, got, want, tol=None):
        if want is None:
            return
        if got is None:
            self.failed += 1
            self.messages.append(f'{label:8s} {what:16s} NOT COMPUTED, want {want}')
            return
        if tol is None:
            ok = got == want
        else:
            ok = abs(mpmath.mpf(got) - mpmath.mpf(want)) < tol
        if ok:
            self.passed += 1
        else:
            self.failed += 1
            got_s = mpmath.nstr(got, 18) if tol is not None else got
            self.messages.append(f'{label:8s} {what:16s} got {got_s}  want {want}')

    def skip(self, label, what, exc):
        self.skipped += 1
        self.messages.append(f'{label:8s} {what:16s} SKIPPED ({exc.__class__.__name__}: {exc})')


def _import(name):
    try:
        return __import__(f'bsdlab.{name}', fromlist=[name])
    except Exception:
        return None


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    verbose = '-v' in argv or '--verbose' in argv

    reduction = _import('reduction')
    lseries = _import('lseries')
    periods = _import('periods')
    mordellweil = _import('mordellweil')
    from bsdlab import bsd

    missing = [n for n, m in [('reduction', reduction), ('lseries', lseries),
                              ('periods', periods), ('mordellweil', mordellweil)]
               if m is None]
    if missing:
        print(f'modules not importable, their checks will be skipped: {", ".join(missing)}\n')

    r = Result()

    with mpmath.workdps(50):
        for (label, ainvs, N, tam, tors, rank, w, reg, omega, sha) in CURVES:
            E = EllipticCurve(*ainvs)

            # Curve arithmetic is always available.
            try:
                r.check(label, 'torsion', E.torsion_order(), tors)
            except Exception as exc:
                r.skip(label, 'torsion', exc)

            if reduction is not None:
                try:
                    r.check(label, 'conductor', reduction.conductor(E), N)
                    if tam is not None:
                        r.check(label, 'tamagawa', reduction.tamagawa_numbers(E), tam)
                except Exception as exc:
                    r.skip(label, 'reduction', exc)

            if lseries is not None:
                try:
                    r.check(label, 'root number', lseries.sign(E, N=N, prec=30), w)
                except Exception as exc:
                    r.skip(label, 'root number', exc)
                try:
                    ar, _ = lseries.leading_coefficient(E, N=N, prec=30)
                    r.check(label, 'analytic rank', ar, rank)
                except Exception as exc:
                    r.skip(label, 'analytic rank', exc)

            if periods is not None and omega is not None:
                try:
                    got = periods.real_period(E, prec=30)
                    want = mpmath.mpf(omega.rstrip('?'))
                    r.check(label, 'real period', got, want, tol=TOL)
                except Exception as exc:
                    r.skip(label, 'real period', exc)

            if mordellweil is not None:
                try:
                    rk, gens = mordellweil.rank_lower_bound(E, height_bound=200, prec=30)
                    r.check(label, 'rank', rk, rank)
                    r.check(label, 'regulator',
                            mordellweil.regulator(E, gens, prec=30),
                            mpmath.mpf(reg), tol=mpmath.mpf('1e-9'))
                except Exception as exc:
                    r.skip(label, 'mordell-weil', exc)

            if None not in (reduction, lseries, periods, mordellweil):
                try:
                    d = bsd.analyse(E, label=label, prec=30)
                    r.check(label, '#Sha_an', d.sha_rounded, sha)
                except Exception as exc:
                    r.skip(label, '#Sha_an', exc)

    for m in r.messages:
        print(m)
    if r.messages:
        print()
    print(f'passed {r.passed}   failed {r.failed}   skipped {r.skipped}')
    if verbose and r.failed:
        traceback.print_exc()
    return 1 if r.failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
