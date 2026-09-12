"""Assembly of the strong Birch--Swinnerton-Dyer formula.

For an elliptic curve E/Q of analytic rank r, strong BSD predicts

    L^(r)(E, 1) / r!  =  Omega * Reg * prod_p c_p * #Sha / |E(Q)_tors|^2

Everything on the right except #Sha is computable.  We therefore solve for #Sha
and call the result the *analytic order of Sha*:

    #Sha_an  =  ( L^(r)(1)/r! ) * |E(Q)_tors|^2 / ( Omega * Reg * prod_p c_p )

This is a prediction, not a computation of the Tate--Shafarevich group.  It is
the number that Sha *would* have to have if strong BSD holds and if the other
five quantities are correct.  Two independent facts make it a useful error
detector rather than a tautology:

  * If Sha is finite its order is a perfect square (Cassels: the Cassels--Tate
    pairing is alternating).  A #Sha_an that is not close to a perfect square
    is therefore evidence of a bug upstream, not of interesting mathematics.
  * #Sha_an must be close to a positive *integer*.  The distance to the nearest
    integer is a direct measure of the accumulated numerical error across the
    L-value, the period, and the regulator.

Both checks are applied here and reported, never silently rounded away.

See docs/GAP_MAP.md for what is actually proven; in particular the algebraic
rank used below is a *lower bound* obtained by point search unless it is
certified by the analytic rank being 0 or 1 (Gross--Zagier + Kolyvagin).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import mpmath

from bsdlab.curve import EllipticCurve

#: Largest deviation from an integer we accept before flagging #Sha_an as
#: numerically untrustworthy.  Deliberately loose: the point is to catch
#: factor-of-two and factor-of-p errors, not to certify precision.
SHA_INTEGRALITY_TOLERANCE = mpmath.mpf('1e-6')


@dataclass
class BSDData:
    """Every quantity in the strong BSD formula for one curve, plus provenance.

    A field set to ``None`` means *not computed*, and ``notes`` says why.  No
    field is ever filled with a plausible-looking placeholder.
    """

    label: Optional[str] = None
    ainvs: Optional[tuple] = None

    conductor: Optional[int] = None
    tamagawa: Optional[dict] = None          # {p: c_p}
    tamagawa_product: Optional[int] = None
    torsion_order: Optional[int] = None

    root_number: Optional[int] = None
    analytic_rank: Optional[int] = None
    leading_coefficient: Optional[mpmath.mpf] = None   # L^(r)(1)/r!

    real_period: Optional[mpmath.mpf] = None
    rank_lower_bound: Optional[int] = None
    generators: Optional[list] = None
    regulator: Optional[mpmath.mpf] = None

    sha_analytic: Optional[mpmath.mpf] = None
    sha_rounded: Optional[int] = None
    sha_integrality_error: Optional[mpmath.mpf] = None
    sha_is_square: Optional[bool] = None

    rank_certified: bool = False
    notes: list = field(default_factory=list)

    def note(self, msg: str) -> None:
        self.notes.append(msg)

    @property
    def consistent(self) -> bool:
        """True only if every automatic check passed.

        Requires: analytic rank equals the algebraic rank lower bound, #Sha_an
        is within tolerance of a positive integer, and that integer is a
        perfect square.
        """
        if None in (self.analytic_rank, self.rank_lower_bound, self.sha_rounded):
            return False
        return (
            self.analytic_rank == self.rank_lower_bound
            and self.sha_integrality_error is not None
            and self.sha_integrality_error < SHA_INTEGRALITY_TOLERANCE
            and self.sha_rounded >= 1
            and bool(self.sha_is_square)
        )


def _is_perfect_square(n: int) -> bool:
    if n < 0:
        return False
    r = mpmath.floor(mpmath.sqrt(n))
    for cand in (int(r) - 1, int(r), int(r) + 1):
        if cand >= 0 and cand * cand == n:
            return True
    return False


def sha_analytic(leading_coeff, torsion_order, real_period, regulator,
                 tamagawa_product):
    """Solve the strong BSD formula for #Sha.

        #Sha_an = L^(r)(1)/r! * |E_tors|^2 / (Omega * Reg * prod c_p)

    All arguments must be already-computed exact or high-precision values.
    Raises on a zero denominator rather than returning infinity.
    """
    denom = mpmath.mpf(real_period) * mpmath.mpf(regulator) * mpmath.mpf(tamagawa_product)
    if denom == 0:
        raise ArithmeticError(
            'Omega * Reg * prod(c_p) vanished; one of the three inputs is wrong. '
            f'Omega={real_period} Reg={regulator} prod_c={tamagawa_product}')
    return mpmath.mpf(leading_coeff) * (mpmath.mpf(torsion_order) ** 2) / denom


def analyse(E: EllipticCurve, label: Optional[str] = None, prec: int = 60,
            height_bound: int = 200) -> BSDData:
    """Run the whole pipeline on one curve and return every intermediate value.

    Never raises for a curve it cannot finish: the failure is recorded in
    ``notes`` and the corresponding fields stay ``None``.  This is so a batch
    run over thousands of curves reports its gaps instead of dying on one.
    """
    from bsdlab import reduction, lseries, periods, mordellweil

    E = E.minimal_model()
    d = BSDData(label=label, ainvs=E.ainvs)

    with mpmath.workdps(prec + 10):
        try:
            d.conductor = reduction.conductor(E)
            d.tamagawa = reduction.tamagawa_numbers(E)
            d.tamagawa_product = reduction.tamagawa_product(E)
        except Exception as exc:
            d.note(f'reduction failed: {exc!r}')
            return d

        try:
            d.torsion_order = E.torsion_order()
        except Exception as exc:
            d.note(f'torsion failed: {exc!r}')
            return d

        try:
            d.root_number = lseries.root_number(E, N=d.conductor)
        except Exception as exc:
            d.note(f'root number unavailable: {exc!r}')

        try:
            r, lead = lseries.leading_coefficient(E, N=d.conductor, prec=prec)
            d.analytic_rank, d.leading_coefficient = r, lead
        except Exception as exc:
            d.note(f'L-series failed: {exc!r}')
            return d

        try:
            d.real_period = periods.real_period(E, prec=prec)
        except Exception as exc:
            d.note(f'period failed: {exc!r}')
            return d

        try:
            rk, gens = mordellweil.rank_lower_bound(E, height_bound=height_bound,
                                                    prec=prec)
            d.rank_lower_bound, d.generators = rk, gens
            d.regulator = mordellweil.regulator(E, gens, prec=prec)
        except Exception as exc:
            d.note(f'Mordell-Weil failed: {exc!r}')
            return d

        # Gross-Zagier + Kolyvagin: analytic rank 0 or 1 certifies the algebraic
        # rank and the finiteness of Sha.  Above that, nothing is certified and
        # the rank below is only a lower bound from naive search.
        d.rank_certified = d.analytic_rank in (0, 1)
        if not d.rank_certified:
            d.note('rank not certified: analytic rank >= 2, Sha finiteness open')
        if d.rank_lower_bound != d.analytic_rank:
            d.note(f'RANK MISMATCH: analytic {d.analytic_rank} vs '
                   f'search lower bound {d.rank_lower_bound}')

        try:
            sha = sha_analytic(d.leading_coefficient, d.torsion_order,
                               d.real_period, d.regulator, d.tamagawa_product)
        except Exception as exc:
            d.note(f'Sha assembly failed: {exc!r}')
            return d

        d.sha_analytic = sha
        nearest = int(mpmath.nint(sha))
        d.sha_rounded = nearest
        d.sha_integrality_error = abs(sha - nearest)
        d.sha_is_square = _is_perfect_square(nearest)

        if d.sha_integrality_error >= SHA_INTEGRALITY_TOLERANCE:
            d.note(f'#Sha_an = {mpmath.nstr(sha, 15)} is not near an integer '
                   f'(error {mpmath.nstr(d.sha_integrality_error, 5)}); '
                   'a quantity upstream is wrong')
        elif not d.sha_is_square:
            d.note(f'#Sha_an rounds to {nearest}, which is not a perfect square; '
                   'Cassels pairing forbids this, so something upstream is wrong')

    return d
