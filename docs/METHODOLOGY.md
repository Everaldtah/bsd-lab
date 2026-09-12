# Methodology: how each BSD invariant is computed

This document describes, invariant by invariant, what the code in `bsdlab/`
computes, which function computes it, how far the result can be trusted, and
when it gives up. It complements `docs/GAP_MAP.md` (what is proven) and the
docstrings (line-level detail). No Sage, PARI/GP, cypari or eclib is used
anywhere in the pipeline; the runtime dependencies are `sympy`, `mpmath` and
`numpy` only. Independence from the Cremona/PARI software lineage is the
entire point of the project: agreement with LMFDB values is treated as
independent evidence, disagreement as a bug hunt.

## What is computed

The pipeline assembles every quantity in the strong BSD formula

    #Sha = ( L^(r)(1)/r! * |E_tors|^2 ) / ( Omega * Reg * prod_p c_p )

for a curve E/Q given by integral a-invariants:

- minimal model, discriminant, j-invariant (`curve.EllipticCurve.minimal_model`);
- conductor N and Tamagawa numbers c_p via Tate's algorithm (`reduction`);
- torsion subgroup order via Lutz-Nagell + Mazur (`curve.EllipticCurve.torsion_order`);
- a_n coefficients, root number, analytic rank, L^(r)(1)/r! (`lseries`);
- real period Omega via a hand-rolled AGM (`periods.real_period`);
- Mordell-Weil rank lower bound, generators, regulator (`mordellweil`);
- the analytic order of Sha and the consistency checks (`bsd.analyse`).

Everything exact is computed in exact integer/Fraction arithmetic; mpmath
enters only where logarithms, exponentials or quadrature are unavoidable.

## Conductor and Tamagawa numbers

**Definition.** The conductor N = prod_p p^f_p, where f_p is the conductor
exponent at each prime of bad reduction, determined by the Kodaira type of
the special fibre of the Neron model. The Tamagawa number c_p is the number
of F_p-rational components of the special fibre. Ogg's formula
f_p = v_p(Delta_min) - m_p + 1 (with m_p the number of geometric components)
links the three and is a theorem.

**Algorithm.** `reduction.local_data` runs Tate's algorithm (Silverman,
Advanced Topics, IV.9; Cremona, section 3.2) on the global minimal model from
`curve.EllipticCurve.minimal_model`, handling good, multiplicative (I_n),
additive (II, III, IV, I_0*, I_n*, IV*, III*, II*) cases, with the p = 2 and
p = 3 branches written out separately where the generic closed forms
degenerate. Non-minimal models trigger the divide-a_i-by-p^i retry loop.
`reduction.conductor` and `reduction.tamagawa_numbers` / 
`reduction.tamagawa_product` aggregate the local data over the bad primes
(listed by `curve.EllipticCurve.bad_primes`). Ogg's formula is asserted for
every prime inside `reduction.local_data`: a violation raises
`ArithmeticError`, since it can only mean an implementation bug.

**Trusted to.** Exact integers, unconditionally: the entire computation is
integer arithmetic. Internal consistency is guarded by the Ogg assertion.

**Fails when.** `reduction.local_data` raises `ValueError` if p is not prime,
`ArithmeticError` if Tate's algorithm fails to terminate within its iteration
bound, if an intermediate divisibility check fails (each such raise indicates
a bug, not a property of the curve), or if Ogg's formula is violated. If
`curve.EllipticCurve.minimal_model` itself fails, the code falls back to the
input model and relies on the p-minimality supplied by the retry loop.

## L-function, root number, analytic rank

**Definition.** L(E,s) = sum a_n n^-s with a_n the Fourier coefficients of
the newform attached to E; a_p = p + 1 - #E(F_p) at good primes, +1/-1 at
split/non-split multiplicative primes, 0 at additive primes, extended
multiplicatively. The completed function Lambda(s) = N^(s/2)(2pi)^-s
Gamma(s) L(E,s) satisfies Lambda(s) = w Lambda(2-s) with w = +-1 the root
number. The analytic rank is the least r with L^(r)(1) != 0.

**Algorithm.** `lseries.an_coefficients` fills a_1..a_nmax by a
smallest-prime-factor sieve, prime powers by the standard recurrences, using
`curve.EllipticCurve.ap` (brute-force point count mod p) at good primes.
Central values and derivatives come from the exact rapidly-convergent
identity in the module docstring: `lseries.l_value`, 
`lseries.l_derivative` (via `lseries._completed_moments`, built on
incomplete-gamma kernels G_r), and `lseries.analytic_rank` / 
`lseries.leading_coefficient`, which scan r of the parity forced by w and
call a value nonzero above the threshold 10^-(prec/2). The root number is
`lseries.root_number`: w = -prod_{p|N} w_p (leading -1 = archimedean place;
split multiplicative w_p = -1, non-split w_p = +1), cross-checked against
`lseries.root_number_numeric`, which recovers w from the Fricke involution
as w = f(i/(2 sqrt(N))) / (4 f(2i/sqrt(N))) and insists the ratio be within
0.01 of +-1. The Halberstadt/Kraus additive local tables are deliberately
not implemented rather than guessed: `_local_root_number` raises
`NotImplementedError` for additive reduction at p = 2, p = 3 and p >= 5.
`lseries.sign` tries the local product and falls back to the numeric sign.

**Trusted to.** a_n and w are exact (integer arithmetic / exact q-series
identity). L^(r)(1)/r! carries roughly `prec` significant digits (default
60) after truncation control; analytic rank is correct to working precision
only — a true leading term below 10^-(prec/2) would be misread as zero, and
raising `prec` is the remedy. When both root number routes are available
they must agree or `lseries.root_number` raises.

**Fails when.** `lseries.root_number` raises `NotImplementedError` at
additive primes and `ArithmeticError` on a local/numeric sign mismatch.
`lseries.root_number_numeric` raises `ArithmeticError` if the Fricke ratio
is not within 0.01 of +-1. `lseries.analytic_rank` raises `ArithmeticError`
if no derivative up to `max_rank` (default 4) is nonzero, or if the first
nonzero derivative has the wrong parity for w. `_resolve_N` raises
`ImportError` if no conductor is supplied and `bsdlab.reduction` is missing.

## Real period

**Definition.** Omega is the BSD/LMFDB real period: the least positive real
period of the invariant differential dx/(2y + a1 x + a3), multiplied by the
number of real components — Omega = w1 when the discriminant is negative
(one component), Omega = 2*w1 when positive (two).

**Algorithm.** `periods.real_period` calls `periods.period_lattice`, which
computes the half-periods of the completed model
Y^2 = 4x^3 + b2 x^2 + 2 b4 x + b6 as w_re = pi / (2 AGM(a,b)) and
w_im = pi i / (2 AGM(a,c)) with a, b, c square roots of root differences
(from `periods._cubic_roots`, mpmath `polyroots`). The AGM itself,
`periods._agm`, is hand-rolled: a plain iterate (a,b) -> ((a+b)/2,
sqrt(ab)) in mpmath arithmetic to a tolerance of 10^-(dps-5). mpmath's
`ellipk` is NOT on the primary path — no CAS elliptic machinery is used;
only `polyroots`, `quad` and elementary mpmath functions appear in the
module. `periods.number_of_real_components` supplies the factor 1 or 2.
`periods.elliptic_logarithm` (Cremona 3.7 numerical inversion) also lives
here, though the height code does not need it.

**Trusted to.** Approximately `prec` digits (default 60; internally
computed at prec + 10 and rounded). The Omega = w1 vs 2*w1 convention is
asserted in the test suite, not merely documented, because a factor of 2
here propagates straight into #Sha.

**Fails when.** `periods._cubic_roots` raises `ValueError` on a singular
curve (repeated root). The complex AGM only converges when the callers'
sign choices keep Re(a * conj(b)) > 0; the branch selection in
`periods._half_periods` handles both root orderings, and a wrong ordering
would show up as a wrong period in the tests rather than an exception.

## Mordell-Weil rank and regulator

**Definition.** The algebraic rank is the rank of E(Q); the regulator Reg is
the determinant of the height pairing matrix
<P_i,P_j> = (hhat(P_i + P_j) - hhat(P_i) - hhat(P_j))/2 on a basis of the
free part, with hhat the Neron-Tate canonical height, normalised so that
hhat(nP) = n^2 hhat(P).

**Algorithm.** `mordellweil.rank_lower_bound` searches for rational points
with `mordellweil.search_points` (every x = a/b^2 in a box, y decided by an
exact integer square test on the discriminant N/b^6), extracts a maximal
independent subset with `mordellweil.independent_points` (height-pairing
Gram determinant pivoting at `INDEPENDENCE_TOLERANCE` = 1e-25 relative), and
returns (rank found, generators). `mordellweil.regulator` computes the Gram
determinant of the generators' height pairings (1 for the empty list).
`mordellweil.canonical_height` is the sum of Silverman (1988) local heights:
the archimedean place by Tate's rapidly convergent theta-series
(`_archimedean_local_height`, ~0.51 * prec terms, translation trick for
|x| < 1/2), the finite places by the valuation formulas on the minimal
model (`_non_archimedean_local_height`). Torsion is decided by the group
law (order <= 16 by Mazur), never by a numerical threshold.

**Trusted to.** The rank returned by `mordellweil.rank_lower_bound` is a
LOWER BOUND from a naive point search, not a proven rank and not a descent:
curves whose generators lie outside the searched box have strictly larger
true rank. This is the weakest link in the pipeline. Heights and regulator
are good to ~`prec` digits; independence of the returned generators is a
numerical judgement at 1e-25, safe for searched points of height O(1).

**Fails when.** Nothing here returns None; the failure mode of the lower
bound is silence (too few points found), which is why `bsd.analyse`
cross-checks it against the analytic rank and records a RANK MISMATCH
note. `mordellweil._to_minimal` raises `ArithmeticError` if a point fails
to map onto the minimal model (a bug signal). Very large height bounds
make `search_points` slow, not wrong.

## Assembling #Sha

**Definition.** The analytic order of Sha is what the strong BSD formula
predicts: #Sha_an = ( L^(r)(1)/r! * |E_tors|^2 ) / ( Omega * Reg *
prod_p c_p ). It is a prediction conditional on strong BSD and on the five
inputs, not a computation of the Tate-Shafarevich group.

**Algorithm.** `bsd.sha_analytic` evaluates the quotient directly from the
five already-computed inputs. `bsd.analyse` runs the whole pipeline for one
curve into a `bsd.BSDData` record: conductor and Tamagawa product
(`reduction`), torsion order (`curve`), root number and analytic rank with
leading coefficient (`lseries`), real period (`periods`), rank lower bound
and regulator (`mordellweil`), then #Sha_an. Rank certification:
Gross-Zagier plus Kolyvagin certify rank and Sha finiteness only for
analytic rank 0 or 1, so `rank_certified = (analytic_rank in (0, 1))`.

**Trusted to.** #Sha_an is a floating-point number whose distance to the
nearest integer measures accumulated error across the L-value, the period
and the regulator; `sha_rounded` is that nearest integer, trustworthy only
when `sha_integrality_error` < 1e-6 (`SHA_INTEGRALITY_TOLERANCE`). Two
structural checks guard it: #Sha must be a positive integer, and by
Cassels' alternating pairing on Sha it must be a perfect square — a
non-square is treated as a hard error signal that something upstream is
wrong, recorded in `notes` (never silently rounded away).

**Fails when.** `bsd.sha_analytic` raises `ArithmeticError` if
Omega * Reg * prod c_p vanishes. `bsd.analyse` itself never raises for a
curve it cannot finish: failures are recorded in `BSDData.notes` and the
corresponding fields stay None. `BSDData.consistent` is False unless the
analytic rank equals the rank lower bound, #Sha_an is integral within
tolerance, and the rounded value is a perfect square >= 1.

## Precision and error control

- Default working precision is 60 significant digits (`lseries.DEFAULT_PREC`,
  the `prec = 60` defaults); internals compute at prec + 10 and round.
- L-series truncation: the smoothed sums carry exp(-2 pi n / sqrt(N)), so
  `n_terms = 10 sqrt(N) + 100` puts the first omitted term near exp(-20 pi)
  ~ 4e-28, below working precision for prec <= 50 or so; raise n_terms for
  higher precision.
- Analytic-rank zero detection: |L^(r)(1)/r!| <= 10^-(prec/2) counts as
  zero. This is a working-precision criterion, not a proof of vanishing.
- Period and regulator accuracy is validated end-to-end by the #Sha
  integrality check: an error anywhere in the chain surfaces as a
  non-integral (or non-square) #Sha_an.
- Structural cross-checks that raise rather than guess: Ogg's formula at
  every bad prime; the local-vs-Fricke root number agreement; the parity of
  the first nonzero L-derivative against w; exact square tests in the point
  search.

## Known limitations

- The algebraic rank is only ever a lower bound from naive search; there is
  no descent, no saturation at primes, no proof that the found generators
  generate. `mordellweil.rank_lower_bound` is the weakest link in the repo.
- Gross-Zagier + Kolyvagin certify rank and Sha finiteness only for analytic
  rank 0 or 1. For analytic rank >= 2 nothing is proven: the reported #Sha
  is CONDITIONAL on the algebraic rank equalling the analytic rank and on
  strong BSD itself, and `rank_certified` is False there.
- The analytic rank can be misread if the true leading term lies below the
  10^-(prec/2) threshold or the rank exceeds `max_rank` (default 4).
- Local root numbers at additive primes are not implemented
  (NotImplementedError at p = 2, 3 and p >= 5); the Fricke numeric sign
  covers them, at the cost of losing the independent two-route cross-check
  for such curves.
- Point counting `curve.EllipticCurve.count_points_mod_p` is O(p) brute
  force; large-conductor curves are slow (correct, but slow).
- #Sha_an is a prediction, not a computation of Sha; the square/integrality
  checks detect internal inconsistency, they do not verify BSD.
