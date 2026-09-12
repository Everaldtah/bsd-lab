# Findings

Results of running `python run.py` over every curve in `data/lmfdb_reference.json`.
All numbers below are produced by `bsdlab/`, which reads nothing from the
reference file; the reference is consulted only afterwards, to score the run.

Run configuration: `--prec 60 --height-bound 200`, 14 curves, conductors 11-5077,
analytic ranks 0 through 3. Raw output: [`data/results.json`](../data/results.json).

**14 / 14 curves fully consistent** — recomputed analytic rank matched the
LMFDB rank, the recomputed rank lower bound matched it too, and #Sha_an came out
integral and a perfect square in every case.

## Per-curve results

| label | N | w | rank | #Sha_an | L^(r)(1)/r! | Reg | Omega | prod c_p | tors |
|---|---|---|---|---|---|---|---|---|---|
| 11a1 | 11 | +1 | 0 | 1 | 0.2538418608 | 1 | 1.2692093042 | 5 | 5 |
| 11a2 | 11 | +1 | 0 | 1 | 0.2538418608 | 1 | 0.2538418608 | 1 | 1 |
| 11a3 | 11 | +1 | 0 | 1 | 0.2538418608 | 1 | 6.3460465213 | 1 | 5 |
| 14a1 | 14 | +1 | 0 | 1 | 0.3302236593 | 1 | 1.9813419560 | 6 | 6 |
| 15a1 | 15 | +1 | 0 | 1 | 0.3501507605 | 1 | 2.8012060846 | 8 | 8 |
| 32a3 | 32 | (*) | 0 | 1 | 0.6555143885 | 1 | 2.6220575542 | 1 | 2 |
| 37a1 | 37 | -1 | 1 | 1 | 0.3059997738 | 0.0511114082 | 5.9869172924 | 1 | 1 |
| 37b1 | 37 | +1 | 0 | 1 | 0.7256810619 | 1 | 2.1770431858 | 3 | 3 |
| 389a1 | 389 | +1 | 2 | 1 | 0.7593165002 | 0.1524601779 | 4.9804251217 | 1 | 1 |
| 571a1 | 571 | +1 | 0 | **4** | 1.7293650250 | 1 | 0.4323412562 | 1 | 1 |
| 571b1 | 571 | +1 | 2 | 1 | 0.9031605169 | 0.1772531402 | 5.0953146197 | 1 | 1 |
| 681a1 | 681 | -1 | 1 | 1 | 1.7148309964 | 0.2251689066 | 3.8078769893 | 2 | 1 |
| 681b1 | 681 | +1 | 0 | **9** | 1.8448152061 | 1 | 0.8199178693 | 4 | 4 |
| 5077a1 | 5077 | -1 | 3 | 1 | 1.7318499001 | 0.4171435587 | 4.1516879830 | 1 | 1 |

(*) 32a3 has additive reduction at 2 with v_2(N) = 5. The local root number there
is deliberately not implemented (see below), so the two-route cross-check is
unavailable for this curve; the Fricke numeric sign was used instead.

## What the integrality check actually shows

`#Sha_an` is a floating-point quotient. If any of the five inputs — the leading
L-coefficient, the real period, the regulator, the Tamagawa product, the torsion
order — were wrong by any amount, the quotient would land somewhere other than a
positive integer. Distance to the nearest integer:

| curve | |#Sha_an - nearest integer| |
|---|---|
| 11a2 | 0 (exact) |
| 11a1, 11a3, 14a1, 15a1, 32a3, 37a1, 37b1 | < 4e-61 |
| 389a1, 681a1 | < 8e-42 |
| 571a1, 571b1, 681b1 | < 9e-39 |
| 5077a1 | 2e-34 |

The tolerance is 1e-6 (`bsd.SHA_INTEGRALITY_TOLERANCE`). The worst curve clears
it by 28 orders of magnitude. Accuracy degrades with rank, as expected: the
regulator is the least precise input, and it enters only for rank >= 1.

## The two curves that matter most

571a1 (#Sha = 4) and 681b1 (#Sha = 9) are the load-bearing tests. On every other
curve in the set #Sha_an = 1, which means a bug that dropped or duplicated a
whole factor would still produce a perfect square and pass the structural checks
silently. On these two it would not: hitting exactly 4 and exactly 9 requires the
period, regulator, Tamagawa product and L-value to all be individually right.

## What is *not* established

- **The algebraic rank is a lower bound, not a rank.** `mordellweil.rank_lower_bound`
  is a naive point search over a height box. There is no descent, no saturation,
  no proof the generators found actually generate the Mordell-Weil group. That
  it agrees with the analytic rank on all 14 curves is evidence the search was
  adequate *here*, nothing more. This is the weakest link in the repo.
- **Rank >= 2 results are conditional.** Gross-Zagier and Kolyvagin certify rank
  and the finiteness of Sha only for analytic rank 0 or 1. For 389a1, 571b1 and
  5077a1, `rank_certified` is False and the reported #Sha is conditional on both
  strong BSD and on the algebraic rank equalling the analytic rank.
- **#Sha_an is a prediction, not a computation of Sha.** Nothing here computes
  the Tate-Shafarevich group. The integrality and square checks detect internal
  inconsistency; they do not verify BSD.
- **Local root numbers at additive primes are not implemented.** The Halberstadt
  and Kraus tables were not derived rigorously here, so `lseries` raises
  `NotImplementedError` rather than guessing. Affected curves fall back to the
  Fricke numeric sign, which is correct but removes the independent cross-check.
- **14 curves is a small sample.** Conductor 5077 is not large.

## Independence claim

No Sage, PARI, cypari or eclib is used anywhere, on any path. The only
dependencies are sympy, mpmath and numpy. The period lattice is computed by a
hand-rolled AGM (`periods._agm`), not `mpmath.ellipk`. The point of agreeing
with LMFDB is that the agreement is between two genuinely separate
implementations; if this repo shared PARI's lineage, agreement would be
tautological rather than evidence.

See [`METHODOLOGY.md`](METHODOLOGY.md) for what each invariant is trusted to and
where each computation fails.
