# BSD Lab — an independent computational attack on Birch and Swinnerton-Dyer

> **Honest framing, up front.** The Birch and Swinnerton-Dyer conjecture is a Clay
> Millennium Prize problem. It is open. This repository does **not** claim a proof,
> and any file here that appears to claim one is a bug — please open an issue.
>
> What this repository *is*: a from-scratch, dependency-light, fully auditable
> computational laboratory for BSD, plus the datasets it produces and a precise map
> of what is actually proven versus what is open.

## Why build it from scratch

Essentially all published numerical evidence for BSD flows through one software
lineage — Cremona's `mwrank`/`eclib` and PARI/GP, surfaced via Sage and the LMFDB.
That software is excellent. It is also a **single point of correlated failure**: if
a subtle bug existed in, say, the period lattice integration or the Tate's-algorithm
branch for additive reduction at 2, every "verification of BSD" downstream would
inherit it and no amount of extra curves would reveal it.

This project implements the whole pipeline independently, in plain Python, from the
definitions:

```
E/Q  ──►  Tate's algorithm  ──►  conductor N, Tamagawa numbers c_p
     ──►  point counting     ──►  a_p, and the L-series coefficients a_n
     ──►  approximate functional equation ──►  L(E,s) near s = 1
     ──►  analytic rank r_an and leading coefficient L^(r)(1)/r!
     ──►  AGM  ──►  real period Ω
     ──►  descent / search  ──►  Mordell-Weil rank, generators, regulator Reg
     ──►  strong BSD  ──►  #Sha_an  =  L^(r)(1)/r! · |E_tors|² / (Ω · Reg · ∏ c_p)
```

The payoff is a genuinely **independent check**. Where we agree with the LMFDB, that
agreement is evidence of a kind the field currently does not have. Where we disagree,
one of the two is wrong, and that is worth knowing either way.

## The three things this repo produces

**1. `bsdlab/` — the library.** Pure Python (`sympy`, `mpmath`, `numpy` only). Every
function states its mathematical definition and its failure modes in the docstring.
No silent fallbacks: a quantity that cannot be computed rigorously is returned as
`None` with a reason, never as a plausible-looking float.

**2. `data/` — the datasets.** Machine-readable results, each row carrying the
provenance and the numerical certainty of every quantity in it. Headline runs:

- *Cross-validation*: our computed `(N, c_p, Ω, r_an, #Sha_an)` against LMFDB values,
  curve by curve, with every discrepancy itemised rather than averaged away.
- *Quadratic twist families*: rank and `#Sha_an` distributions over twists of a fixed
  curve, testing Goldfeld's conjecture (density ½ rank 0, ½ rank 1) and the
  observed-but-unproven growth rate of Sha.
- *Sha spectrum*: which integers actually occur as `#Sha`, and at what frequency.
  Cassels' pairing forces a perfect square; deviation from a square is a hard error
  signal that something upstream miscomputed.

**3. `docs/` — the gap map.** What the Gross–Zagier / Kolyvagin machinery actually
buys (rank ≤ 1, and only then), what Skinner–Urban and the p-converse theorems
added, and precisely where rank ≥ 2 breaks every known method. Written so the open
boundary is a list of concrete statements, not a vibe.

## Status

Working pipeline. The 12-curve self-test passes (100 checks, 0 failures) and the
batch run reproduces rank and #Sha for all 14 reference curves, conductors
11-5077, including both nontrivial-Sha cases (571a1 #Sha = 4, 681b1 #Sha = 9).
See `docs/FINDINGS.md` for the results and `docs/METHODOLOGY.md` for how each
number is computed and how far it is trusted.

## Reproducing

```bash
pip install sympy mpmath numpy
python -m bsdlab.selftest        # verifies the library against known curves
python run.py --help             # batch driver -> data/results.json
python run.py --prec 60          # reproduces docs/FINDINGS.md
```

Every dataset in `data/` is regenerable from a single command recorded in its header.

## Method note

The mathematics is orchestrated by Claude Opus and executed by a fleet of
Qwen3-27B agents; every mathematical claim and every line of library code is
verified against known values before it is committed. The self-test suite in
`bsdlab/selftest.py` is the contract: it pins our pipeline to curves whose
invariants are independently documented, and nothing ships while it fails.

## Licence

MIT.
