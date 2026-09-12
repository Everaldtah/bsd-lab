# A plan of attack for the general case

This document turns the boundary map of [`BARRIER.md`](BARRIER.md) into a
program: ranked routes, concrete milestones, and the experiments this lab can
run that no one else can. It is grounded in the data collected so far — the
14-curve full-formula verification ([`FINDINGS.md`](FINDINGS.md)), the
39,790-curve enumeration, and the 559-curve local-data verification against
the LMFDB oracle.

It is not a proof sketch. No accepted program for general BSD exists. What
follows is the strongest attack we can specify today, with the lab's
computational leverage attached to the points where it can falsify or
corroborate candidate constructions.

## 1. The target, stated precisely

`BSD(E/Q)` for every elliptic curve over Q, unconditionally:

1. **Rank equality**: `rank E(Q) = ord_{s=1} L(E,s)`.
2. **The exact formula** at `r = rank`:

   ```
   L^(r)(E,1)/r!  =  (prod_p c_p . Omega . Reg(E) . #Sha) / |E_tors|^2
   ```

"General case" means: all analytic ranks `r >= 2` (already theorems for
`r <= 1`, by D4 below), and all primes `p` in the formula at once, including
`p = 2, 3`, supersingular and additive reduction.

## 2. What our data establishes, and why it matters strategically

Each fact below was recomputed from scratch and checked against an
independent oracle; each has a direct consequence for the attack.

- **The formula is rigid.** 14/14 curves consistent; `#Sha_an` lands on a
  positive integer with 28 orders of magnitude to spare over the `1e-6`
  tolerance. Every one of the five inputs is individually pinned — a wrong
  period, regulator, Tamagawa product or L-value anywhere destroys
  integrality. *Consequence:* any proposed identity in Route B (a
  "Gross–Zagier for the second derivative") can be tested numerically to
  several digits in minutes; false constructions die instantly. This is the
  lab's most important property as an instrument.

- **The hard analytic primes are the active ones.** The only curves in our
  set with nontrivial Sha are `571a1` (`#Sha = 4 = 2^2`) and `681b1`
  (`#Sha = 9 = 3^2`) — both rank 0, both *certified* as to rank and
  finiteness, and both failing at exactly the primes where the known p-part
  theorems carry hypotheses (`p = 2`, `p = 3`). The general case's weakest
  flank is demonstrably live in the data, not a technicality.

- **The rank >= 2 frontier is calibrated.** `389a1` (rank 2, the smallest
  conductor at rank 2), `571b1` (rank 2) and `5077a1` (rank 3) have
  regulators computed from independently found points agreeing with
  `L^(r)(1)/r!` to `1e-9`. *Consequence:* these three curves are the
  reference targets for any construction claiming to produce rank-2 points:
  success means reproducing `Reg = det <P_i, P_j>` numerically, and we can
  score it.

- **Parity holds everywhere we can check it, and never suffices.** Root
  numbers came out consistent with rank parity on all 14 curves (and the
  whole enumerated set where implemented). No case exists in our data where
  parity alone pins the rank — as expected: it is a congruence mod 2
  (BARRIER.md §3.3). *Consequence:* any rank-equality proof must pass
  through a mechanism that upgrades the congruence to an equality; there is
  no empirical shortcut around this.

- **The local side is conquered.** 559 curves verified against the oracle
  with zero mismatches, Ogg's formula asserted (never violated) at every
  prime of every curve, 39,790 curves enumerated. Conductor, Tamagawa
  numbers, Kodaira types: solved, in an independent implementation.
  *Consequence:* all residual difficulty in the general case is global —
  the L-function side and the Mordell–Weil side. Nothing about the local
  factors blocks any route below.

## 3. Architecture of a proof: the load-bearing decomposition

Full BSD for all `E/Q` is equivalent to the conjunction of five claims:

- **D1 Modularity**: every `E/Q` is attached to a weight-2 newform.
  **DONE** (Wiles et al.; Breuil–Conrad–Diamond–Taylor).
- **D2 Parity**: `corank_p Sel_{p^inf}(E/Q) = ord L(E,1) (mod 2)` for every
  p. **DONE** (Nekovář; Dokchitser–Dokchitser 2010).
- **D3 p-part of the formula**: for each prime p, the p-part of the exact
  formula, given rank equality. **PARTIAL**: analytic rank 0 and 1, `p >= 5`
  good ordinary, residual irreducibility hypotheses (Kato 2004; Skinner–Urban
  2014; Jetchev–Skinner–Wan 2017).
- **D4 Rank equality**: **DONE** for `ord L <= 1` (Gross–Zagier 1986 +
  Kolyvagin 1990). **OPEN** for `ord L >= 2`.
- **D5 Sha finiteness**: same boundary as D4 — Kolyvagin gives it exactly
  where his Euler system runs. **OPEN** at rank >= 2.

The closure is: D3 for all p simultaneously + D4 + D5 imply the full formula
(the p-parts control each prime of `#Sha` in turn; D5 makes the product
finite; D4 makes `r` meaningful on both sides).

So the general case needs exactly three advances:

- **A.** D3 without hypotheses, at every p (including 2, 3, supersingular,
  additive) — *Route A*.
- **B.** D4 at rank 2 — *Route B*, the genuine barrier.
- **C.** D5 at rank 2 — expected to fall out of B by the same mechanism by
  which Kolyvagin's Euler system delivers both rank and finiteness at rank 1.

## 4. Route A — the p-part program (the tractable flank)

**Goal.** For `E/Q` of analytic rank `<= 1`: prove the exact formula one
prime at a time, for *every* prime, with no hypotheses. This route has no
conceptual wall — every obstruction is a technical one on a known machine —
which is exactly why it is where publishable progress is realistic.

Milestones, in dependency order:

- **A1. Supersingular primes `p >= 5`.** The Skinner–Urban argument needs
  p-adic families with good ordinary vanishing; at supersingular primes the
  correct objects are Pollack's `±` p-adic L-functions and Kobayashi's
  `±` Selmer groups. Target: the main conjecture for each sign branch,
  hence the p-part on both `±` sides, glued. The lab's contribution: our
  census marks how common supersingular `p = 2, 3` is in the enumerated
  range, i.e. how far A1 alone goes before A2 becomes mandatory.
- **A2. The primes `p = 2` and `p = 3`, ordinary.** The hardest local
  analytic primes of the Iwasawa machine. Note the empirical anchor from
  §2: `571a1` and `681b1` say that when Sha is nontrivial, it is *at these
  primes* (2 and 3) that it is. Route A without A2 would certify almost
  every curve in our set and miss precisely the interesting ones.
- **A3. Additive reduction at p / p dividing the conductor.** Requires the
  ramified-level Iwasawa theory. The same local difficulty is visible on
  the computational side: at `32a3` (additive at 2, `v_2(N) = 5`) even our
  two-route root-number cross-check breaks and falls back to the numeric
  sign (FINDINGS.md) — the local p=2 analytic invariants are the shared
  bottleneck of both the theoretical and the computational attack.
- **A4. Reducible residual representations.** Drop irreducibility: the
  reducible cases degenerate toward CM (handled classically by Rubin's
  Euler system) and Eisenstein level-lowering cases. Combine.
- **A5. Assembly.** Conjunction over all p of A1–A4 plus D4-at-rank-<=1
  gives: *BSD holds unconditionally for every curve of analytic rank <= 1.*
  That is the headline deliverable of this route — a theorem that would be
  recognized as the completion of the rank 0/1 chapter of the conjecture.

**Risk profile.** High effort, zero conceptual risk. Each of A1–A4 is an
active subfield; none requires an object that does not currently exist.
This is also the route where our lab contributes only marginally (data, not
proofs) — the leverage is all in Route B.

## 5. Route B — the rank-2 object (the real problem)

### 5.1 The specification

What must exist, for the barrier to fall: a construction attaching to `E/Q`
with `ord L(E,1) = 2` two cohomology classes `P_1, P_2` together with

```
L''(E,1)/2!  =  (explicit local correction factors) . det <P_i, P_j>_hat
```

(where `<,>_hat` is the Néron–Tate height pairing), *and* an Euler-system
argument converting `det != 0` into: `P_1, P_2` span a rank-2 subgroup,
`#Sha[E]` is bounded by the index, Sha is finite. This is the rank-2
instance of the machine whose rank-1 instance is Gross–Zagier + Kolyvagin
and whose rank-0 instance is Kato. Its absence is BARRIER.md §3.1–3.2
restated positively.

### 5.2 Candidate constructions, ranked

**(a) Quadratic twists of Heegner points — available, provably insufficient,
run first anyway.** For rank-2 `E`, twists `E ⊗ χ` with flipped sign have
odd analytic rank; Gross–Zagier points exist on them, but their traces to
`E(Q)` are torsion (the rank-2 information lives in `E` over the twist's
quadratic field, canceling in cohomology). *Why bother:* the obstruction is
precise, and it defines the shape of the replacement — points over **real**
quadratic fields, where no cancellation occurs. *Lab experiment L6(i):* the
twist census on the 39,790-curve set makes the obstruction quantitative
(rank distribution of `E ⊗ χ_d` for the rank-2 reference curves across all
fundamental `d` with `|d| <= D`).

**(b) Stark–Heegner / Darmon points over real quadratic fields.** The
conjectural construction (Darmon, 2000s): rigid-analytic integration on
`X_0(N)` produces points on `E` over ring fields of real quadratic `K`,
conjecturally satisfying a GZ-type formula for `L'(E/K, 1)` — and for
suitable real `K`, `L(E/K, s) = L(E,s) L(E ⊗ χ_K, s)` picks up rank-2
behavior of `E`. Missing: the conjecture itself (global compatibility of
the p-adic integration), norm relations, any Euler system. *First step is
numerical, and it is decisive:* compute candidate points for `389a1`, run
the height pairing, compare `det <P_i,P_j>` against
`L''(389a1, 1)/2! . |tors|^2 / (prod c_p . Omega)`. Our `Reg(389a1) =
0.1524601779...` to 10+ digits is exactly the oracle this test needs. If
the determinant misses, the construction as conjectured is wrong — the
cheapest possible kill.

**(c) Beilinson–Flach Euler systems (`E ⊗ E'`).** The one *existing*
Euler system naturally attached to a Rankin–Selberg product
(Lei–Loeffler–Zerbes): proven classes, proven norm relations. The bottom
class is pairing-shaped rather than point-shaped — the only machine in the
literature whose bottom object is 2-dimensional in the relevant sense.
Missing: any argument converting nonvanishing of `L''(E,1)` into
nonvanishing of a Beilinson–Flach element for a well-chosen partner form
`E'`. *Step:* the lab can search, among forms of the right levels in our
database, for partner pairs where the Rankin–Selberg L-value degenerates to
`L''(E,1)` — the arithmetic side would then be a height pairing of two BF
elements, testable by the same regulator oracle.

**(d) Diagonal cycles / triple products (Darmon–Rotger).** Proven formula:
generalized Heegner cycles on triple products have p-adic heights governed
by triple-product central L-derivatives. If `E` of analytic rank 2 appears
as one factor of a triple `(g, h, E)` whose triple-product derivative
encodes `L''(E,1)`, this is a *theorem-shaped* source of the needed
2-dimensional class (the cycle's projections to `E`). Missing: the partner
selection, and everything global. *Step:* L6(ii), a partner hunt over the
enumerated database.

**(e) Two-variable p-adic L-functions with a second-derivative regulator.**
Bertolini–Darmon-style Hida-family anticyclotomic machinery: a 2-variable
p-adic L-function whose mixed second derivative at the trivial character
should equal a p-adic height pairing of *two* classes. The p-adic world has
historically been more forgiving than the archimedean one. Missing: the two
classes — i.e. (b) or (d) supplies them and (e) supplies the formula; they
pair naturally.

### 5.3 What rank >= 3 needs

Not assumed to follow. Each rung so far (rank 0: Kato; rank 1:
Gross–Zagier) required a new object; `5077a1` (rank 3, in our verified set)
is the reference target one rung beyond the barrier. The honest plan: once
the rank-2 mechanism of 5.2 is visible, its generalization — r-dimensional
wedge of classes against `L^(r)/r!` — is dictated rather than guessed. The
lab's `5077a1` row (`Reg = 0.4171435587...`) is the calibration target for
that rung.

## 6. Route C — corank-to-rank at 2 (closing the loop)

A p-converse theorem at rank 1 exists (Skinner): Selmer corank 1 plus
`L' != 0` forces rank 1. The rank-2 analogue would read: corank 2 plus
`L'' != 0` forces rank 2. Chasing the proof shape, the missing step is
exactly the vanishing of `Sha[2]`'s corank contribution — D5 — which is
Route B's second deliverable (§5.1). Route C is thus Route B restated as a
clean implication, listed for completeness and to close off the illusion of
an independent shortcut.

## 7. The lab's program (what we build, in order)

The lab cannot prove any of A/B, but it can be the sharpest available
*instrument* pointed at them, while producing unconditional deliverables of
its own. Ordered by dependency, not by glamour:

- **L1. Rigorous error bounds (interval certificates).** Replace the `1e-6`
  integrality heuristic with: explicit truncation bounds for the
  L-series at `s = 1 + delta`-style evaluation points, Archimedean error
  control on the AGM period, and certified error on the regulator heights.
  Deliverable: the 14-curve table becomes a list of machine-checkable
  certificates (`#Sha_an ∈ Z` as an interval statement, not a heuristic).
  Everything below depends on trusting numbers farther than we can
  currently *prove* they are right.
- **L2. Descent.** Implement 2-descent over Q (2-coverings à la Birch–
  Swinnerton-Dyer's original descent, exact, no oracles) plus the
  isogeny `φ`-descent for curves with rational torsion. This converts
  `mordellweil.rank_lower_bound` — the documented weakest link — into a
  two-sided bound `rank <= dim Sel^2 − dim E(Q)[2]`. Deliverable: on the
  14 reference curves the algebraic rank becomes known independently of
  L; for `389a1`, `571b1`, `5077a1` the remaining conditionality shrinks
  to Sha + the exact formula only.
- **L3. The regulator oracle for candidate points.** A one-call harness:
  `test_generators(E, [P_1, P_2])` → the height-pairing determinant vs
  `L^(r)(1)/r! · |tors|² / (prod c_p · Omega)`, with L1-grade error bars.
  This is the lab's unique leverage on the actual Millennium barrier:
  every Route-B candidate (b), (c), (d), (e) dies or lives here, at a cost
  of minutes, before anyone spends proof-years on it.
- **L4. Scale the certified set.** 14 → ~500 curves, full formula,
  prioritizing curves with `#Sha_an > 1` (each pins another prime's
  p-part; `571a1` and `681b1` showed these are where formula errors
  surface). Depends on L1 (else we are producing heuristics at scale).
- **L5. Local root numbers at additive primes.** Derive the Halberstadt
  and Kraus tables in-house (from the same Tate's-algorithm data we now
  compute exactly), restoring the two-route sign cross-check for the
  `32a3`-type additive cases. Closes the last gap between what we compute
  and what the local theory knows.
- **L6. Experiments in service of Route B.**
  - *L6(i) Twist census*: rank distribution of `E ⊗ χ_d` for the three
    frontier curves across fundamental discriminants, quantifying the §5.2(a)
    obstruction (predicts: never rank 2 in the twist family with sign +1
    trace-recovering two independent rational points).
  - *L6(ii) Triple-product partner hunt*: over the enumerated database,
    search partner pairs `(g, h)` for the frontier curves whose
    triple-product central derivative degenerates to `L''(E,1)` (§5.2(d)).
  - *L6(iii) Rigidity map*: how far the integrality clearance survives
    precision reduction per curve/rank — tells theorists at which conductor
    and rank the formula's numerical rigidity is actually testable.

## 8. Sequencing and milestones

- **Phase 1 (lab, deterministic):** L1 → L2 → L5, then L4. All
  engineering; all unconditional deliverables; nothing speculative.
- **Phase 2 (lab, parallel):** L3 as soon as L1 lands; L6 experiments on
  the 39,790-curve base. Publishable on their own (the certified table and
  the twist census are contributions regardless of Route B's fate).
- **Phase 3 (theory, long-horizon):** Route A is the profession's
  engineering flank — tractable, publishable at every milestone.
  Route B is the genuine attempt; its next step after Phase 2 is *chosen by
  what L3/L6 falsify*: candidates that fail numerically are dead cheaply;
  survivors earn theory investment.
- **Rung definitions of success:**
  - R1 = BSD unconditional for all analytic rank ≤ 1 curves (Route A,
    profession-scale).
  - R2 = Sha finite for one rank-2 curve (would be the first in history;
    Route B byproduct).
  - R3 = the general case.

## 9. Honest risk statement

The probability that this plan closes R3 is effectively zero on any
planned timescale — it is a Clay Millennium Problem and no accepted program
exists; §5.2 is a ranked list of leads, not a ladder. The plan's value is
in what it secures unconditionally regardless:

1. L1–L5 are deterministic deliverables of a fully independent BSD
   recomputation, stronger than anything currently in the repo.
2. L3/L6 make this lab the cheapest decisive falsifier for candidate
   rank-2 constructions — a role with real expected value, since most
   conjectural constructions in §5.2 are wrong and nobody else can kill
   them in minutes.
3. The dataset itself (39,790 curves, 559 verified, Sha distribution)
   is empirical evidence about the shape of the conjecture — the same role
   the original Birch–Swinnerton-Dyer computations played in *forming* it.

We are not close to the general case. We are exactly as close as the
boundary allows, the boundary is documented to the statement level
(BARRIER.md §5), and every step this plan takes is either unconditional
progress or a sharper instrument at the wall.

