# Route B / W3b — Stark–Heegner (Darmon) points on 389a1 over real quadratic K

Status: SPEC, not implemented. Targets PLAN.md §5.2(b). Oracle: Reg(389a1) =
0.1524601779... (FINDINGS.md, `bsdlab.oracle.test_generators`).

Contents:
1. The choice of K — sign analysis, genus decomposition, primary/backup fields
2. The construction — the rigid integral, the p-adic Abel–Jacobi map, the two
   classes, what is THEOREM vs CONJECTURAL, sources
3. The implementation plan — module by module, precision, cost, hardest part
4. The kill test — exact comparison, tolerance, decision rule
5. Honesty — what this can never prove, and exactly where the conjecture enters
6. References — verified against arXiv/published sources, 2026-09-13

One-paragraph summary. 389a1 has prime conductor N = 389 with split
multiplicative reduction, so the only implementable Darmon-type construction
in the literature (Bertolini–Darmon and successors) forces the auxiliary
inert prime to be p = 389 itself, and the genus-character sign condition
forces K = Q(sqrt(D1*D2)) with D1, D2 < 0 both fundamental. For such K the
conjectural Stark–Heegner points are p-adic objects mapped into E(Q_389)
via the Tate uniformisation; rationality of genus combinations over the
genus field is a THEOREM (Bertolini–Darmon 2007, hypotheses removed by Mok
2018) but NO Gross–Zagier-type formula for complex Néron–Tate heights
exists anywhere in the literature — the rank-2 statement this document
tests is OUR extension, and section 5 says precisely what a miss kills.

## 1. The choice of K

**Root-number algebra (all of it THEOREM; formulas standard, cf. Rohrlich).**
For (d,N) = 1 the quadratic twist satisfies `w(E ⊗ χ_d) = w(E)·χ_d(−N)`, so

```
w(E/K) = w(E)·w(E ⊗ χ_K) = χ_K(−N).
```

389a1 has N = 389 prime, split multiplicative reduction at 389 (a_389 = +1),
w(E) = +1, rank 2. The Bertolini–Darmon setting (section 2) for prime
conductor forces p = 389 to be **inert** in K: `χ_K(389) = −1`. With K real,
`χ_K(−1) = +1`, hence

```
w(E/K) = χ_K(−389) = χ_K(−1)·χ_K(389) = −1  (always, in this setting).
```

Since `L(E/K,s) = L(E,s)·L(E ⊗ χ_K,s)` and `ord L(E,1) = 2`, we get
`ord L(E/K,1) = 2 + ord L(E ⊗ χ_K,1) ≥ 3`, odd. Two consequences:

- The Heegner-point-style statement "trace point non-torsion iff
  `ord L(E/K,1) = 1`" predicts the class-group *trace* is torsion here. The
  rank-2 information cannot come from the trace; it must come from the
  individual class-group-graded points (section 2).
- `rank (E ⊗ χ_K)(Q)` is odd, generically 1. Verify numerically before use.

**The genus sign condition.** Bertolini–Darmon's rationality theorem needs a
genus character `χ = χ₁χ₂` of K (Kronecker characters of fundamental
discriminants D₁, D₂ with `d = D₁·D₂`, `(D₁,D₂) = 1`, `(d,389) = 1`) with

```
χ₁(−M) = −w_M,   where M = N/p = 1, w_M = +1  ⟹  χ₁(−1) = −1  ⟹  D₁ < 0.
```

Both D_i < 0 keeps `d = D₁D₂ > 0` (K real) while the genus field
`H_χ = Q(√D₁, √D₂)` is CM biquadratic. Then
`w(L(E/K,χ,·)) = χ₁(−389)χ₂(−389) = (−1)(−1)·χ_d(389) = −1`:
`L′(E/K,χ,1)` is the central object, and
`L(E/K,χ,s) = L(E ⊗ χ₁,s)·L(E ⊗ χ₂,s)` makes it a *product of one odd-rank
and one even-rank twist factor* — the mechanism section 2 exploits.

**Primary recommendation: K = Q(√21)** (d = 21 = (−3)·(−7)).
- χ_{−3}(389): 389 ≡ 2 (mod 3) ⟹ −1. χ_{−7}(389): 389 ≡ 1 (mod 4),
  (7/389) = (389/7) = (4/7) = 1 ⟹ +1. Product χ_{21}(389) = −1: **inert** ✓.
- w(E ⊗ χ_{−3}) = χ_{−3}(−389) = (−1)·(−1) = +1 (even rank, expect 0);
  w(E ⊗ χ_{−7}) = χ_{−7}(−389) = (−1)·(+1) = −1 (odd rank, expect 1);
  w(E ⊗ χ_{21}) = −1 (expect rank 1). All three predictions must be *checked*
  numerically (step K-verify below), not assumed.
- h⁺(Q(√21)) = 2 (narrow; no unit of norm −1), ordinary class number 1:
  the narrow class group itself supplies the two classes of the mission.

**Backups** (same verification protocol): Q(√33) = (−3)(−11), Q(√57) =
(−3)(−19). Selection among {21, 33, 57} is by K-verify output: prefer
analytic rank of E ⊗ χ_d exactly 1 and smallest d.

**Step K-verify (script, no new theory).** For d ∈ {21, 33, 57}:
(i) `bsdlab.twists.kronecker(d, 389) == −1`; (ii) build `E_d =
twists.twist(E, d)`, check `lseries.root_number` and `root_number_numeric`
agree = −1; (iii) `lseries.analytic_rank(E_d)` = 1 and for the genus factors
E_{D_i}: ranks 0 and 1 as above; (iv) narrow class number of Q(√d) by
counting reduced indefinite binary quadratic forms of discriminant d
(≈ 40 lines, cycle-detection on continued-fraction-style reductions —
classical, algorithm in e.g. Buchmann–Vollmer). All four checks use
machinery that already exists in this repo or is 40 lines of integer
arithmetic.

## 2. The construction

### 2.1 Setting (as in Bertolini–Darmon; Mok 2018 states it exactly)

E/Q of conductor N = p·M, p odd prime, **split multiplicative** at p; K real
quadratic of discriminant d = D₁·D₂ (D_i < 0 fundamental, pairwise coprime,
coprime to N); **p inert in K, every q | M split in K**; genus character
χ = χ₁χ₂ with the sign condition of section 1. For us: p = 389, M = 1
(the "M split" condition is vacuous — this is exactly the case Mok 2018
covers after removing Bertolini–Darmon's extra hypothesis that E have
multiplicative reduction at some second prime).

### 2.2 The p-adic Abel–Jacobi recipe (CONJECTURAL object, implementable)

1. **Optimal embeddings.** ψ : R = Z[√d-ish order] → M₂(Z), optimal at the
   Γ-level used; ψ-classes up to Γ-conjugacy ↔ narrow ideal classes of K.
   For d = 21: two classes, ψ₁, ψ₂. Each ψ(K^×) ⊂ GL₂(Q) has two fixed
   points α_ψ, β_ψ ∈ P¹(Q_p) on the boundary of the p-adic upper half
   plane H_p (p inert ⟹ the torus is anisotropic mod p, fixed points are
   Galois-conjugate over Q_p(√d)/Q_p). [Classical; Darmon 2001, §4.]
2. **The rigid cocycle.** Φ_f : {paths in H_p} → C_p, the *rigid Eichler–
   Shimura cocycle* attached to the newform f of E: the unique
   Γ-equivariant rigid-analytic function with Φ_f(r) = the p-adic avatar
   of the modular symbol 2πi∫_{i∞}^{r} 2πi f dz for r ∈ P¹(Q).
   Computable by overconvergent modular symbols (Greenberg 2006-style;
   made algorithmic by Guitart–Masdeu 2013, arXiv:1307.2556, whose public
   code is the reference implementation). This is the one genuinely
   p-adic subroutine — see §3.
3. **The period.** J_ψ = Φ_f evaluated on the "geodesic" (rigid path
   class) from α_ψ to β_ψ ∈ C_p^× (multiplicative presentation; the
   cocycle takes values in C_p^× modulo the period lattice of the
   Tate curve). [CONJECTURAL well-definedness up to the predicted
   lattice — this is the global-compatibility input, §5.]
4. **Tate uniformisation.** E has split mult reduction at 389, so
   E(Q_p) ≅ Q_p^×/q_E^Z with Tate parameter q_E ∈ 389·Z_389. Recover q_E
   to p-adic precision by Newton-solving j(q_E) = j_E where
   j(q) = q^{-1} + 744 + 196884q + … (j_E ∈ Z is exact input).
   Then P_ψ := η_Ta(J_ψ) ∈ E(Q_p) ⊂ E(Q_p̄). [THEOREM as a *local*
   recipe; global meaning of P_ψ is the conjecture.]

### 2.3 The two classes and the statements they should satisfy

The two narrow ideal classes of Cl⁺(Q(√21)) give P_{ψ₁}, P_{ψ₂} ∈ E(Q_p).

| statement | status | source |
|---|---|---|
| genus combination P_χ = χ-weighted Σ P_ψ is, up to rational scalar, a global point in E(H_χ)⊗Q on the χ-eigencomponent | **THEOREM** | Bertolini–Darmon (as cited by Mok, arXiv:1801.09285, Thm 1.1 = BD Thm 4.3); hypothesis removed by Mok Thm 1.2 |
| P_{ψ} themselves lie in E(H_c⁺) with Galois-equivariance | CONJECTURE | Darmon's Stark–Heegner conjecture; numerical evidence: Guitart–Masdeu 1307.2556, Damm-Johnsen 2301.08977 Tables 1–2 |
| P_χ non-torsion ⟺ L′(E/K,χ,1) ≠ 0 | CONJECTURE | BD/Mok intro; L(E/K,χ,s) = L(E⊗χ₁,s)·L(E⊗χ₂,s), sign −1 (section 1) |
| p-adic GZ: weight-derivative of twisted triple product p-adic L = log_E(P_ϕ)·L_p(E,s)·(explicit) | **THEOREM** | Hsieh–Yamana, arXiv:2002.11858, Thm A |
| *complex Néron–Tate height formula for P_ψ* (any GZ-type identity for ĥ) | **DOES NOT EXIST** in the literature | §5 |
| our rank-2 extension under test (below) | **UNVERIFIED — ours** | this document |

**Our extension (the thing actually tested; UNVERIFIED).** Since L(E/K,χ,·)
= L(E⊗χ₁)·L(E⊗χ₂) with ranks (predicted) 0 and 1, the genus machinery
produce points on the twists, not on E(Q). The only route to E(Q)-points is
the individual P_ψ: conjecturally in E(H_c⁺), Galois-equivariant, so the
norm Nm(P_{ψ_i}) = Σ_{σ∈Gal} σ P_{ψ_i} ∈ E(K), and its (P + P̄)/2 component
lands in E(Q)⊗Q. **Hypothesis K2 (unverified, falsifiable):** for the two
narrow classes ψ₁, ψ₂, the points π_i = (Nm P_{ψ_i} + conjugate)/2 ∈ E(Q)⊗Q
are Q-linearly independent and det ⟨π₁,π₂⟩_NT = c_d · Reg(E/Q) for an
explicitly *predictable* rational c_d (predicted form: c_d = (index)²;
see §4). No literature supports or refutes K2; that is precisely why it is
worth one cheap numerical shot.

## 3. The implementation plan

New module `bsdlab/darmon.py` (+ `bsdlab/qadic.py` for fixed-precision
p-adics: class Qp as (mantissa: int, valuation) mod p^n — plain integer
arithmetic, ~120 lines; mpmath is *not* used on this path). Reused as-is:
`periods.period_lattice`, `mordellweil` heights, `oracle.test_generators`,
`lseries` + `twists` for step K-verify, `certify` error bars.

- **M0. Classical modular symbol engine** (~150 lines, mpmath):
  {r,s}_f := 2πi ∫_{r}^{s} 2πi f(z) dz for r,s ∈ P¹(Q), via term-by-term
  integration of the q-expansion (sympy polynomial arithmetic) between
  cusps, consistent with `periods`. Cross-check: {0,∞}_f must recover
  L(E,1)/Ω-ish rational periods to dps-40. Doubles as the L6(i) twist-census
  feed. *Not p-adic. Do this first — it is the classical avatar of Φ_f and
  the calibration for everything after.*
- **M1. Step K-verify** (~60 lines): section 1 checks on d ∈ {21,33,57};
  narrow class group by cycle-counting on reduced indefinite forms (~40
  lines, Buchmann–Vollmer Ch. 6–7). Output picks the primary K.
- **M2. Qp class + Tate parameter** (~80 lines): Hensel square roots in
  Q_p(√d) (needed for the fixed points α_ψ, β_ψ ∈ P¹(Q_p̄)), and q_E by
  Newton on j(q) = j_E, both to n = 60 p-adic digits (q_E ∈ 389Z_389).
- **M3. Overconvergent modular symbols — THE HARDEST SUBROUTINE** (~700
  lines, weeks of care): distributions on Z_p^× of tame level 1 for p = 389
  to overconvergence, the U_p operator as an integer matrix mod p^n
  (entries from Gauss-sum/M moments — the Greenberg 2006 / Guitart–Masdeu
  1307.2556 algorithm; their public `darmonpoints` code is the structural
  reference, not a dependency). Then Φ_f on the path α_ψ→β_ψ, J_{ψ_i} in
  C_p^×, P_{ψ_i} = η_Ta(J_{ψ_i}) ∈ E(Q_389). Precision: n = 60 digits.
  **This step cannot be avoided by an archimedean substitute**: the whole
  content of the conjecture is the p-adic integration (§5). Wall-clock:
  U_p on a ~O(p)-dimensional space to power ~n² ≈ hours in CPython;
  acceptable.
- **M4. Recognition + kill test** (~120 lines): match P_{ψ_i} ∈ E(Q_p)
  against the Mordell–Weil lattice: solve P ≈ m·P₁ + n·P₂ (mod p^n) for
  integers m,n by p-adic lattice reduction (LLL on 3×3, exact integer
  mod p^n linear algebra); then §4.

Estimated totals: ~1100 new lines, of which M3 is two-thirds; dps-60
archimedean, n-60 p-adic; wall-clock dominated by M3 (hours, not days,
once correct). Dependencies: none beyond sympy/mpmath/numpy + an LLL
(~80 lines, standard integer LLL) — no Sage, no PARI, consistent with the
repo's independence rule.

## 4. The kill test

Target value (FINDINGS.md, 10+ digits, L1-certified):
`Reg(389a1) = L''(E,1)/2! · |tors|² / (∏c_p · Ω) = 0.1524601779836…`
with Ω = 4.9804251217…, ∏c_p = 1, tors = 1, L''/2! = 0.7593165002….

The test has three escalating gates; each is decisive on its own.

- **Gate A (recognition).** From M4, attempt P_{ψ_i} ≡ m_i·P₁ + n_i·P₂
  (mod 389^60) where P₁, P₂ are our known generators of E(Q). No integer
  solution ⟹ the candidate points do not lie in the predicted lattice
  *even allowing the full known Mordell–Weil group* ⟹ **K2 dead**, and the
  failure is informative: it localizes to either the globality conjecture
  or (more likely) our normalization of J_ψ, which must then be audited
  against the one *proven* case available — the genus combination, whose
  rationality is a theorem. Two failed normalization audits kill the
  experiment.
- **Gate B (structure).** With (m_i, n_i) recovered, predict
  `c_d = det[[m₁,m₂],[n₁,n₂]]²` and check
  `det ⟨π₁,π₂⟩_NT = c_d · Reg(E/Q)` to the L1 error budget (`certify`
  intervals, relative width < 1e-9 achievable at dps 60). This must hold
  *identically* — it is linear algebra, not conjecture — once Gate A
  passes; failure means a bug, not a kill.
- **Gate C (the oracle).** `oracle.test_generators(E, [π₁, π₂])`:
  the height-pairing determinant of the *independently constructed*
  points vs Reg. Decision rule: after convergence in the p-adic digit
  count (n = 40, 50, 60 stable) and the overconvergence parameter,
  **relative deviation > 1% ⟹ the construction as conjectured (K2 + the
  globality input it needs) is dead for this (E, K)**; agreement to 1e-3
  or better ⟹ corroboration — the strongest output this lab can produce
  (PLAN.md §5.2(b): "the biggest event this lab could produce").

## 5. Honesty section

1. **This computation can never prove BSD, nor any part of it.** A Gate-C
   match shows a *construction* reproduces points it predicts; it produces
   zero theorems. The missing GZ-type height formula would still have to
   be *proven*, and no such statement exists to prove (§2 table, row 5).
2. **Exactly where the conjecture enters the algorithm.** M3's Φ_f is a
   rigid cocycle defined up to periods, and the pairing J_ψ ∈ C_p^×/q^Z is
   well-defined *only conjecturally* (Darmon's "global compatibility of
   the p-adic integration"): different path representatives for the
   "geodesic" α_ψ→β_ψ differ by lattice elements, and the conjecture
   asserts these lie in the *predicted* lattice. The proven theorem covers
   only the genus-character combination (Bertolini–Darmon via Mok). So a
   Gate-A failure does not even falsify Darmon's conjecture — the
   individual P_ψ may be global over H_c⁺ while the E(Q)-projection K2
   asserts simply is not there (E(Q) ⊗ Q need not meet the H_c⁺-trace
   image nontrivially at all: no theory predicts it does).
3. **What a kill kills.** Precisely: hypothesis K2 of §2.3 for this curve
   and field class. That is still worth having — K2 is this document's
   best-justified candidate for the "two classes" the rank-2 barrier
   needs, and the falsification costs hours.
4. **UNVERIFIED items** (marked in text): the exact bibliographic form of
   the original Bertolini–Darmon reference (accessed only through Mok's
   citation); Darmon's 2001 original paper's exact hypotheses (our §2.1
   is as restated by Mok and Hsieh–Yamana); no numerical constant c_d
   prediction beyond "index²" is claimed.

## 6. References

Verified against arXiv (fetched 2026-09-13; titles/authors/ids as listed
there, claims above checked against abstracts/intros, not memory):

- C. P. Mok, *On a theorem of Bertolini–Darmon about rationality of
  Stark–Heegner points over genus fields*, arXiv:1801.09285. States BD Thm
  4.3 exactly (§1 of this doc); removes the "multiplicative at a second
  prime" hypothesis — the fact that prime conductor 389 is admissible
  rests on this.
- X. Guitart, M. Masdeu, *Overconvergent cohomology and quaternionic
  Darmon points*, arXiv:1307.2556. The M3 algorithm; public code.
- X. Guitart, M. Masdeu, *Computation of ATR Darmon points on
  non-geometrically modular elliptic curves*, arXiv:1204.6680. Template
  for recognition (Gate-A style `J_τ =? 4z` at 5·10⁻¹¹) — but note: ATR
  setting (E over real quadratic F, archimedean double integrals), NOT
  our setting; cited for method, not mathematics.
- X. Guitart, M. Masdeu, M. H. Şengün, *Darmon points on elliptic curves
  over number fields of arbitrary signature*, arXiv:1404.6650. The
  general framework: choice of non-split place v, cocycle Δ_ψ, additive
  pairing, uniformization at v; Conjectures 1.2–1.3 = our §2.2 items 3–4.
- M.-L. Hsieh, S. Yamana, *Restriction of Eisenstein series and
  Stark–Heegner points*, arXiv:2002.11858. Theorem A: p-adic GZ formula
  in the p|N, p inert, N/p split setting (matches §2.1). This is the
  strongest *proven* statement in the neighborhood of this experiment.
- H. Damm-Johnsen, *Modular algorithms for Gross–Stark units and
  Stark–Heegner points*, arXiv:2301.08977. Recent explicit computations
  (tables D < 100); confirms field-of-definition conjectures numerically;
  no height formula — consistent with our §2 table row 5.
- M. Longo, S. Vigni, arXiv:1105.3721; M. Longo, V. Rotger, S. Vigni,
  arXiv:1004.3424 (GZ-type formulas for *quaternionic* Darmon points on
  Shimura curves — requires N⁺ > 1, not our case; cited for context).
- M. Fornea, L. Gehrmann, *Plectic Stark–Heegner points*,
  arXiv:2104.12575 (+ arXiv:2203.15998): the literature's closest object
  to "two classes / higher rank"; conjectural, p-adic; listed as the
  fallback construction if K2 dies.
- H. Darmon, *Rational Points on Modular Elliptic Curves*, CBMS Reg.
  Conf. Ser. 101 (2004) — standard reference for the program. [Book;
  not fetched today — treat section-level claims as UNVERIFIED.]
- Bertolini–Darmon original (genus-field rationality), accessed only via
  Mok's citation [BD1]. **Exact title/venue UNVERIFIED** — a worker
  implementing this must read it before M3.
