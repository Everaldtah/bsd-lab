# Route B, exhaustively: every known rank-2 mechanism, tested

*Research note, 2026-09-23. Numerics are from this lab; the no-go lemma is
proved inline; everything else is cited. No open case of BSD is proved here.*

## 0. Summary

We went through every mechanism in the literature that has been proposed
as a source of rank-≥2 information for `E/Q`. They fall into two groups:

- **Point constructions:** Heegner points, twists of Heegner points,
  Stark–Heegner/Darmon points, ATR points. These are governed by a
  **first** derivative `L'(E/K, χ, 1)`. By the No-Go Lemma (§1) they can
  never produce a non-torsion point of `E(Q)` once `ord_{s=1} L(E,s) ≥ 2`.
  Experiment RB1 confirms this for 43 of 43 Heegner points on `389a1`
  (rank 2) and `5077a1` (rank 3), to 31 digits. The same program
  reproduces the known non-torsion Heegner points on `37a1` (rank 1) as
  exact integer multiples of the generator.
- **Class constructions:** Kolyvagin systems, generalised Kato /
  diagonal-cycle classes, Beilinson–Flach elements, and derived `p`-adic
  heights. These produce **Selmer classes, one prime at a time**. Some do
  see rank 2: Castella–Hsieh prove `κ_p ≠ 0 ⇒ dim Sel(Q, V_pE) = 2`. None
  are uniform in `p`, and none connect to the *complex* `L''(E,1)`.

Route B therefore reduces to three precisely stated missing theorems,
MT1–MT3 in §4. MT4, the Route A completion, is needed only for the exact formula. They are the whole remaining content of rank equality and
Sha finiteness at rank ≥ 2. Route B does not produce rational points at
rank 2, and it doesn't need to: the generators of `389a1` are already
known. What it has to produce is a **bound**.

## 1. The First-Derivative No-Go Lemma

**Setting.** Let `K` be a quadratic field and `H/K` a ring class field with
`G = Gal(H/K)`. Suppose a construction `Π` produces, for each character
`χ` of `G`, an element

    c_χ ∈ (E(H) ⊗ C)^χ

and satisfies a first-derivative formula

    ⟨c_χ, c_χ⟩ = λ_χ · L'(E/K, χ, 1),    λ_χ ≠ 0

for a height pairing `⟨ , ⟩` that is non-degenerate on `E(H) ⊗ R`.

**Lemma.** If `ord_{s=1} L(E,s) ≥ 2`, then `Π` contributes nothing to
`E(Q) ⊗ Q`. Precisely, for every `K`, `H`, `χ`, the image of `c_χ` under
`Tr_{H/Q}` is zero in `E(Q) ⊗ C`.

**Proof.**

- *Case χ ≠ 1.* `Tr_{H/K}` is the projection onto the trivial isotypic
  component. It kills `(E(H) ⊗ C)^χ`, so the trace is `0`.
- *Case χ = 1.* `L(E/K, s) = L(E,s) · L(E^K, s)`, and `L(E,s)` vanishes to
  order ≥ 2. So `L'(E/K,1) = 0`, hence `⟨c_1, c_1⟩ = 0`. Non-degeneracy
  (Néron–Tate) gives `c_1 = 0` in `E(K) ⊗ R`, so `c_1` is torsion. ∎

**Remarks.**

1. *p-adic versions.* For `p`-adic heights, non-degeneracy is Schneider's
   conjecture, which is open. The `p`-adic lemma therefore only puts `c_1`
   in the kernel of `h_p`. This is weaker, not stronger. The point
   constructions do not escape the lemma by going `p`-adic.
2. *Why rank 1 works and rank 2 does not.* At rank 1,
   `L'(E/K,1) = L'(E,1) · L(E^K,1)`, and `K` can be chosen with
   `L(E^K,1) ≠ 0`. That is the whole of Gross–Zagier + Kolyvagin. At rank 2
   *every* term of the product rule for `L'(E/K,1)` contains a factor that
   vanishes.
3. *What the lemma does not kill.* Kolyvagin's **derivative classes**
   `κ_n`. These are not points of the form `c_χ`; they are images of
   `D_n y_n` under a Kummer map mod `p`, and they live in `Sel_p(E/Q)` or
   `Sel_p(E/K)` with no first-derivative height formula. It also does not
   kill generalised Kato classes, whose governing formula (Castella–Hsieh)
   is a *second-order* leading-term identity for an anticyclotomic
   `p`-adic L-function. This is exactly where §3's live candidates sit.

## 2. Experiment RB1: the lemma, measured

Driver: `tests_research/rb1_heegner_no_go.py`. Data: `data/rb1_heegner.json`.

**Method.**

1. For every fundamental `D < 0`, `|D| ≤ 160`, satisfying the Heegner
   hypothesis, take one Heegner form `[N a', B, C]` per class of
   `Cl(D)`, with `B ≡ β mod 2N`.
2. Sum `φ(τ) = Σ a_n/n · qⁿ` over the orbit. The `a_n` are computed from
   scratch by point counting.
3. Reduce modulo the Néron lattice from the lab's own AGM. The Manin
   constant is 1 for these optimal curves.
4. `z_K ∈ (1/m)Λ` for small `m` is equivalent to `y_K` being torsion.

No Sage, PARI or eclib is used.

| curve | rank | `w` | discriminants | torsion `y_K` | non-torsion `y_K` |
|---|---|---|---|---|---|
| `37a1` (control) | 1 | −1 | 22 | 4 (`D = −95, −104, −107, −139`) | 18, each an exact multiple `k·P` of `P = (0,0)`: `k ∈ {±1, ±2, ±3, ±4, ±6}`, residual ≤ 1e−27 |
| `389a1` | 2 | +1 | 21 | **21** | 0 |
| `5077a1` | 3 | −1 | 22 | **22** | 0 |

**Reading.**

- On `389a1` and `5077a1`, every `z_K` lies in `Λ` itself (`m = 1`,
  residual ≤ 2e−31). The torsion subgroups are trivial, so `y_K = O`,
  exactly as the lemma forces.
- On `5077a1` this is the mechanism Buhler–Gross–Zagier (Math. Comp. 1985)
  used to *prove* `ord L(E,s) = 3`: `y_K` torsion ⇒ `L'(E,1) = 0` by
  Gross–Zagier, and `L'''(E,1) ≠ 0` numerically. RB1 reproduces it outside
  the PARI/Sage lineage.
- On `37a1` the torsion cases `D = −95, −104, −107, −139` are where Gross–Zagier
  predicts `L(E^D,1) = 0`. They are the rank-1 shadow of the same lemma.
  See §2.1 for the independent L-value check.

### 2.1 Cross-check on 37a1: Gross–Zagier in exact arithmetic

For every control discriminant we computed `L_alg = L(E, χ_D, 1)·√|D|/Ω⁻`
**exactly**, as a rational number, from the lab's modular symbols
(`modsym.algebraic_twisted_l`). This is independent of the `q`-series
evaluation. On all 22 discriminants:

    L_alg(37a1 ⊗ χ_D)  =  2·k²,     where y_K = k·P.

This is the Gross–Zagier formula for 37a1, with its constant `2`, recovered
with zero error. In particular, `y_K` is torsion **iff** `L(E^D,1) = 0`,
with 4 of 4 torsion cases and 18 of 18 non-torsion cases matching. The
instrument that reports 43/43 torsion on the rank-2 and rank-3 curves is
therefore validated in both directions on the rank-1 control.

*A filtering note.* A first run admitted `D = −111`, which is ramified at
37. The Heegner hypothesis needs `gcd(D, N) = 1`. At that `D` the law gave
`4` rather than `2k² = 8`, which exposed the filter bug. The driver now
excludes ramified `D`. No rank-2 or rank-3 row was affected, because no
admitted `D` there shares a factor with 389 or 5077.

## 3. The candidates, one by one

The verdicts use these terms:

- **KILLED:** provably cannot supply rank ≥ 2 over `Q`.
- **ALIVE:** supplies genuine rank-2 information, with a precisely known
  gap.
- **DORMANT:** no mechanism connecting it to `L''(E,1)` is known.

### (a) Heegner points and quadratic twists — KILLED as a point source

- **Reason:** §1 and §2.
- **Twist variant:** a Heegner point on an odd-rank twist `E^D` lies in
  the minus part `E(K)^−`. It never traces to `E(Q)`.
- **Twist census:** the lab's census over `|d| ≤ 200` records the rank
  distribution of the twist family (`data/twist_census.json` once complete;
  partial state in `data/twist_census_partial.json`). It locates the twists
  where the family itself stays at rank ≥ 2. By the same lemma, the
  Heegner points attached to those twists are torsion too.

### (a′) Kolyvagin systems — ALIVE, the most concrete route to Sha finiteness at rank 2

- **Mechanism:**
  - `κ_n ∈ H¹(K, E[p])` is built from derivative operators on Heegner
    points of conductor `n`.
  - Kolyvagin's structure theorem (1991) reads the Selmer structure off the
    vanishing orders `ν(n)` of the non-zero `κ_n`.
  - W. Zhang (2014) proved Kolyvagin's conjecture that the system is
    non-trivial, under hypotheses (`p ≥ 5`, large residual image,
    ramification conditions).
  - Jetchev–Lauter–Stein (2009) computed explicit Kolyvagin classes for
    `389a1`.
  - C.-H. Kim's higher Gross–Zagier/Kolyvagin-system refinement
    (arXiv:2203.12161) has no low-rank assumption.
- **What it gives at rank 2:**
  - Take `389a1` and `K` with `rank E^K(Q) = 1`, so `rank E(K) = 3`.
  - Then `Sha(E/K)[p] = 0` iff a class with `ν(n) = 2` (two Kolyvagin
    primes) is non-zero mod `p`.
  - Each such `p` is *computable*.
- **Gap:** the statement must hold for all but finitely many `p` at once.
  This is MT1 below, in its Kolyvagin form.

### (b) Stark–Heegner / Darmon points — KILLED for Route B, even if their conjectures are proved

- **Mechanism:**
  - Points over ring class fields of **real** quadratic `K` come from
    `p`-adic integration on `H_p × H`.
  - They are conjecturally rational.
  - Conjecturally, they are non-trivial exactly when `L'(E/K,ψ,1) ≠ 0`
    (Bertolini–Darmon; see the diagonal-class note arXiv:2207.01310, which
    restates this).
- **Why it dies:**
  - That governing quantity is a first derivative, so the lemma applies
    verbatim.
  - The trivial-`ψ` component is torsion, because
    `L(E/K,s) = L(E,s) · L(E^K,s)` again.
  - The non-trivial `ψ` components vanish under the trace.
- **Consequence:** the `K2` rank-2 height test in `ROUTE_B_DARMON.md`,
  which W3b estimated at a 10–15% chance of matching, is predicted to
  **fail by theorem**. We therefore do **not** build the 700-line
  overconvergent-symbols stage (M3) for Route B. Darmon points remain
  interesting for rank-1-over-`K` questions, which are not ours.

### (c) Beilinson–Flach elements — DORMANT

- **Mechanism:**
  - These are classes in `K₁` of a product of modular curves, with syntomic
    regulators tied to Hida's `p`-adic Rankin L-function
    (Bertolini–Darmon–Rotger, Lei–Loeffler–Zerbes).
  - With an Eisenstein partner they degenerate to `L(f,s) · L(f,χ,s)`.
    Bertolini–Darmon–Rotger II used exactly this for BSD in analytic rank 0
    for Hasse–Weil–Artin L-series.
- **What it gives at rank 2:** nothing known.
  - The known derivative formulas for Beilinson–Flach classes (Gross–Stark
    type) are exceptional-zero *first* derivatives.
  - No formula connects them to `L''(E,1)`.
- **Status:** we record it as dormant, not killed. A partner `g` with
  `L(f ⊗ g, s)` factoring through `L(f,s)²` would be needed. The
  self-product `f ⊗ f` does factor through `Sym²` and `ζ`, not `L(f,s)²`,
  so the obvious partner fails.

### (d) Diagonal cycles / generalised Kato classes — ALIVE, the strongest per-prime tool

- **Mechanism:**
  - Classes come from Gross–Kudla–Schoen diagonal cycles on triple
    products `(f, g, h)`, with `g, h` weight-one theta series
    (Darmon–Rotger; Bertolini–Seveso–Venerucci).
  - For `E` of positive even analytic rank:
    - `κ_p(E) ≠ 0 ⇒ dim Sel(Q, V_pE) = 2` (Castella–Hsieh, Forum Math.
      Sigma 2022).
    - The CM case: Castella, arXiv:2204.09608.
    - The link to Iwasawa main conjectures: Castella, arXiv:2312.01481,
      Thm 5.2.3. It also shows `κ_p ≠ 0` iff `loc_p` is non-zero on the
      Selmer group.
- **Gaps:**
  - (i) `κ_p ≠ 0` is itself conjectural in general (Darmon–Rotger).
  - (ii) `κ_p` is detected by a `p`-adic leading term, not by the complex
    `L''(E,1)`. That is MT2.
  - (iii) `dim_{Q_p}` statements give `Sha[p^∞]` finite at that `p`, not
    `Sha[p] = 0` for almost all `p`. That is MT1.

### (e) Two-variable / anticyclotomic p-adic L-functions and derived heights — ALIVE per prime

- **Mechanism:**
  - The leading term of an anticyclotomic `p`-adic L-function is expressed
    through derived `p`-adic heights and an enhanced regulator
    (Castella–Hsieh).
  - On the cyclotomic side, Banwait (arXiv:2609.08431) gives a criterion
    for rank-2 CM curves: outside an explicit set of primes,
    `L_p''` normalised is a unit **iff** `Reg_p` is a unit and
    `Sha[p^∞] = 0`.
  - This lab's E4 table computes exactly the cyclotomic quantity
    `𝔅_p = L_p^(r)/(r!·ε_p·Reg_p)`, Sage-free. It is `1` at every resolved
    prime for `37a1`, `389a1`, `571b1` and `5077a1`
    (`data/e4_beta_table*.json`).
- **Gap:** each such criterion is per prime, which is MT1. Equating
  `ord L_p` with `ord L` is MT2.

### (f) Higher Gross–Zagier (Yun–Zhang) — the right shape, over the wrong base

- **Mechanism:**
  - Over function fields, `L^(r)` is the self-intersection of
    Heegner–Drinfeld cycles on moduli of shtukas with `r` legs
    (Yun–Zhang, Annals 2017/2019).
  - This has been extended to deeper level (Bieker, arXiv:2607.07531).
- **Gap:**
  - No number-field analogue of "`r` legs" exists. `Spec Z` is not a curve
    over a field, so there is no `Spec Z ×_{F_1} Spec Z` to host the
    cycles.
  - This is MT3. It is the only route in the literature whose *output
    shape* matches `L^(r)(E,1)/r! ↔ det(height matrix)` for every `r`.

## 4. The missing theorems

Every live route ends in one of the following statements. None is known.
Each is stated for a fixed `E/Q` with `ord_{s=1} L(E,s) = r ≥ 2`.

**MT1 — Uniformity in `p`.** `Sha(E/Q)[p] = 0` for all but finitely many
primes `p`. Equivalent forms that would suffice:

- (Kolyvagin form) for all but finitely many `p`, some `κ_n` with
  `ν(n) = r_K − 1` is non-zero mod `p`;
- (Global Rational Shadow, Note 1) `𝔅_p` is a `p`-independent rational
  number;
- (anti-concentration, Note 3) the weighted toric periods are non-zero
  mod `p` for almost all `p`.

MT1 plus per-prime finiteness at the remaining `p` (available from (d)/(e)
under hypotheses, or from descent at `p = 2, 3`) gives **Sha finite for
`E`**. That would be the first rank-2 example ever: rung R2.

**MT2 — Complex ↔ `p`-adic order.** `ord_{s=1} L(E,s) = ord_{T=0} L_p(E,T)`,
and `L''(E,1) ≠ 0 ⇒ κ_p(E) ≠ 0`. Without MT2, the per-prime Selmer
theorems of (d)/(e) cannot be *triggered* by the complex L-function. It is
needed for rank equality in general, not for individual curves, where the
vanishing order is computable up to `r = 3`.

**MT3 — A higher-derivative Gross–Zagier formula over `Q`.**
`L^(r)(E,1)/r! = (explicit local factors) · det⟨P_i, P_j⟩` for classes
`P_i` built from special cycles. This is the number-field analogue of
Yun–Zhang.

- It is the only known *shape* that handles all `r` at once.
- It is also required to even **prove** vanishing orders `≥ 4`. No
  elliptic curve is currently proven to have analytic rank `≥ 4`. Known
  techniques prove `L^(k)(E,1) = 0` only for `k = 0` (modular symbols),
  for `k` odd or even by sign, and for `k = 1` through Gross–Zagier (the
  BGZ trick RB1 reproduces).

**MT4 — Route A completion.** The `p`-part of the exact formula at every
`p`, including `p = 2, 3`, supersingular, additive and reducible cases.
See `PLAN.md` §4. It is needed for the exact formula, not for rank
equality.

**The logical closure.**

- For each individual rank-2 curve: MT1 gives Sha finite.
- For all curves of rank ≤ 2: MT2 + MT1 + a uniform MT1 mechanism give
  rank equality with Sha finite.
- MT3 is needed beyond rank 3, and is the natural source of MT1/MT2 in
  every rank.
- Adding MT4 gives the exact formula: full BSD.

## 5. What we did *not* do, and why

- **Darmon-point numerics (M3, W3b):** superseded by §3(b). The outcome is
  determined by the lemma.
- **Beilinson–Flach partner hunt:** no target formula exists to test
  against, so any hit would be uninterpretable.
- **Generalised Kato classes on `389a1`:** this is the right next
  computation, but it needs Lauder's overconvergent-projection algorithm
  for `p`-adic iterated integrals of weight-one forms. It is specified in
  `GENERAL_SOLUTION_PLAN.md` §5 (L7).

## References

- Buhler, Gross, Zagier, *On the conjecture of Birch and Swinnerton-Dyer for an elliptic curve of rank 3*, Math. Comp. 44 (1985).
- Gross, Zagier, *Heegner points and derivatives of L-series*, Invent. Math. 84 (1986).
- Kolyvagin, *On the structure of Selmer groups*, Math. Ann. 291 (1991).
- W. Zhang, *Selmer groups and the indivisibility of Heegner points*, Camb. J. Math. 2 (2014).
- Jetchev, Lauter, Stein, *Explicit Heegner points: Kolyvagin's conjecture and non-trivial elements in the Shafarevich–Tate group*, J. Number Theory 129 (2009).
- C.-H. Kim, *A higher Gross–Zagier formula and the structure of Selmer groups*, arXiv:2203.12161.
- Castella, Hsieh, *On the nonvanishing of generalised Kato classes for elliptic curves of rank 2*, Forum Math. Sigma (2022), arXiv:1809.09066.
- Castella, *Generalised Kato classes on CM elliptic curves of rank 2*, arXiv:2204.09608.
- Castella, *Nonvanishing of generalised Kato classes and Iwasawa main conjectures*, arXiv:2312.01481.
- Banwait, *Second derivatives of p-adic L-functions and the Shafarevich–Tate group of rank-two CM elliptic curves*, arXiv:2609.08431.
- *Stark–Heegner points and diagonal classes*, arXiv:2207.01310.
- Bertolini, Darmon, Rotger, *Beilinson–Flach elements and Euler systems I, II*, J. Alg. Geom. 24 (2015).
- Bertolini, Seveso, Venerucci, *Diagonal classes and the Bloch–Kato conjecture*, Münster J. Math. 13 (2020); *Reciprocity laws for balanced diagonal classes*, Astérisque 434 (2022).
- Yun, Zhang, *Shtukas and the Taylor expansion of L-functions*, Ann. of Math. 186 (2017); II, 189 (2019).
- Bieker, *A higher Gross–Zagier type formula for moduli of shtukas at deeper level*, arXiv:2607.07531.
