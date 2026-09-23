# BSD-Lab research-agent campaign — operational plan

Owner: GLM-5.2 orchestrator + GLM-5.2 headless workers (all agents on glm-5.2,
z.ai plan). Strategic plan: [`PLAN.md`](PLAN.md) (routes A/B/C, lab program
L1–L6). This document is the **execution** plan: who runs what, in what order,
with what acceptance. Started 2026-09-13.

## 1. Architecture

- **Orchestrator** (this session, GLM-5.2): scopes, dispatches, verifies,
  integrates, commits. Does no heavy work unless the two-failure rule fires.
- **Workers**: headless GLM-5.2 via the `glm` launcher. 8192-token output cap
  → every task file bakes in the discipline: skeleton first, one function per
  edit (≤120 lines), ≤25-line final report. **Concurrency 2** (z.ai limits).
- Workers are stateless: `scratch/p_glm_*.md` task files carry the full spec,
  file whitelist, acceptance values, and report format.
- **Verification at boundaries**: the orchestrator re-runs acceptance checks,
  reads the diff, spot-checks one number. Two failed delegations of the same
  task → orchestrator implements it directly.

## 2. Waves

**W1 — engineering foundations (dispatched 2026-09-13):**
- **W1a = L2 descent via 2-isogeny** → `bsdlab/descent.py`. Turns the rank
  lower bound into a two-sided bound; computes Sel^φ images and Sha[φ].
  Full spec: `scratch/prompt_descent.txt`.
- **W1b = L1 interval certificates** → `bsdlab/certify.py`. Rigorous error
  budget (L-tail via Hasse |a_n|≤2√n, AGM last-iterate gap, regulator
  dps-50/80 diff ×10 safety) → certified `#Sha_an` intervals for the
  14-curve reference set.

**W2 — instrument + local theory:**
- **W2a = L3 regulator oracle**: `test_generators(E, [P1,P2])` → height-pairing
  determinant vs `L^(r)(1)/r!·|tors|²/(∏c_p·Ω)` with W1b-grade error bars.
  This is the falsifier every Route-B candidate faces.
- **W2b = L5 root numbers at additive primes**: derive the Halberstadt/Kraus
  tables in-house from our exact Tate data; restore the two-route sign
  cross-check for the 32a3 class.

**W3 — research wave 1 (creative maths):**
- **W3a = L6(i) twist census** on 389a1/571b1/5077a1: rank distribution of
  E⊗χ_d for fundamental |d| ≤ 200; quantify the PLAN §5.2(a) trace
  obstruction with data.
- **W3b = Route-B construction spec**: `docs/ROUTE_B_DARMON.md` — an
  implementable numerical recipe for Stark–Heegner/Darmon points on 389a1
  over a real quadratic field, ending in the height-pairing kill test
  against Reg(389a1) = 0.1524601779….

**W4 — research wave 2:**
- **W4a = L6(ii) triple-product partner hunt** over the 39,790-curve base.
- **W4b = Route-A status** (web): `docs/ROUTE_A_STATUS.md` — what is actually
  proven at p = 2, 3 today, cited; where A2's gap sits.
- **W4c = constructions agent**: generate ≥5 NEW falsifiable rank-2
  constructions, each scored by the W2a oracle; survivors earn theory push.

**W5 — scale: L4**, 14 → ~500 certified curves, prioritizing `#Sha_an > 1`.

## 3. Decision gates

- **G1 (after W1)**: algebraic rank known independently of L for reference
  curves with 2-torsion; certified #Sha intervals at rank ≤ 1. Update
  FINDINGS.md; commit.
- **G2 (after W2)**: the lab can kill any proposed rank-2 construction in
  minutes — the campaign's strategic weapon.
- **G3 (after W3/W4)**: a construction survives its kill test → all-hands
  theory push on it; none does → PLAN §5.2 pruned by experiment; publish the
  negative results + census data as the deliverable.
- Rungs R1/R2/R3 as in PLAN.md §8. Risk stance unchanged from PLAN.md §9:
  guaranteed deliverables are L1–L5 + the falsifier role; R3 itself remains
  a Millennium-grade long shot, attacked honestly.

## 4. Log

- 2026-09-13 — campaign opened. W1a + W1b dispatched to GLM-5.2 workers.
- 2026-09-13 — **twist bug fixed (orchestrator).** `twists.twist` applied
  the d-scaling to raw `E.a2/a4/a6` without first passing through
  `a13_zero_model`: every twist had j = 21952/9 (wrong; 389a1 d=1 came out
  rank 0!). Root cause of the conductor crashes at d ≡ 1 mod 4 with odd
  p|d in `_tate_step` (p ≥ 5) also fixed by routing through the global
  short model y² = x³ − 27c4·x − 54c6 before the star-type step. Verified:
  N(389a1⊗χ_d) = 389·∏p² for d = 5, −7, 13, −155, 184, −4 + regressions
  11/37/32/14. Commit 165ac15. Selftest 140/0.
- 2026-09-13 — **W2b + W3b wave dispatched** (root numbers at additive
  primes; Darmon/Stark–Heegner spec). Both GLM slots busy. NOTE: wave
  launched before the twist fix landed — W3b's numbers must be re-verified
  against the fixed `twists.py` at boundary.
- 2026-09-13 — **W3a census running** (twist census, 369 fundamental
  |d| ≤ 200 twists of 389a1/571b1/5077a1; analytic rank + root number
  each; resumable). Also W5 driver (`certify_scale.py`) and the E4 task
  ladder staged: W4p1 (padic + MST p-heights) → W4p2 (exact rational
  modular symbols) → W4p3 (MTT p-adic L). E4 driver
  (`scratch/e4_beta_table.py`) written against the W4 specs — the
  𝔅_p = L_p^(r)/(r!·E_p·Reg_p) table testing GRS 𝔅_p = 1.
- 2026-09-13 — **W2b landed** (`bsdlab/rootnumbers.py`; GLM-5.2 hit quota,
  finished on fallback). Orchestrator verification: 32a3 table route = +1 =
  numeric route; 14 reference curves → 13 two-route agreements, 1 honestly
  blank (undecoded additive p=2 row → `P2_TABLE_MISSING`, no guess), 0
  disagreements. Selftest 140/0. Not yet wired into `lseries` — the census
  keeps the numeric fallback.
- 2026-09-13 — **W3b landed** (`docs/ROUTE_B_DARMON.md`). Recommends
  K = Q(√21) (backups √33, √57), p = 389 inert, rigid p-adic (Bertolini–
  Darmon/Mok) setting. Key honest finding: **no Gross–Zagier-type complex
  height formula exists for Stark–Heegner points**; the rank-2 test (K2) is
  our own unverified extension, worker-estimated 10–15% to match. M3
  (overconvergent modular symbols, ~700 lines) is the bottleneck.
- 2026-09-13 — **Spec defect caught in W4p1 (orchestrator review).** MST
  p-adic heights depend on the p-adic E2(E, ω); a wrong E2 shifts h_p by
  c·log_E(P)² — still quadratic, so the bilinearity acceptance could never
  detect it and every 𝔅_p would be silently wrong. Split into W4p1a (padic +
  Kedlaya Frobenius + E2, gated by charpoly = x² − a_p x + p and the Serre
  congruence E_{p+1} ≡ E_2 mod p) and W4p1b (heights). Wave W4p1a + W4p2
  dispatched to GLM-5.2.
- 2026-09-13 — **W4p2 landed: exact modular symbols (`bsdlab/modsym.py`)**,
  after an orchestrator fix. The worker passed all 16 anchor values exactly
  (Merel Heilbronn matrices; T_l verified by exact rational asserts). But an
  **out-of-sample gate** (`tests_research/oos_modsym.py`: modsym vs 40-digit
  numeric L(E_d,1)/Ω(E_d) on twists of 11a1/37a1/389a1 the worker never saw)
  found **10/60 wrong by factors of 2 or 4**. Cause: the normaliser C had been
  fitted to the anchors, but the Hecke eigenvector's scale is arbitrary per
  curve *and* sign. Fix: `period_scales` fixes the scale once per
  (curve, sign) against E's own Ω⁺ (least real period) / Ω⁻ (least imaginary
  period), and `algebraic_twisted_l` returns exact L(E,χ_d,1)·√|d|/Ω^±. The
  twist's own-period value multiplies by a lattice index in {1/2, 1, 2}
  (identified to 1e-12). Result: **68/68** (≈62 independent of the
  calibration twists), anchors 16/16. Lesson for every future worker task:
  a constant pinned to anchors must be validated off-anchor.
- 2026-09-13 — **W4p1a attempt 1 died on the 8192 output cap** after delivering
  `bsdlab/padic.py` (A1 passes: ring identities mod p^20, Hensel sqrt,
  log∘exp). Resumed for Frobenius/E2 only.
- 2026-09-13 — **W4p3 spec rewritten (orchestrator).** v1 had the MTT measure
  wrong (α⁻ⁿλ·C_p instead of the p-stabilised α⁻ⁿ[a/pⁿ]⁺ − α⁻⁽ⁿ⁺¹⁾[a/pⁿ⁻¹]⁺),
  demanded L_p''(0) be a unit (p-adic BSD only predicts nonzero), and dropped
  the Euler factor (1−1/α)² from the interpolation check. v2 follows
  Stein–Wuthrich 2013 §3 and gates on theorems: distribution relation,
  trivial- and quadratic-character interpolation.
- 2026-09-13 — **W4 p-adic ladder complete, externally anchored.** GLM-5.2 hit
  quota twice and the TPU fallback tunnel is down, so the orchestrator
  finished the ladder. Every stage is checked against published Sage/MAGMA
  doctest values, not just internal consistency:
  - `bsdlab/frobenius.py` (Kedlaya Frobenius + Katz E2; worker died twice).
    Charpoly = x² − a_p x + p on 6 curves × 6 primes. E2(37a1, 5) equals
    Sage/MAGMA to **5^20**. Katz congruence E_{p+1} ≡ E2·E_{p−1} mod p holds.
    (A first version of that gate omitted the Hasse invariant E_{p−1} and
    wrongly failed 20 correct values.)
  - `bsdlab/padicl.py` (MTT p-adic L). The worker draft had four defects, all
    fixed. The serious one: `unit_root` used a_p − D for the other square
    root, giving a wrong α for about half of all (E, p), with
    sign-alternating Riemann sums for 389a1 at p = 5. Also a PNum addition
    that dropped a summand, [a/p^(n+1)] in place of [a/p^(n−1)], and w1 in
    place of the Néron period. Now matches Sage's L_5(389a, T) digit for
    digit (series(2) and series(3)) and Sage's 37a measure to 5^9, and
    satisfies MTT interpolation for 11a1 at p = 7, 13.
  - `bsdlab/pheights.py` (MST heights). h_5(37a1) and **Reg_5(5077a1)**
    (rank 3, anomalous) equal Sage to 5^10 on the first run. Side finding:
    the commonly quoted 5077a1 points (−3,0), (−2,3), (−1,3) span an
    index-2 sublattice (classical regulator 4×); a true basis is
    (0,−3), (−1,−4), (1,−1).
- 2026-09-13 — **E4 pilot:** 𝔅_p = 1 for 37a1 (p = 5, 7 mod p^4; p = 11 mod
  11^3) and 389a1 (mod 5^4, 7^3, 11^2, 13). As expected, this reproduces
  Stein–Wuthrich 2013, now from an independent Sage-free implementation.
  The full runs (p ≤ 60; 5077a1 p ≤ 40) write `data/e4_beta_table*.json`.
- 2026-09-23 — **Route B exhaustive pass (orchestrator)**
  → `docs/research/06_route_b_exhaustive.md` + `docs/GENERAL_SOLUTION_PLAN.md`.
  - **First-Derivative No-Go Lemma** proved: any construction with a
    `⟨c,c⟩ ∝ L'(E/K,χ,1)` height formula contributes nothing to `E(Q)`
    once `ord L ≥ 2`.
  - **RB1** (`tests_research/rb1_heegner_no_go.py`): Heegner points are
    torsion on 21/21 discriminants for `389a1` and 22/22 for `5077a1`.
    This reproduces the Buhler–Gross–Zagier 1985 mechanism Sage-free.
    Control `37a1`: `L_alg(E^D,1) = 2k²` exactly on 22/22 (Gross–Zagier).
    A ramified `D = −111` leaked through the first filter and was caught by
    the `2k²` law.
  - **Verdicts:**
    - Heegner/twists: killed.
    - Stark–Heegner/Darmon: killed by the same lemma. W3b's M3 stage is
      cancelled.
    - Beilinson–Flach: dormant.
    - Kolyvagin systems, generalised Kato / diagonal classes, derived
      heights: alive, per prime.
  - **Remaining content:** MT1 (uniformity in `p`), MT2 (complex ↔
    `p`-adic order), MT3 (higher Gross–Zagier over `Q`), MT4 (Route A).
  - **Other runs:**
    - E4 extended: `37a1` complete to `p ≤ 59`, including the 4
      previously unresolved primes; `389a1` and `571b1` complete to
      `p ≤ 59`. That is 42/42 cells with `𝔅_p = 1`.
    - Twist census: `389a1` complete (123 twists, rank 0/1/2 = 55/60/8,
      no parity violations); `571b1` complete (rank 0/1/2 = 46/61/16, no
      parity violations); `5077a1` 34/123. The run was stopped by the host
      for low memory. It is resumable with `python scratch/census_w3a.py 200`.
- 2026-09-13 — **W1a landed.** Worker died twice on the 8192-token output
  cap (third attempt produced 276 of ~280 lines before dying); two-failure
  rule fired, orchestrator finished and certified it. Acceptance 7/7 PASS:
  32a rank certified 0; n=5 and n=7 rank certified 1 (algebraic rank now
  independent of L); 14a1/15a1/681b1 honest two-sided bounds. Selftest
  140/0. Frontier curves 389a1/571b1/5077a1 have no rational 2-torsion —
  2-isogeny descent does not apply to them (honest None); general
  2-descent remains future work. W2a dispatched into the freed slot.
