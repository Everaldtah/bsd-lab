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
- 2026-09-13 — **W1a landed.** Worker died twice on the 8192-token output
  cap (third attempt produced 276 of ~280 lines before dying); two-failure
  rule fired, orchestrator finished and certified it. Acceptance 7/7 PASS:
  32a rank certified 0; n=5 and n=7 rank certified 1 (algebraic rank now
  independent of L); 14a1/15a1/681b1 honest two-sided bounds. Selftest
  140/0. Frontier curves 389a1/571b1/5077a1 have no rational 2-torsion —
  2-isogeny descent does not apply to them (honest None); general
  2-descent remains future work. W2a dispatched into the freed slot.
