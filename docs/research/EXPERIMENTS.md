# Falsification experiments for the research notes

Specifications, plus results where an experiment has been run (RB1, E4). Every experiment
uses the lab's reference curves (`389a1`, `571b1`, rank 2; `5077a1`, rank 3;
`37a1`, rank 1) and reports certified error bars in the sense of PLAN.md L1.
The priority column orders them by kill value per unit of work.

| id | note | what it tests | priority | new lab capability needed |
|---|---|---|---|---|
| E3 | 4 | geodesic barycenters vs Mordell–Weil (BG) | **1** | numerical `F(z)` via `Γ_0(N)`-reduction; closed-geodesic quadrature; PSLQ (mpmath) |
| E4 | 2 | `ord L_p = ord L` and `p`-independence of `𝔅_p` | **2** | modular symbols → `p`-adic L-function; `p`-adic heights (Mazur–Stein–Tate) |
| E1 | 1 | GRS: `𝔅_p = 1` at every ordinary `p < 200` | 3 | same as E4 (E1 is E4b at larger scale) |
| E2 | 3 | Kolyvagin weighted orbit sums `W_{ℓ,q}` mod `p` | 4 | Brandt matrices (Pizer), Gross points, level-raised mod-`p` eigenforms |
| E5 | 5 | truncated Euler-product volumes | 5 | fast `a_p` to `10^5` (baby-step giant-step) |

## E3 — Barycentric geodesic points (Note 4)

- **Input.** `389a1`, and ten fundamental `D > 0` with narrow class number
  one and `(D/389) = 1`.
- **Step 1 (sanity).** Check that `I(Q_D) = 0` to 40 digits. This validates
  the use of the real-quadratic toric period formula. If it fails, the
  definition of `B(Q)` needs the trace version.
- **Step 2.** Compute `B(Q_D) = (1/ℓ_Q)∫_{C_Q} F ds` to 40 digits. Before
  evaluating `q`-series, move each quadrature node into
  `Im z ≥ √3/(2N)` using `Γ_0(389)` and `w_389`, and add the period
  corrections `ω(γ)`.
- **Step 3.** Run PSLQ on `{B(Q_D), z(P_1), z(P_2), ω_1, ω_2}` with
  coefficient bound `10^6`.
- **Kill criteria.**
  - BG-strong dies if no relation appears for any `D`.
  - BG-span dies if every relation lies on one line.
- **Extra.** Tabulate the traces for `D < 200` and test the BG-modular
  prediction.

## E4 — Realization independence (Note 2)

- **E4a (calibration).** For `37a1`, at every good ordinary `p < 100`, check
  that `L_p'(E)/(E_p·h_p(P)) = 1`.
- **E4b (rank 2).** For `389a1` and `571b1`, at every good ordinary
  `p < 100`: check `L_p'(E) = 0` and `L_p''(E) ≠ 0`, and check
  `L_p''/(2E_p·Reg_p) = 1`.
  - **Kill:** any `p` with `L_p'' = 0`, or a quotient `≠ 1`, beyond
    precision.
- **E4c (formal-kernel test).** Write the first 50 Gross–Zagier
  generating-series coefficients for `37a1` in the formal symbol space `𝓢`,
  realize them at `∞` and at `p = 5, 7`, and confirm that one formal identity
  gives both numerical ones.
  - **Kill:** a needed correction that cannot be written in `𝓢`.

## E1 — Global Rational Shadow (Note 1)

Run E4b over all good ordinary `p < 200`, and add `5077a1` using `L_p'''`
and the `3 × 3` regulator `Reg_p`.

- **Prediction:** `𝔅_p = 1` identically.
- **Deliverable:** the first table of `p`-adic BSD quotients computed outside
  the Sage/PARI lineage, which is the lab's founding thesis.
- **Also record** the primes where `Reg_p` is not a unit. Note 1 §4 predicts
  roughly `Σ 1/p` of them, which quantifies how often unit criteria are
  silent.

## E2 — Kolyvagin anti-concentration (Note 3)

- **Setup.** Take `389a1`; `K` with the Heegner hypothesis and `E^K` of
  analytic rank 1; `p ∈ {5, 7, 11, 13}`; Kolyvagin primes `ℓ < 2000`;
  admissible primes `q < 200`.
- **Compute.** For each `q`: the Brandt module of level 389 for `B_{q,∞}`,
  the mod-`p` eigenvector `φ_q`, and the Gross point `x_0` for `O_K`. For each
  `ℓ`: the `ℓ + 1` neighbours with their `T_ℓ`-cyclic order, then
  `W_{ℓ,q} ∈ F_p`.
- **Predictions.**
  - some `W ≢ 0` for each `p`;
  - a vanishing fraction of about `1/p`;
  - no excess correlation across different `ℓ`.
- **Kill criteria.**
  - All `W ≡ 0` at a `p` with `Sha[p] = 0` known. This refutes the explicit
    formula or Expectation 3.1.
  - A strongly structured vanishing pattern. This refutes the independence
    model.
- **Cross-check.** Compare with the Jetchev–Lauter–Stein Kolyvagin classes
  for `389a1`.

## E5 — Truncated volumes (Note 5)

- **Compute.**
  `V(x) = (∏_{p≤x} #Ẽ(F_p)/p)^{−1} · (log x)^r · √2 e^{rγ}/r!` against
  `L^(r)(E,1)/r!`, for `x ≤ 10^5`, on `389a1`, `571b1` and `5077a1`.
- **Deliverable.** A convergence profile plus a rigidity map. This is a
  diagnostic, not a kill test.

## RB1 — First-derivative no-go (Note 6) — **RUN, 2026-09-23**

- **Driver:** `tests_research/rb1_heegner_no_go.py`.
- **Data:** `data/rb1_heegner.json`.
- **Method:** Heegner points `y_K` for every fundamental `D`,
  `|D| ≤ 160`, satisfying the Heegner hypothesis, via `φ(τ) = Σ a_n/n qⁿ`
  mod the lab's Néron lattice, at 30 digits.

**Results.**

- `389a1`: 21/21 torsion.
- `5077a1`: 22/22 torsion.
- Control `37a1`: 18 non-torsion points, each exactly `k·P`, and 4 torsion
  points.
- On all 22 control discriminants, the exact modular-symbol value satisfies
  `L_alg(E^D,1) = 2k²`. This is Gross–Zagier with zero error.

**Verdict.** The lemma holds as predicted. Point constructions are dead
for Route B. See Note 6 §2.

## E4 status (2026-09-23)

`𝔅_p = 1` at every resolved good ordinary prime, with no failures:

- `37a1`: all 12 good ordinary `p ≤ 59`. The formerly unresolved
  `p = 13, 47, 53, 59` are now resolved at higher level.
- `389a1`: all 15 good ordinary `p ≤ 59`.
- `571b1`: all 15 good ordinary `p ≤ 59` (first run for this curve).
- `5077a1`: `p ≤ 17`.

That is 42/42 cells for `p ≤ 59` on the three curves of rank ≤ 2.

See `data/e4_beta_table*.json`. This reproduces Stein–Wuthrich 2013. It
is not new evidence for GRS uniformity (Note 1 §8.1).
