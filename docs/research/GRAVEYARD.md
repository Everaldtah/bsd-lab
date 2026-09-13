# Graveyard — ideas developed and then killed

Each entry gives the idea, why it looked promising, and the argument that
kills it (or shows it is already known). These are kept on purpose: they
mark dead ends for whoever tries next.

## G1. Congruence-rescaled `p`-adic heights ("isotropy conjecture")

**Idea.**
- Level-raise the rank-2 form `f` at primes `ℓ_k`, giving forms `g_k ≡ f`
  mod `p^k` with flipped sign, so each `g_k` has odd rank.
- Their Heegner points `y_{g_k}` give classes `c_k ∈ Sel_{p^k}(E/Q)`.
- The `p`-adic L-functions satisfy
  `L_p(g_k) ≈ (Euler factor at ℓ_k) · L_p(f)` mod `p^k`, and `L_p(f)` vanishes
  to order 2. So `L_p'(g_k) ≡ 0` mod `p^k`.
- By `p`-adic Gross–Zagier, `h_p(c_k) ≡ 0`: the classes would be isotropic
  for the `p`-adic height form. Rescaling by `p^{−k}` might then recover
  `L_p''(f)`.

**Kill.** Level raising changes the local condition at `ℓ`, and
`Sel_p(g) = ker(loc_ℓ) ⊂ Sel_p(E)` (W. Zhang's rank-lowering). The class of
`y_g` mod `p` therefore spans the line
`L_ℓ = ker(E(Q)/p → E(F_ℓ)/p)`. That line is fixed by the global points and
`ℓ`, not by heights. As `ℓ` varies (Chebotarev), `L_ℓ` runs over **every**
line in `E(Q)/p ≅ F_p²`. Isotropy of every line forces the height form to
vanish mod `p` identically, which is false for generic `p`.

The flaw in the idea: the height congruence ignores the change of the local
height term at `ℓ`. Any corrected version must carry `ℓ`-local terms, and then
it no longer isolates `L_p''(f)`.

## G2. Hessian / pencil regulator in the Hida × cyclotomic plane

**Idea.** At a rank-`r` point of the two-variable `p`-adic L-function, the
degree-`r` Taylor form in `(dk, ds)` should equal
`det(ds·H_cyc + dk·H_wt)`, a pencil of heights.

**Kill (not new).** The algebraic side is Nekovář's generalised Rubin formula
for Selmer complexes, and the analytic side follows from the two-variable
Skinner–Urban main conjecture. Venerucci and Bertolini–Darmon treat the
exceptional and weight-direction cases. It is also purely per-prime, so it
does nothing for uniformity (Note 1).

## G3. Mixed archimedean–Mazur–Tate jets

**Idea.** Take the group-ring-valued complex function
`L_M(s) = Σ_χ L(E,χ,s) e_χ ∈ C[G_M]` and conjecture a mixed leading term in
`(s−1)^a · I^b/I^{b+1}`, interpolating the Mazur–Tate refined conjecture
(`a = 0`) and classical BSD (`a = r`).

**Kill.** Over `R` or `C` the augmentation filtration of a finite group ring
collapses: `I/I² ≅ G_M` is finite, so `I ⊗ R = I² ⊗ R`. Mazur–Tate leading
terms need **integral** coefficients. Derivatives `L'(E,χ,1)` have no common
integral structure, since each twist carries its own transcendental period.
A mixed statement with `0 < a < r` has no meaning.

## G4. The wedge of points as a zero-cycle on `E^r`

**Idea.** Realise `P_1 ∧ ⋯ ∧ P_r` as the cycle
`Σ ± (P_{σ1}, …, P_{σr})` in `CH_0(E^r)`, and look for an L-function formula
for its "height".

**Kill.** The cycle lies in the deepest step `F^r` of the Bloch–Beilinson
filtration, which is conjecturally zero (torsion) over number fields. Its
Albanese-kernel part is expected to be torsion, so there is nothing to
measure. This is also why the Gross–Schoen/Kudla-type cycles that do carry
information live on *modular* varieties, not on powers of `E`.

## G5. Heegner points from several CM fields landing in `E(Q)`

**Idea.** Combine Heegner points from `K_1, K_2`, or work over the
biquadratic field `K_1K_2`, to reach the `E(Q)`-part of a rank-2 curve.

**Kill.** Complex conjugation acts on `y_K` by `−ε(E)` modulo torsion. When
`ε(E) = +1` (even rank), every Heegner point lies in the `−1` eigenspace,
`E^K(Q)`. Over `K_1K_2`, the sign of the real-quadratic twist `E^{D_1D_2}`
is `χ_{D_1D_2}(−N) ε(E) = +1`, and the twists `E^{D_i}` are odd. No
combination reaches `E(Q)`. Heegner points never see even rank. (This is the
precise form of PLAN.md §5.2(a).)

## G6. `(log p)^r`-weighted arithmetic intersection pairings

**Idea.** Replace `log p` by `(log p)^r` in the local Gross–Zagier symbols to
match `∂_s^r` of the kernel.

**Kill.** This fails the product formula, so it does not descend to divisor
classes (example: `h = 2/3`). Kept as the motivation for Note 4 §1: extra legs
need *independent* product formulas.

## G7. Derived deformation rings / eigenvariety jets

**Idea.** Read `rank E(Q)` off the tangent or jet structure of a derived
Galois deformation ring (Galatius–Venkatesh), or off an eigenvariety at a
CAP/Eisenstein point (Bellaïche–Chenevier), with intersection multiplicity
equal to `ord L`.

**Kill.** Derived deformation rings of `ρ_E` have homotopy given by
*adjoint* Selmer groups (`Sym²`), not `H^1_f(Q, V_pE)`. Venkatesh's derived
Hecke action likewise sees the adjoint motive. The Eisenstein/CAP jet version
reproduces the Skinner–Urban divisibility: a *lower* bound on Selmer corank in
terms of the `p`-adic vanishing order, one prime at a time. It gives no upper
bound, no uniformity, and nothing at the archimedean place.

## G8. An archimedean Kato Euler system at the central point

**Idea.** Use Beilinson–Kato elements in `K_2` of modular curves, which are
global motivic objects, together with their archimedean regulator, to get
p-independent control of `L^(r)(E,1)`.

**Kill.** The Beilinson regulator of the Kato elements computes the
*non-critical* value `L(E,2)` (equivalently `L'(E,0)`). The critical value
`L(E,1)` is reached only through `p`-adic explicit reciprocity (the dual
exponential). There is no archimedean regulator interpretation at `s = 1`,
so Kato's system is global at its base but `p`-adic at the point we need.
