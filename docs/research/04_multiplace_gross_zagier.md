# Note 4 — Legs over Q: the product-formula obstruction and geodesic barycenters

**Target obstacles:** new machinery (3) and exact matching (2).
**Status:** §1 is an explanation, §3 is an elementary proposition (proved),
§4 is a conjecture that we expect to be **cheap to kill**.

## 1. Why number-field Gross–Zagier stops at the first derivative

The Gross–Zagier height is a sum of local symbols:

```
⟨y, T_m y⟩ = Σ_{ℓ finite} i_ℓ(y, T_m y) · log ℓ  +  g_∞(y, T_m y).
```

It is well defined on divisor classes because, for a principal divisor
`div(h)`, the local symbols are `log|h(D)|_v` and the **product formula**
`Σ_v log|x|_v = 0` kills their sum. The product formula is *linear* in the
logarithms.

Any other weighting fails. For example
`Σ_ℓ i_ℓ · (log ℓ)²` does not vanish on principal divisors: for
`h(D) = 2/3` it gives `(log 2)² − (log 3)² ≠ 0`. So no quadratic or higher
weighting of places descends to the Mordell–Weil group.

**What Yun–Zhang do over function fields.** Their `r`-th derivative is an
intersection number on a moduli space of shtukas with `r` legs. The
intersection points lie over `r`-**tuples** of points of the curve, weighted
by `∏ deg(x_i)`. The `r` legs carry `r` independent Frobenii, and each
coordinate has its own "degree of a principal divisor is 0" law. **`r` legs =
`r` independent product formulas.** Over `Q` there is exactly one product
formula, so there is only one leg.

**What the analytic side does anyway.** The Gross–Zagier kernel is built from
an Eisenstein series whose Fourier coefficients contain divisor sums
`Σ_{d|n} d^{−s}`. Differentiating `r` times produces `(log d)^r`, and

```
(log d)^r = Σ_{ℓ_1,…,ℓ_r} ∏_i ord_{ℓ_i}(d) · log ℓ_i ,
```

a sum over `r`-tuples of primes. **The kernel has `r` legs; the arithmetic
side cannot.** For a curve with ranks `(r(E), r(E^K)) = (2,1)` the target is
`L'''(E/K,1)/3! = (L''(E,1)/2) · L'(E^K,1)`, which by BSD is proportional to
the rank-3 regulator of `E(K)`. That is a three-legged quantity.

## 2. Where extra independent legs could come from over Q

| candidate leg | its "product formula" | status |
|---|---|---|
| **Kolyvagin primes** `ℓ` (finite tori `T_ℓ ≅ Z/(ℓ+1)`) | A class-sum vanishes mod `p` because `p | a_ℓ` (Note 3 §2) | Known in spirit: this is Darmon's refined Mazur–Tate conjecture for Heegner points (Invent. 1992) and Bertolini–Darmon derived heights. Values are mod `ℓ+1`, never archimedean. |
| **p-adic legs** (several legs at one prime) | The Scholze–Weinstein moduli of local shtukas with several legs, and Drinfeld's lemma for diamonds | Local only. A mixed-characteristic *higher* arithmetic fundamental lemma for CM cycles with `r ≥ 2` legs at `p` is, as far as we found, not formulated. We list it as a direction and do not develop it here. |
| **archimedean legs from real tori** `T(R)` of real quadratic fields (closed geodesics) | The cycle integral vanishes when `L(E,1) = 0` | **New proposal, §§3–4.** |

The third row is the only candidate that lives at the archimedean place,
which is where `L''(E,1)` lives (wall W2 in the [README](README.md)).

## 3. The archimedean Kolyvagin derivative is a geodesic barycenter (proved)

**Setup.**
- `f ∈ S_2(Γ_0(N))` is the newform of `E`, and
  `F(z) = 2πi ∫_{i∞}^{z} f(w) dw`. The map `F mod Λ_E` is the modular
  parametrization `X_0(N) → E(C)`, up to sign and the Manin constant.
- `Q = [a,b,c]` is a form with `N | a` and discriminant `D > 0`, `D` not a
  square. `S_Q` is the geodesic joining its roots, `γ_Q` generates its
  stabiliser in `Γ_0(N)`, and `C_Q = ⟨γ_Q⟩\S_Q` is the closed geodesic, of
  length `ℓ_Q = 2 log ε_D`.
- The cycle integral is `I(Q) = ∫_z^{γ_Q z} f(w) dw`, which does not depend
  on `z`.

**When `I(Q) = 0`.** For `E` with `L(E,1) = 0` and `D` of narrow class number
one with `(D/p) = 1` for every `p | N`, the real-quadratic toric period
formula (Popa, Compositio 2006) gives `|I(Q)|² ∝ L(E,1) L(E^D,1) = 0`. In
that case `F` is `γ_Q`-periodic on `S_Q`, and we define

```
B(Q) := (1/ℓ_Q) ∫_{C_Q} F(z) ds(z)  ∈ C,        P(Q) := B(Q) mod Λ_E  ∈ E(C),
```

where `ds` is hyperbolic arc length.

**Proposition 4.1.**
1. `P(Q)` depends only on the `Γ_0(N)`-class of `Q`.
2. **(Sawtooth identity.)** Take `z_0 ∈ S_Q`, let `θ ∈ [0,1)` be the
   normalised arc-length coordinate from `z_0`, and let `B̄_1(θ) = θ − 1/2`.
   Then
   ```
   ∫_{z_0}^{γ_Q z_0} B̄_1(θ(z)) · 2πi f(z) dz  =  F(z_0) − B(Q).
   ```
3. **(Split-torus degeneration.)** For the non-closed geodesic from `0` to
   `i∞`, the arc-length coordinate is `s = log y`, and for every `k`
   ```
   ∫_{0}^{i∞} (log(z/i))^k f(z) dz  =  i · Λ^(k)(E,1),     Λ(E,s) = (2π)^{−s}Γ(s)L(E,s).
   ```

*Proof.*
1. Replacing `Q` by an equivalent form moves `S_Q` by some `γ ∈ Γ_0(N)`.
   Since `F(γz) = F(z) + ω(γ)` with `ω(γ) ∈ Λ_E`, and `ds` is invariant,
   `B` changes by an element of `Λ_E`.
2. Integrate by parts with `d(F(z) − F(z_0)) = 2πi f(z) dz`. The boundary
   term is `[(θ − ½)(F(z) − F(z_0))]_0^1 = ½ (F(γ_Q z_0) − F(z_0)) = 0`,
   because `I(Q) = 0`. What remains is
   `−∫_0^1 (F(z(θ)) − F(z_0)) dθ = F(z_0) − B(Q)`, using `ds/ℓ_Q = dθ`.
3. Substitute `z = iy`, `dz = i dy`, `log(z/i) = log y`, and use the Mellin
   formula of [Note 2](02_motivic_transfer.md) §2. ∎

**Reading.** The discrete Kolyvagin derivative
`W = Σ_t log_g(t) φ(x_t)` weights a torus orbit by the discrete sawtooth.
Its Fourier coefficients are `1/(ζ − 1)`, the discrete analogue of
`1/(2πik)`. Part 2 is the exact archimedean counterpart: the sawtooth-weighted
orbit integral is the *displacement of the base point from the barycenter
`P(Q)`*. Its shift-invariance is the analogue of `Σ_t φ(x_t) ≡ 0`. Part 3 says
that the **L-derivatives themselves are the arc-length moments of the one
degenerate (split) torus**, so closed real-quadratic geodesics are the natural
"legs" on which to look for rank-2 information.

## 4. The conjecture (expected to be killable)

> **Conjecture BG (barycentric geodesic points).** Let `E/Q` have
> `L(E,1) = 0`, and let `D > 0` be as in §3.
> - **BG-strong:** `P(Q_D) ∈ E(Q) ⊗ Q`. Equivalently, `B(Q_D)` is a
>   `Q`-linear combination of the elliptic logarithms of Mordell–Weil
>   generators plus a period.
> - **BG-span:** for rank-2 `E`, the points `P(Q_D)` span `E(Q) ⊗ Q` as `D`
>   varies.
> - **BG-modular (weaker, independent):** the traces `Σ_{[Q]} B(Q)` over
>   discriminant `D`, extended beyond narrow class number one, are the
>   Fourier coefficients of the holomorphic part of a weight-1/2 mock modular
>   form whose shadow is the Shintani/Kohnen lift of `f`. This is the
>   weight-0 Eichler-integral analogue of the Duke–Imamoḡlu–Tóth theory of
>   cycle integrals of `j`.

**The motivation is honest but thin.**

- In favour: the exact identity of §3 matches the discrete Kolyvagin class,
  and Kolyvagin classes *are* Mordell–Weil elements mod `p`. Darmon's `p`-adic
  Stark–Heegner points are built from integrals along these same geodesics.
- Against: no mechanism is known that forces a non-holomorphic average
  (`ds`) of a holomorphic lift to be algebraic. The most likely outcome is
  **BG-strong false** and **BG-modular true-or-false**.

**Why a kill is still worth having.** A kill of BG-strong would sharpen W2:
it would show that archimedean analogues of Kolyvagin derivatives need
genuinely non-holomorphic height data (Green functions), not just arc-length
averages. That would point the higher-derivative problem toward
Kudla-program derivatives of Eisenstein series rather than toric periods.

**If BG-span survived,** it would be the first construction of a rank-2
Mordell–Weil group from `L`-function-side geometry over `Q`. That is a
direct attack on BARRIER.md §5 item 1.

## 5. Novelty and nearest literature

- **Nearest work.** Yun–Zhang (function fields); Darmon, *Integration on
  `H_p × H`* (Annals 2001) and Darmon–Vonk (`p`-adic, same geodesics);
  Duke–Imamoḡlu–Tóth on cycle integrals of `j` (Annals 2011) and on
  geometric invariants and modular cocycles for real quadratic fields;
  Bruinier–Funke–Imamoḡlu, regularised theta lifts and periods of modular
  functions; Bruinier–Ono (Annals 2010), where weight-1/2 harmonic Maass forms
  detect `L'(E^D,1)`; Popa (real-quadratic toric periods).
- **What we did not find:** the "r legs = r independent product formulas"
  diagnosis applied to number fields; arc-length barycenters of the modular
  parametrization along closed geodesics as candidate points of `E`; and the
  sawtooth identity (Proposition 4.1(2)) linking them to Kolyvagin
  derivatives. These are proposed, not established.

## 6. Kill test

This is experiment **E3** in [`EXPERIMENTS.md`](EXPERIMENTS.md), on `389a1`
(rank 2, `N` prime).

1. **Choose discriminants.** Take ten fundamental `D > 0` with narrow class
   number one and `(D/389) = 1`. Check numerically that `I(Q_D) = 0` at the
   working precision. This is itself a check of the toric period formula.
2. **Compute barycenters.** Compute `B(Q_D)` to 40 digits. Evaluate `F(z)`
   by its `q`-series after moving each point up with `Γ_0(389)` and `w_389`.
   The minimum height is about `√3/(2·389)`, so a few thousand `a_n` suffice.
   Integrate along `C_Q` with Gauss–Legendre quadrature.
3. **Search for relations.** Run PSLQ on
   `{B(Q_D), z(P_1), z(P_2), ω_1, ω_2}`, where `z(P_i)` are the elliptic
   logarithms of the generators and `ω_i` the periods, with coefficient bound
   `10^6`.
   - **BG-strong kill:** no relation for any `D` at 40 digits.
   - **BG-span kill:** all relations found lie on one line in
     `E(Q) ⊗ Q`.
4. **BG-modular.** Tabulate traces for all `D < 200`, fit against a basis of
   the weight-1/2 mock modular space. **Kill:** no fit.
