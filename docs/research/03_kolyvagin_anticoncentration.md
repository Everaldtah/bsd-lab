# Note 3 — Kolyvagin anti-concentration: a second route to uniformity

**Target obstacles:** control in arbitrary rank (1) and new machinery (3).
**Status:** §1 assembles known theorems (the one "Expectation" is flagged).
§2 gives an explicit formula whose precise form must be re-derived before
use. §§4–6 are new proposals.

## 1. Rank lowering by a parity-chosen twist (any rank `r ≥ 2`)

Let `E/Q` have rank `r ≥ 2` and root number `(−1)^r`. Choose an imaginary
quadratic field `K` such that:

- every prime dividing `N` splits in `K` (the Heegner hypothesis), so the sign
  of `L(E/K,s)` is `−1`;
- `ord_{s=1} L(E^K, s) = r' ∈ {0,1}`.

Parity forces `r' ≡ r + 1 (mod 2)`. Such `K` exist by the twist non-vanishing
theorems used in the classical rank ≤ 1 proof (Bump–Friedberg–Hoffstein;
Murty–Murty). Then:

- **The twist is fully known.** Gross–Zagier plus Kolyvagin over `Q` give
  `rank E^K(Q) = r'` and `Sha(E^K)` finite, hence `Sha(E^K)[p] = 0` for
  `p ≫ 0`.
- **The Heegner point is torsion.** `ord L(E/K) = r + r' ≥ 2`, so
  `L'(E/K,1) = 0` and `y_K` is torsion. Hence `κ_1 = 0` for `p ≫ 0`.
- **Sign bookkeeping.** Complex conjugation acts on `κ_n` by
  `ε·(−1)^{ν(n)+1}`, where `ν(n)` is the number of prime factors of `n`. The
  classes with `ν(n) = r − 1` therefore lie in the `E(Q)`-eigenspace.

**Kolyvagin's structure theorem** (Math. Ann. 1991) combined with **W. Zhang's
proof of Kolyvagin's conjecture** (Camb. J. Math. 2014, valid for `p` outside
an explicit finite set) says: with `ν_p` the minimal `ν(n)` for which
`κ_n ≠ 0` in the `p^∞` sense, one eigenspace of `Sel_{p^∞}(E/K)` has corank
`ν_p + 1` and the other has corank at most `ν_p`. Our eigenspaces have
coranks `r + t_p` (where `t_p` is the corank of `Sha(E)[p^∞]`) and `r' ≤ 1`.
So

```
ν_p = r − 1 + t_p ,     and     Sha(E)[p^∞] finite  ⟺  ν_p = r − 1.
```

> **Expectation 3.1 (mod-p refinement, to be checked against Kolyvagin 1991
> and Zhang 2014 §§9–10).** For `p` outside an explicit finite set,
> `Sha(E)[p] = 0` iff some `κ_n` with `ν(n) = r − 1` is nonzero in
> `H^1(K, E[p])`.

**Consequence.** In every rank, finiteness of `Sha` for a curve of known rank
comes down to the **uniform-in-p non-vanishing mod p of Kolyvagin classes at
depth exactly `r − 1`**. Kolyvagin's machine controls every rank one prime at
a time. Uniformity is the only wall, as [Note 1](01_global_rational_shadow.md)
found by a different route. Jetchev–Lauter–Stein (2009) computed such classes
for `389a1` at one small prime.

## 2. An explicit form of the depth-1 class (`r = 2`)

Fix `p`, a Kolyvagin prime `ℓ` (inert in `K`, with `p | ℓ+1` and `p | a_ℓ`),
and an **admissible** prime `q` in the sense of Bertolini–Darmon (Annals 2005):
`q ∤ NpD_K`, `q` inert in `K`, `p ∤ q² − 1`, and `p | (q+1)² − a_q²`. For such
`q`, `H^1_fin(K_q, E[p])` is 1-dimensional. By Chebotarev, `κ_ℓ ≠ 0` iff
`loc_q κ_ℓ ≠ 0` for some admissible `q`.

Here is the computation.

- **Reduction to supersingular points.** Since `q` is inert in `K`, Heegner
  points reduce mod `q` to supersingular points. These correspond to the finite
  set `X_q = Pic` of an Eichler order of level `N` in the definite quaternion
  algebra `B` ramified at `{q, ∞}`.
- **The level-raised eigenform.** By Ihara's lemma and Ribet, the reduction
  map `Z[X_q]^0 → Ẽ(F_{q²})/p` factors through a mod-`p` eigenform
  `φ_q : X_q → F_p`. Its Hecke eigenvalues are those of `f` mod `p`, level-raised
  at `q`.
- **The torus action.** The `G_ℓ = Gal(K[ℓ]/K[1])`-conjugates of `y_ℓ` reduce
  to the `ℓ + 1` neighbours `x_t` of `x_0 = red(y_1)` in the `ℓ`-isogeny tree.
  The nonsplit torus `T_ℓ = (O_K ⊗ F_ℓ)^×/F_ℓ^× ≅ Z/(ℓ+1)` permutes these
  neighbours simply transitively.

Tracing through `Pic(O_K)`, the expected form of the localized class is

```
loc_q κ_ℓ  ≐  W_{ℓ,q} := Σ_{σ ∈ Pic O_K}  Σ_{t ∈ T_ℓ}  log_g(t) · φ_q(σ · x_t)   ∈ F_p
```

where `log_g : T_ℓ → Z/(ℓ+1) → Z/p` is the discrete logarithm to a generator
`g`, and `≐` means "up to `F_p^×`". Two consistency checks pass:

- **Independent of the choice of base point for `log`.** The unweighted sum is
  `Σ_t φ_q(σx_t) = (T_ℓ φ_q)(σx_0) = a_ℓ φ_q(σx_0) ≡ 0 (mod p)`, because
  `p | a_ℓ`. Shifting `log` by a constant therefore changes nothing.
- **Independent of the generator up to a unit.** Replacing `g` by `g^u`
  multiplies `W` by `u^{-1}`.

**Fourier form.** For characters `χ` of `T_ℓ` with `ζ = χ(g) ≠ 1`, one has
`Σ_{i=0}^{ℓ} i ζ^i = (ℓ+1)/(ζ − 1)`. Let `P_χ` be the toric period of a
characteristic-0 lift of `φ_q`. Then

```
W = Σ_{χ ≠ 1} P_{χ̄}/(χ(g) − 1)  +  (ℓ/2) · P_1 .
```

Only characters of `p`-power order have non-unit `χ(g) − 1`, so `W mod p`
is the **first Taylor coefficient at `χ = 1`** of `χ ↦ P_χ` along the
`p`-part of the ring class group. By the Gross–Waldspurger formula,
`|P_χ|² ≐ L^alg(g/K, χ, 1)`. Non-vanishing of `κ_ℓ` is thus the non-vanishing
mod `p` of a *derivative* of algebraic twisted central L-values of a
level-raised form. This is an anticyclotomic Kurihara-type number, in the
spirit of C.-H. Kim (arXiv:2203.12161).

## 3. Why analytic equidistribution cannot decide it

The natural tools — Duke/Michel–Venkatesh equidistribution of CM points, and
the Ramanujan bound `|λ| ≤ 2√ℓ` for the `ℓ`-Hecke operator on `X_q` — control
a weighted sum `Σ_t w(t) F(x_t)` for real-valued `F`: they give the main
term plus an error `O(√ℓ)`. `W mod p` depends on the counts
`N_{c,v} = #{t : log t ≡ c (mod p), φ_q(x_t) = v}` **modulo `p`**, and an
error of `O(√ℓ)` cannot determine a residue class mod `p`. The residues of
these counts are *arithmetic* invariants, invisible to archimedean mixing.

The two known sources of *exact* mod-`p` information are:

1. **Chebotarev in a fixed field.** Useless here: `κ_ℓ` lives over `K[ℓ]`,
   which varies with `ℓ`, so `ℓ ↦ κ_ℓ` is not a Frobenius function in any
   fixed extension.
2. **Cornut–Vatsal's simultaneous surjectivity** (via `p`-adic Ratner). This
   is exact, but it applies to *vertical* towers `ℓ_0^k`, not to the
   *horizontal* family of Kolyvagin primes with `p | ℓ + 1`.

This gap is precisely why uniformity is open.

## 4. How much independence is needed (Borel–Cantelli count)

If `Sha(E)[p] ≠ 0`, Expectation 3.1 forces `W_{ℓ,q} ≡ 0` for **every**
Kolyvagin `ℓ` and admissible `q`. Delaunay's heuristics give
`Prob(p | #Sha) ≈ p^{−(2r+1)}` for rank `r`, which is `p^{−5}` at rank 2.
That is summable, which is consistent with `Sha` being finite.

A proof needs a mechanism giving at least **two effectively independent**
`F_p`-tests per prime, so that the failure probability is `≲ p^{−2}` and
still summable. One test per prime is not enough: `Σ 1/p` diverges, as in
Note 1 §4. The Kolyvagin data offer infinitely many pairs `(ℓ, q)`. The task is
to prove enough independence.

## 5. Proposal: inverse Littlewood–Offord ⇒ hidden symmetry ⇒ contradiction

We propose an argument by contradiction, in four steps.

1. **Ultraproduct setting.** Suppose `Sha(E)[p] ≠ 0` for infinitely many `p`.
   Pass to a non-principal ultrafilter concentrated on them (Note 1 §5). The
   graphs `X_{q_p}`, the forms `φ_{q_p}` and the primes `ℓ_p` become
   non-standard objects over `F_U`, and every `W ≡ 0` holds exactly.
2. **Inverse theorem.** The weights `log_g(t)` take every residue mod `p`
   equally often. They are maximally unstructured. Inverse Littlewood–Offord
   theorems over `F_p` (Tao–Vu; Nguyen–Vu) say that if a linear form with such
   weights vanishes on *many* independent configurations, then the
   coefficients must carry strong additive structure. Here the coefficients
   are the values of `φ_q` along `T_ℓ`-orbits.
3. **Structure ⇒ symmetry.** Vanishing for all Kolyvagin `ℓ` should force
   `φ_q`, on `T_ℓ`-orbits of the CM point, to factor through a coarse
   invariant of the orbit (a sub-torus quotient). A mod-`p` eigenform
   invariant under such sub-torus quotients for a positive density of `ℓ`
   commutes with extra correspondences. That is expected to force the image of
   `ρ̄_{E,p}` to be dihedral or CM.
4. **Contradiction.** Serre's open image theorem makes `ρ̄_{E,p}` surjective
   for `p ≫ 0` when `E` has no CM. (CM curves need the separate theory of
   Castella, arXiv:2204.09608.)

**Weakest step:** step 3. It is unproven and may be false as stated. The
honest research question is:

> *Classify mod-`p` quaternionic eigenforms whose discrete-log-weighted
> `T_ℓ`-orbit sums vanish for a positive density of Kolyvagin primes `ℓ`.*

If the only such forms are those attached to non-surjective `ρ̄`, uniformity
follows for every rank-2 non-CM curve.

## 6. Rank `r`: a multilinear version

At depth `ν(n) = r − 1` with `n = ℓ_1 ⋯ ℓ_{r−1}`, the class localizes to

```
W_{n,q} = Σ_σ Σ_{t_1,…,t_{r−1}}  ∏_i log_{g_i}(t_i) · φ_q(σ · x_{t_1…t_{r−1}}),
```

a multilinear form over the product of the torus orbits. The vertices
`x_{t_1…t_{r−1}}` are neighbours in a product of `ℓ_i`-trees. The
combinatorial engine becomes **polynomial anti-concentration** (Meka–Nguyen–Vu
2016) on products of Ramanujan graphs, and the Borel–Cantelli requirement is
unchanged. This is the note's arbitrary-rank proposal. It needs no rank-`r`
Euler system: only depth-`(r−1)` Kolyvagin classes, which exist in every rank.

## 7. Novelty and nearest literature

- **Nearest work.** Kolyvagin 1991 (structure); W. Zhang 2014 (conjecture);
  Bertolini–Darmon 2005 (admissible primes, reciprocity); Cornut–Vatsal 2005
  (Ratner surjectivity); Jetchev–Lauter–Stein 2009 (explicit classes on
  `389a1`); C.-H. Kim arXiv:2203.12161 (Kolyvagin systems vs Kurihara numbers,
  arbitrary rank); Sweeting arXiv:2012.11771 and Angurel arXiv:2605.26917
  (ultra-Kolyvagin systems over `p^k`, fixed `p`); Delaunay 2001.
- **What we did not find stated anywhere:** the uniform-in-`p` problem posed
  through the explicit weighted orbit sum `W`; the argument that archimedean
  equidistribution cannot decide it; the Borel–Cantelli independence count;
  and the inverse-Littlewood–Offord ⇒ symmetry ⇒ open-image architecture,
  including its multilinear rank-`r` form. Existing ultraproducts run over
  powers of a *fixed* `p`; ours runs over *all primes*.

## 8. Kill test

This is experiment **E2** in [`EXPERIMENTS.md`](EXPERIMENTS.md), on `389a1`.

- **Setup.** Take `K` with the Heegner hypothesis and `E^K` of analytic rank 1
  (the smallest such `|D_K|`, found by the lab). Use `p ∈ {5, 7, 11, 13}`,
  Kolyvagin `ℓ < 2000` and admissible `q < 200`.
- **Compute.** Build Brandt matrices (Pizer's algorithm), extract `φ_q`,
  locate the CM point, and evaluate `W_{ℓ,q}`.
- **Predictions.**
  - For each `p`, some `W ≢ 0`. This is consistent with `Sha(389a1)[p] = 0`,
    which our table and 2-descent support.
  - The fraction of vanishing pairs is `≈ 1/p`.
  - There is no correlation between `W_{ℓ,q}` and `W_{ℓ',q}` beyond chance
    (the independence needed in §4).
- **Kill criteria.**
  - All `W ≡ 0` at some `p` where `Sha[p] = 0` is independently known. This
    would refute the explicit form in §2 or Expectation 3.1 (not BSD).
  - A vanishing fraction far above `1/p` with visible structure. This would
    refute the independence model.
- **Cross-check.** Compare with the Jetchev–Lauter–Stein classes at `p = 3`.
