# Note 5 — The exact formula as a volume identity on Bloch's group

**Target obstacle:** exact matching (2).
**Status:** the most speculative note. §3 contains a consistency check that
we believe is new and informative. §4 explains where the approach may be
circular.

## 1. Bloch's reformulation

Bloch (Invent. 58, 1980) considers the semi-abelian variety

```
0 → T → X → E → 0,        T = Hom(E^∨(Q)/tors, G_m),
```

the extension of `E` by a split torus of dimension `r` classified by the
tautological element of `Ext¹(E, T) = Hom(X^*(T), E^∨(Q))`. He shows that
(given finiteness of `Sha`) the full BSD formula for `E` is equivalent to a
**Tamagawa number formula** for `X`:

```
τ(X) = #Pic(X)_tors / #Sha¹(X).
```

This is exactly the shape of **Ono's formula** for algebraic tori (Ann. of
Math. 1963), `τ(T) = #Pic(T)/#Sha¹(T)`. Ono's formula is proved, and the
analytic input it needs is the class number formula: Artin L-functions at
`s = 1`.

## 2. The proposal: Manin–Peyre harmonic analysis on `X`

For tori and for vector groups, Tamagawa-type constants have been computed by
**harmonic analysis of height zeta functions**: Batyrev–Tschinkel for toric
varieties, Chambert-Loir–Tschinkel for equivariant compactifications of vector
groups. The method is Poisson summation on `X(A)/X(Q)`:

- the **trivial character** contributes an Euler product of local height
  integrals, whose poles and leading constant carry the Tamagawa volume;
- **non-trivial characters** carry cohomological corrections. In Peyre's
  constant these are the Brauer-group and `Sha`-type terms.

**Proposal.** Run the same programme on an equivariant compactification of
Bloch's `X`.

- **Local integrals at good `p`.** Integrate the height over `X(Q_p)`,
  stratifying by reduction. The `E`-part should produce the factors
  `#Ẽ(F_p)/p = L_p(E,1)^{−1}`, which are exactly the convergence factors that
  bring in `L(E,s)`.
- **The dual side.** The Pontryagin dual of `E(A)/closure(E(Q))` is controlled
  by Poitou–Tate duality through Selmer groups of the dual, so `Sha` enters
  here, as it does in Ono's formula.
- **The archimedean and bad places** give `Ω` and `∏ c_p`.

## 3. A consistency check: why Bloch's torus is needed (√Reg versus Reg)

Count points of `E` alone. The Mordell–Weil lattice `(E(Q)/tors, ĥ)` gives

```
#{P : ĥ(P) ≤ H} ~ vol(B_r) · H^{r/2} / √Reg .
```

That is a **square root** of the regulator. A lattice-point count on `E` can
never produce the `Reg¹` that appears in BSD, whatever the normalisation.

In `X`, the rational points over `P ∈ E(Q)` form a torsor under `T(Q)`, and
the character lattice of `T` is `E^∨(Q) ≅ E(Q)`. The covolume of `X(Q)` in the
norm-one adelic points is governed by the **pairing** between `E(Q)` and its
dual lattice through the biextension height. The Gram determinant of that
pairing is `det⟨P_i, P_j⟩ = Reg` to the first power. Poisson summation between
a lattice and its dual turns `Reg^{−1/2} · Reg^{−1/2}` from the two
theta-function asymptotics into a single `Reg^{±1}`, which is exactly BSD's
exponent.

**Conclusion.** The regulator in BSD is the covolume of a *pairing of two
lattices*, not a point count. Any volume-theoretic proof must go through an
object like Bloch's `X` that carries both lattices. This rules out a family of
naive "count rational points of `E` in adelic boxes" approaches before anyone
builds them.

## 4. What would be proved, and the circularity risk

Split the exact formula into two statements.

- **(H) Harmonic identity.** With a Tamagawa measure regularised by
  **truncated Euler products** `∏_{p≤x}`, Poisson summation gives an Ono-type
  identity. It expresses the regularised volume of `X(Q)\X(A)` through `Ω`,
  `∏ c_p`, `Reg`, `|tors|` and a Selmer/`Sha` term.
- **(A) Euler-product asymptotic.** `∏_{p≤x} #Ẽ(F_p)/p ~ C (log x)^r`. By
  Goldfeld (C. R. Acad. Sci. 1982) this implies the Riemann hypothesis for
  `L(E,s)`, `r = ord_{s=1} L`, and `C = √2 · e^{rγ} · r!/L^(r)(E,1)`, in the
  normalisation of Goldfeld and Conrad (Canad. J. Math. 2005, which explains the
  `√2`).

Then **(H) + (A) ⇒ the BSD formula.** The attraction is that (H) is pure
arithmetic harmonic analysis and (A) is pure analysis. Neither involves Euler
systems.

**Circularity risks, stated plainly.**

1. **Bloch's own theorem.** Bloch's equivalence uses the full `L`-function as
   convergence factors, so `τ(X) = ` Ono ratio *is* BSD. The new content has
   to come from defining the measure by truncated products and proving (H)
   independently. Whether this is possible is unknown.
2. **Finiteness of `Sha`.** Ono's proof for tori uses that `Sha¹(T)` is
   finite, which is known. For `X`, the Poitou–Tate term involves the full
   Selmer group. If `Sha` had a divisible part, the dual side of (H) would be
   infinite. (H) could therefore at best prove "`Sha` finite ⟺ the analytic
   side is finite". That is not nothing, but it is not the formula.
3. **Difficulty of (A).** (A) is at least as hard as the Riemann hypothesis
   for `L(E,s)`, with extra conditions (Kuo–Murty 2005). This route trades BSD
   for GRH-type analysis. That is a different wall, not a lower one.

## 5. Novelty and nearest literature

- **Nearest work.** Bloch 1980; Ono 1963; Peyre (the Manin constant);
  Batyrev–Tschinkel; Chambert-Loir–Tschinkel; Goldfeld 1982; Conrad 2005;
  Kuo–Murty 2005.
- **What we did not find:**
  1. Manin–Peyre height-zeta harmonic analysis on Bloch's semi-abelian group
     as a strategy for the exact formula;
  2. the √Reg-versus-Reg argument that any volume proof needs a pairing of
     two lattices;
  3. the explicit split of the exact formula into a harmonic identity (H) and
     an Euler-product asymptotic (A).

## 6. Diagnostic experiment

This is experiment **E5** in [`EXPERIMENTS.md`](EXPERIMENTS.md). It is not a
kill test: this note's kill criterion is the theoretical circularity in §4.

- **Compute.** For `389a1`, `571b1` and `5077a1`, compute
  `V(x) = (∏_{p≤x} #Ẽ(F_p)/p)^{−1} · (log x)^r · √2 e^{rγ}/r!` for `x` up to
  `10^5`, using point counts or baby-step giant-step `a_p`. Compare with
  `L^(r)(E,1)/r!` from our certified table.
- **Deliverables.**
  - the convergence profile, which is known to oscillate with the low zeros
    of `L(E,s)`;
  - the ratio at `x = 10^5`;
  - a "rigidity map" (PLAN.md L6(iii)) of how far truncated volumes are from
    the exact leading term at each rank.
