# Note 2 — Realization-independent generating series (how one would prove GRS)

**Target obstacles:** exact matching (2) and new machinery (3).
**Status:** a program. §2 records elementary facts, §3 is a conjecture, and
§4 explains why rank 1 already follows this pattern.

## 1. The problem this note solves

[Note 1](01_global_rational_shadow.md) reduced uniformity in `p` to **GRS**:
the ratios `𝔅_p = L_p^*/(E_p·Reg_p)` are one rational number, independent of
`p`. The archimedean analogue `𝔅_∞ = L^(r)(E,1)/(r!·Ω·Reg)` should be that same
rational. A proof of GRS must produce an identity that holds **in every
realization at once**: archimedean, and `p`-adic for every `p`. It must do so
without first knowing `#Sha`.

This note makes one claim. **Rank 1 already has such an identity, and it is a
modular generating series. The generalisation to rank `r` needs an
`r`-linear generating series** — the object of
[Note 4](04_multiplace_gross_zagier.md).

## 2. The shape of the leading terms (elementary)

Write `f(τ) = Σ a_n q^n` and `Λ(E,s) = (2π)^{-s} Γ(s) L(E,s)`. Then

```
Λ(E,s) = ∫_0^∞ f(iy) y^{s-1} dy,       so
Λ^(k)(E,1) = ∫_0^∞ f(iy) (log y)^k dy,
L^(r)(E,1) = 2π · Λ^(r)(E,1)            when ord_{s=1} L = r.
```

The last line holds because every lower derivative of `L` vanishes, so Leibniz
leaves one term. On the `p`-adic side, with `μ_f` the Mazur–Swinnerton-Dyer
measure,

```
L_p^(r) = ∫_{Z_p^×} (log_p x)^r dμ_f(x).
```

Both leading terms are an **r-th power of a single logarithm** integrated
against an object attached to `f`. At `∞` the logarithm is `log y`; at `p` it
is `log_p`. The regulators have a different shape: `Reg = det⟨P_i,P_j⟩`, where
each height is a *sum over all places* of local terms. Away from the
archimedean place and `p`, these local terms are rational multiples of
`log ℓ`, respectively of `log_p ℓ`.

So in each realization, `log ℓ` and `log_p ℓ` are the two realizations of one
Kummer period `[ℓ] ∈ Q^× ⊗ Q`. The heights are therefore
realization-compatible **provided the non-Kummer local terms are too**: Green
functions at `∞`, Coleman–Gross local heights at `p`. The comparison between
these two is the biextension picture of Hain (Duke 1990), Scholl ("Height
pairings and special values of L-functions", 1994), Nekovář (1993), and Besser
(the Coleman–Gross and Nekovář `p`-adic heights agree).

The L-side has no such decomposition, because `log y` is not a Kummer period
of any algebraic number. **The L-derivative is the obstruction to a direct
motivic comparison.** The only known way around it is a *kernel identity*
that rewrites the L-derivative as a sum of local symbols.

## 3. The proposal: a formal calculus of local symbols

Let `𝓢` be the free `Q`-vector space on these formal symbols:

- `[ℓ]` for primes `ℓ` (Kummer symbols);
- `G(x, y)` for pairs of CM or Heegner divisors on `X_0(N)`, at level `m`
  (formal Green/Coleman symbols, one per unordered pair and Hecke index);
- `E_v` for the Euler-factor symbols at the auxiliary place `v`.

Each place defines a **realization** `ρ_v : 𝓢 → C_v`:

| symbol | `ρ_∞` | `ρ_p` |
|---|---|---|
| `[ℓ]` | `log ℓ` | `log_p ℓ` (Iwasawa branch, `log_p p = 0`) |
| `G(x,y)` | Gross–Zagier archimedean Green function `g_s` at `s = 1` | Coleman–Gross local height at `p` (ordinary splitting) |
| `E_v` | `1` | `(1 − 1/α_p)^{±1}`-type factors |
| holomorphic projection `Hol` | Sturm/Gross–Zagier `Hol` | Hida ordinary projector `e` |

> **Conjecture RIGS (Realization-Independent Generating Series).** Suppose a
> `q`-expansion identity `Σ_m A_m q^m = Hol(𝒦)` holds **formally** in
> `𝓢[[q]]`, where the `A_m ∈ 𝓢` are arithmetic intersection symbols and `𝒦`
> is an analytic kernel whose Fourier coefficients are also written in `𝓢`.
> Then applying `ρ_v` gives a valid identity for every place `v`, modulo an
> explicitly bounded Euler-factor correction.

**Why this proves GRS.** If the `f`-isotypic projection of such an identity
equates "`L`-leading term" with "`𝔅 ×` (formal regulator)", and the rational
number `𝔅` is produced *formally* (it lives in `Q`, not in any `C_v`), then
`𝔅_∞ = 𝔅_p = 𝔅` for all `p`. This is GRS, obtained without computing `Sha`.
As a by-product, `ord L_p = ord L` holds wherever the formal regulator
realizes to a nonzero number.

## 4. Rank 1 already works this way

In rank 1:

- **Archimedean**: Gross–Zagier compute `⟨y_K, T_m y_K⟩ = Σ_ℓ (local) [ℓ] +
  G(·,·)`. They show the generating series equals the holomorphic projection
  of `θ_K · ∂_s E(·, s)` at `s = 1`. The `f`-component gives
  `L'(E/K,1) = c · ĥ(y_K)`.
- **p-adic**: Perrin-Riou (Invent. 1987) proves the same identity with
  `p`-adic local heights and the ordinary projector, giving
  `L_p'(E/K) = c_p · E_p · h_p(y_K)`.
- **The constants match.** `c` and `c_p` differ exactly by the realization of
  the same formal constants. Therefore
  `L'(E/K,1)/ĥ(y_K)` and `L_p'(E/K)/(E_p·h_p(y_K))` are the *same rational
  multiple* of the respective periods. That is GRS in rank 1, obtained
  exactly as RIGS predicts.

Perrin-Riou's proof is not literally a formal substitution. It re-does the
local computations `p`-adically. **RIGS proposes to prove the substitution
principle once**, for a class of kernels, so that every future Gross–Zagier-type
identity is automatically valid in every realization. Nekovář (Heegner cycles),
Disegni (Shimura curves) and Kobayashi (supersingular) re-proved
`p`-adic analogues one by one. That repetition is the evidence that a
uniform principle exists.

## 5. What rank `r` then needs

GRS in rank `r` follows from RIGS once there is a formal identity whose
`f`-isotypic part equates the `r`-th leading term with an `r`-linear
determinant of formal height symbols. The kernel must be `r`-linear in
Kummer symbols: it has to contain terms like `[ℓ_1][ℓ_2]⋯[ℓ_r]`, because
`det⟨P_i,P_j⟩` expands into products of `r` local symbols at possibly
different places. [Note 4](04_multiplace_gross_zagier.md) shows that
`∂_s^r` of the Gross–Zagier kernel produces exactly such products, and
explains why they do not descend to Mordell–Weil through the product formula.
That is the remaining hard point.

## 6. Novelty and nearest literature

- **Nearest work.** Compatibility of BSD with `p`-adic BSD is expected from
  Fontaine–Perrin-Riou and equivariant-TNC formalisms (fundamental lines,
  Kato's zeta elements). Scholl's blended extensions give a motivic reading of
  the first derivative. Perrin-Riou, Nekovář, Disegni and Kobayashi prove
  `p`-adic Gross–Zagier formulas case by case.
- **What we did not find stated anywhere:**
  1. GRS as the *specific* target that uniformity requires (Note 1);
  2. a formal symbol calculus in which Gross–Zagier-type kernel identities
     are proved once and then realized at every place (RIGS);
  3. the observation that rank `r` needs a kernel `r`-linear in Kummer
     symbols across places, tying the uniformity problem to the higher
     Gross–Zagier problem.
- **Main risks.**
  - Holomorphic projection and ordinary projection may not be formally
    parallel beyond rank 1. The known `p`-adic proofs contain genuinely
    `p`-adic estimates (for example the boundedness of `p`-adic Green
    functions) that have no formal counterpart.
  - Supersingular `p` need signed projectors.

## 7. Kill tests

This is experiment **E4** in [`EXPERIMENTS.md`](EXPERIMENTS.md).

- **E4a (calibration, rank 1).** For `37a1` with generator `P`, compute
  `L_p'(E)/(E_p·h_p(P))` at every good ordinary `p < 100`. The prediction is
  the same rational number `#Sha·∏c/|tors|² = 1` at every `p`. This is not
  new mathematics; it validates the lab's `p`-adic pipeline.
- **E4b (rank 2, the real test).** For `389a1` and `571b1`:
  - check `ord L_p = 2` (so `L_p' = 0` and `L_p'' ≠ 0`) at every good
    ordinary `p < 100`;
  - check that `L_p''/(2·E_p·Reg_p)` equals `𝔅_∞ = L''(E,1)/(2·Ω·Reg) = 1`.
  - **Kill:** any `p` with `L_p'' = 0`, or `𝔅_p ≠ 𝔅_∞`, beyond certified
    precision.
- **E4c (formal-kernel test for RIGS).** Compute the first 50 coefficients of
  Gross–Zagier's generating series for `37a1` symbolically in `𝓢` (Kummer
  symbols plus Green symbols). Realize them at `∞` and at `p = 5, 7` and check
  that one formal identity gives both numerical identities. **Kill:** a
  coefficient where the `p`-adic realization needs a correction not
  expressible in `𝓢`.
