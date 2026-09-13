# Note 1 — The Global Rational Shadow principle

**Target obstacles:** control in arbitrary rank (1) and exact matching (2).
**Status:** a reformulation plus a conjecture. Proposition 1 and Lemma 2 are
proved below; everything else is a proposal.

## 1. How uniformity in `p` is actually obtained today

Showing that `Sha(E/Q)` is finite means two things: `Sha[p^∞]` is finite for
every `p`, **and** `Sha[p] = 0` for all but finitely many `p`. The second part
is always the harder one. In ranks 0 and 1 it comes from one **global
rational number** whose prime factors control every `p` at once:

| rank | the global rational shadow | how it controls all `p` |
|---|---|---|
| 0 | `L(E,1)/Ω ∈ Q^×` (Manin–Drinfeld) | Kato's Euler-system divisibility bounds `#Sha[p^∞]` by the `p`-part of this number, one prime at a time. Only finitely many primes divide its numerator. |
| 1 | Heegner index `[E(K) : Z·y_K] ∈ Z` | Kolyvagin: for `p` outside an explicit finite set, `#Sha(E/K)[p^∞]` divides the `p`-part of the index squared. |

The per-prime machine is different in each row. The **uniformity** always
comes from a single number that does not depend on `p`.

## 2. What rank 2 has: per-prime leading-term control

Let `p` be a prime of good ordinary reduction with unit root `α_p`. Let
`L_p(E,T)` be the Mazur–Swinnerton-Dyer `p`-adic L-function and write
`E_p = (1 − 1/α_p)²`. Let `Reg_p` be the `p`-adic regulator (the
Schneider/Mazur–Tate cyclotomic height determinant, with the standard
normalisation by `log_p` of a topological generator). The `p`-adic BSD
conjecture of Mazur–Tate–Teitelbaum predicts

```
ord_{T=0} L_p(E,T) = r,     L_p^*(E) = E_p · Reg_p · #Sha · ∏c_v / |E(Q)_tors|²
```

Put together Schneider (Invent. 1985), Perrin-Riou, and Nekovář (*Selmer
complexes*, Astérisque 310) with the Skinner–Urban main conjecture. For every
`p` where the Skinner–Urban hypotheses hold, `Reg_p ≠ 0`, and
`rank E(Q) ≥ r`, one expects:

```
L_p^(r) ≠ 0   ⇒   corank Sel_{p^∞}(E/Q) = r,   Sha[p^∞] finite,   and
ord_p #Sha[p^∞] = ord_p L_p^* − ord_p E_p − ord_p Reg_p − ord_p(∏c_v/|tors|²).
```

The "unit" case of this ("normalised `L_p''` is a unit iff `Reg_p` is a unit
and `Sha[p^∞] = 0`") is what arXiv:2609.08431 proves for rank-2 CM curves
outside an explicit set of primes. The general valuation form must be checked
carefully against the control theorem and the finite-submodule issues. It is
stated here as the expected per-prime input.

Define the **p-adic BSD quotient**

```
𝔅_p(E) := L_p^*(E) / (E_p · Reg_p)  ∈  Q_p .
```

## 3. The principle

> **Conjecture GRS(E) (Global Rational Shadow).** There is one `𝔅 ∈ Q^×`
> such that `𝔅_p(E) = 𝔅` for every good ordinary `p` with `Reg_p ≠ 0`. The
> same `𝔅` also appears for the signed leading terms at supersingular `p`
> (Pollack `±`, Kobayashi, Sprung `♯/♭`).
>
> **Weak GRS(E).** Outside a finite set of primes, `𝔅_p(E)` takes only
> finitely many values, all of them in `Q^×`.

- **BSD + p-adic BSD ⇒ GRS**, with `𝔅 = #Sha·∏c_v/|tors|²`. GRS says nothing
  about `Sha`. It only asserts that the `p`-adic leading terms of different
  primes are *shadows of one rational number*. It is a comparison statement
  between realizations, which is what makes it approachable by motivic means
  ([Note 2](02_motivic_transfer.md)).
- **Weak GRS + the per-prime input of §2 ⇒ `Sha[p] = 0` for all but finitely
  many `p`.** This holds at the ordinary primes, and at the supersingular
  primes once signed versions exist (infinitely many supersingular primes
  exist by Elkies 1987, so this is required, not optional). The reason: finitely
  many rationals have finitely many prime factors in total.
- **The route to a historic first.** For `389a1`, rank 2 is known (points plus
  2-descent). GRS plus the per-prime input would leave finitely many primes to
  treat by explicit computation (per-prime criteria, descent). Result:
  **`Sha(389a1)` finite**, which would be the first rank-2 curve for which this
  is known.

## 4. A Borel–Cantelli no-go for "one test per prime"

This heuristic explains why something like GRS is *needed*, not just
convenient.

**Model.** Say a per-prime argument decides `Sha[p] = 0` from a `p`-adic
quantity `X_p` (for example "`Reg_p` is a unit"). Suppose `X_p` behaves like
an independent random element of `Z_p`, failing to be a unit with probability
`≍ 1/p`. Since `Σ 1/p = ∞`, the second Borel–Cantelli lemma says the argument
is silent at **infinitely many** primes almost surely. It can never yield "for
all `p ≫ 0`".

**Evidence that the model is real.** arXiv:2609.08431 finds `Reg_p`
non-unit at 3 of 8,050 tested (curve, `p`) pairs. The ordinary primes of a CM
curve are the split primes, so a `1/p` model predicts about
`5 × ½ Σ_{5<p<30000} 1/p ≈ 5` such pairs over five curves. That is the same
order as what was observed. The same model predicts infinitely many non-unit
primes per curve, so **unit criteria alone cannot give uniformity**.

**The two escapes.**

1. **A global shadow**: `X_p` is a fixed rational number, so the failures are
   deterministic and finite (Note 1, Note 2).
2. **At least two independent tests per prime**: failure probability
   `≍ 1/p²`, which is summable, so the first Borel–Cantelli lemma leaves only
   finitely many failures ([Note 3](03_kolyvagin_anticoncentration.md)).

Rank 1 already follows escape 1: the Heegner index is one integer.

## 5. An ultraproduct reformulation (proved)

Let `𝒫` be the set of primes, `U` a non-principal ultrafilter on `𝒫`, and
`F_U = ∏_U F_p`. This is a pseudo-finite field of characteristic 0. Define the
`F_U`-vector space `Sel_U(E) := ∏_U Sel_p(E/Q)`.

**Proposition 1.** Let `E/Q` have `rank E(Q) = r`. The following are
equivalent:
1. `Sha(E/Q)[p] = 0` for all but finitely many `p`;
2. `dim_{F_U} Sel_U(E) = r` for **every** non-principal ultrafilter `U`.

*Proof.* For `p > 7`, `E(Q)[p] = 0` (Mazur), so
`dim_{F_p} Sel_p(E/Q) = r + dim Sha[p] =: r + s_p`. "`V` has a basis of
exactly `d` vectors" is a first-order statement in the two-sorted structure
`(F_p, V)`, so by Łoś, `dim Sel_U = d` iff `s_p = d − r` for `U`-most `p`.
If (1) holds, the set `{p : s_p = 0}` is cofinite, so it lies in every
non-principal `U`, which gives (2). If (1) fails, the set
`S = {p : s_p ≥ 1}` is infinite. Pick a non-principal `U ∋ S`. Then
`dim Sel_U ≥ r + 1` (or `Sel_U` is infinite-dimensional), and (2) fails. ∎

**Corollary.** `Sha(E/Q)` is finite iff (2) holds and `Sha[p^∞]` is finite for
each of the finitely many `p` with `Sha[p] ≠ 0`. The ultraproduct is blind to
that finite set, which must be handled separately. For finiteness this need
not be effective, but identifying the set in practice needs per-prime
criteria.

**Lemma 2 (finitely many rational values).** Let `x_p ∈ Q_p` for each `p`.
Suppose that for every non-principal `U` there is `b_U ∈ Q` with
`{p : x_p = b_U} ∈ U`. Then there is a finite set `B ⊂ Q` with `x_p ∈ B` for
all but finitely many `p`.

*Proof.* Let `T = {p : x_p ∉ Q}`. If `T` were infinite, a non-principal
`U ∋ T` would contain no set `{x_p = b}`. So `T` is finite. For `b ∈ Q` let
`S_b = {p : x_p = b}`; these sets are disjoint. Suppose infinitely many
`S_b` are non-empty. Choose one prime from each to get an infinite set `P'`,
and take a non-principal `U ∋ P'`. Some `S_b ∈ U`, so `S_b ∩ P' ∈ U`. But
`S_b ∩ P'` is a single prime, which contradicts `U` being non-principal. ∎

So **Weak GRS ⟺ for every `U`, the ultraproduct `[𝔅_p]_U ∈ ∏_U Q_p` is a
"standard" rational number.** Why this helps:

- **Łoś transfer of per-prime theory.** For non-CM `E`, the hypotheses of
  Kolyvagin's structure theorem and of Skinner–Urban hold for `U`-most `p`
  (Serre's open image theorem). The whole per-prime theory therefore becomes
  a single theory over `F_U`, or over `∏_U Z_p`, in characteristic 0.
- **Uniformity becomes rigidity.** The question "is this non-standard number
  built from standard global data actually standard?" is answered in model
  theory by *uniform definability*. Concretely: if `L_p^*` and `E_p·Reg_p` are
  computed by one algebraic recipe uniform in `p` (modular symbols, global
  points, Coleman integrals), and their ratio satisfies an identity with no `p`
  in it, standardness follows. A motivic period identity (Note 2) is exactly
  that kind of identity.

**Warnings.**
- The ramification of `E[p]` at `p` moves with `p`. Galois cohomology does
  not commute with ultraproducts in general, so we use ultraproducts only at
  the level of finished Selmer groups and of numbers.
- The abstract isomorphism `∏_U F̄_p ≅ C` (a Lefschetz principle) depends on
  the axiom of choice and respects no topology. **Nothing here relies on it.**
  A *canonical* structure-preserving comparison would itself be a major new
  theorem.

## 6. Arbitrary rank

Everything above works for any `r`, using `L_p^(r)/r!` and the `r × r`
regulator `Reg_p`. For a curve whose rank `r` is known, **finiteness of
`Sha` reduces to three inputs, none tied to `r = 1`:**

- **(a)** per-prime main conjectures, ordinary and signed (Skinner–Urban;
  Wan; Sprung; Castella–Çiperiani–Skinner–Sprung), including at `p = 2, 3`
  (Route A of `PLAN.md`);
- **(b)** `Reg_p ≠ 0` at every `p`. This is a `p`-adic Schanuel-type
  statement and open, but under GRS it is equivalent to `L_p^(r) ≠ 0`;
- **(c)** Weak GRS.

This is the note's answer to obstacle 1. The arbitrary-rank control theorem
is **not** a rank-`r` Euler system. It is the per-prime theory that already
exists, plus one comparison-of-realizations statement. The difficulty moves
out of Iwasawa theory and into the theory of periods.

What this does **not** give is `rank E(Q) ≥ ord L(E,s)`, the existence of
points. That is the subject of [Note 4](04_multiplace_gross_zagier.md), and
[Note 2](02_motivic_transfer.md) handles `ord L_p = ord L`.

## 7. Novelty and nearest literature

- **Nearest work.** Stein–Wuthrich (Math. Comp. 2013) bound `Sha[p]` prime by
  prime for rank ≥ 2 curves using Iwasawa theory. arXiv:2609.08431 gives the
  rank-2 CM unit criterion. The uniformity mechanisms in rank 0 and 1 are
  classical.
- **What we did not find stated anywhere:** the explicit reduction of rank-2
  finiteness to a *p-independence of p-adic BSD quotients* statement; the
  Borel–Cantelli no-go showing that unit criteria are silent at infinitely
  many primes; and the ultraproduct form (Proposition 1 and Lemma 2) as a
  bridge from per-prime theory to uniformity.
- **Main risk.** GRS may be as hard as BSD. The only route to it we can see
  is motivic (Note 2). Even so, GRS isolates exactly what is missing, and it
  is falsifiable (below).

## 8. Kill test

This is experiment **E1** in [`EXPERIMENTS.md`](EXPERIMENTS.md). Compute
`𝔅_p` for `389a1`, `571b1` and `5077a1` at every good ordinary `p < 200`.

- GRS predicts `𝔅_p = 1` exactly in all three cases. For each curve,
  `#Sha = ∏c_v = |tors| = 1` in our certified table.
- **Kill criterion:** any `p` where `𝔅_p` differs from `1` beyond the
  certified `p`-adic precision. That would also refute `p`-adic BSD, which
  makes it unlikely, but the lab's independence from the Sage/PARI lineage is
  exactly what makes the check worth running.

### 8.1 Correction (2026-09-13, orchestrator literature check)

This experiment has largely been run already. Stein–Wuthrich, *Algorithms
for the arithmetic of elliptic curves using Iwasawa theory*, Math. Comp. 82
(2013) 1757–1792, **Theorem 12.3**: for `389a1` and all 5,005 good ordinary
`p < 48,859` except `p = 16,231`, `Sha(E/Q)[p] = 0`, and the `p`-adic BSD
conjectural order of `Sha` is `≡ 1 (mod p)`. In our notation that is
`𝔅_p ≡ 1 mod p` at those primes. **Theorem 1.1** gives the same for
1,534,422 pairs `(E, p)`: rank ≥ 2, conductor ≤ 30,000, `5 ≤ p < 1000`.
Remark 12.4: at `p = 16,231`, `ord_p Reg_p = 3`, not 2.

Consequences, stated plainly:

1. E1 at `p < 200` **cannot be a historic first**. It is an independent
   reproduction of known numerics by non-Sage code. That still has value:
   it validates our p-adic machinery against a published external oracle.
   A mismatch would almost certainly mean a bug in our code, not a
   counterexample to p-adic BSD.
2. **No finite computation tests GRS's new content.** GRS adds uniformity
   in `p` to p-adic BSD, and a table over finitely many primes is consistent
   with GRS true or false. So E1 can only ever kill GRS through a
   p-adic-BSD failure, which Stein–Wuthrich's data make very unlikely.
3. The normalisation to use (Balakrishnan–Müller–Stein, Math. Comp. 85
   (2016), §1): `T = κ(γ)^{s-1} − 1` and `Reg_γ = Reg_p / (log_p κ(γ))^r`.
   The leading `T`-coefficient is compared with `E_p · Reg_γ · #Sha ·
   ∏c_v / |tors|²`.
4. `p = 16,231` for `389a1` is a free extra anchor: the pipeline must find
   `ord_p Reg_p = 3` there, if it can reach that prime.
