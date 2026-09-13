# The barrier at rank 2

This document states, as precisely as we can, what is actually proven about BSD
and where each known method stops. It exists because "BSD is open" is true but
useless: the interesting fact is that BSD is *closed* in analytic rank 0 and 1,
by two different machines, and that both machines fail at rank 2 for reasons
that are structural rather than technical.

Nothing in this repository proves anything below. This is a map of the
literature, written so that the boundary is a list of statements.

## 1. What is unconditionally proven

**Theorem (Gross–Zagier 1986 + Kolyvagin 1990, with Waldspurger /
Bump–Friedberg–Hoffstein / Murty–Murty, and modularity by Wiles et al.).**
Let `E/Q` be an elliptic curve. If `ord_{s=1} L(E,s) <= 1`, then

- `rank E(Q) = ord_{s=1} L(E,s)`, and
- `Sha(E/Q)` is finite.

That is the entire unconditional rank statement. Note what it is *not*: it says
nothing when the order of vanishing is 2 or more, and it does not give the full
BSD formula — only the rank and the finiteness of Sha.

The mechanism, in one paragraph, because the shape of it is the whole point.
Pick an imaginary quadratic field `K` satisfying the Heegner hypothesis for the
conductor `N`. The modular parametrisation `X_0(N) -> E` carries CM points on
`X_0(N)` to a point `y_K ∈ E(K)`, the Heegner point. Gross–Zagier computes

```
L'(E/K, 1)  =  (explicit nonzero constant) · ĥ(y_K)
```

so `y_K` is non-torsion exactly when `L'(E/K,1) ≠ 0`. Kolyvagin then builds an
Euler system out of the Heegner points over ring class fields and shows that if
`y_K` is non-torsion, `rank E(K) = 1` and `Sha(E/K)` is finite, with the
Selmer group bounded by the index of `y_K`. The nonvanishing results supply a
`K` for which the twisted L-function behaves, which is what converts a statement
about `E/K` into one about `E/Q`.

## 2. The p-part of the formula

Rank and finiteness are not the same as the BSD *formula*. For that:

- **Kato (2004)**, via the Euler system of Beilinson–Kato elements, gives one
  divisibility of the Iwasawa Main Conjecture for `E/Q` at good ordinary `p >= 5`
  under hypotheses on the mod-`p` Galois representation. This already yields
  `Sha` finite and `rank = 0` when `L(E,1) ≠ 0`.
- **Skinner–Urban (2014)** supply the reverse divisibility under a further
  ramification hypothesis, giving the Main Conjecture and hence the `p`-part of
  the BSD formula in analytic rank 0.
- **W. Zhang (2014)** proves Kolyvagin's conjecture on the nonvanishing of the
  Heegner point Euler system under hypotheses, and **Jetchev–Skinner–Wan (2017)**
  give the `p`-part of the formula in analytic rank 1.
- **Bhargava–Skinner–Zhang (2014)** combine these with rank distribution results
  to show a positive proportion — at least 66% — of elliptic curves over `Q`,
  ordered by height, satisfy BSD.

Every one of these carries hypotheses: `p >= 5`, good ordinary reduction,
irreducibility or surjectivity of the residual representation, ramification
conditions. **The full BSD formula — all primes at once, including `p = 2, 3`,
supersingular and additive reduction — is not known for a single curve of rank 1,
let alone in general.** What is known is a collection of `p`-parts.

## 3. Where rank >= 2 breaks each method

### 3.1 There is no source of points

Gross–Zagier produces exactly one point. That is not a limitation of the proof;
it is what the construction is. A CM point on `X_0(N)` maps to a single point of
`E`, and the formula equates its height to a first derivative. When
`ord_{s=1} L(E,s) >= 2`, the Heegner point `y_K` is torsion, `ĥ(y_K) = 0`, and
the formula is the true but empty statement `0 = 0`.

No proven construction produces two independent points of infinite order from
L-function data. Darmon's Stark–Heegner points and their relatives are
conjectural. There is no known analogue of Gross–Zagier expressing `L''(E,1)` as
a pairing on a two-dimensional space of points.

### 3.2 Euler systems degenerate, not weaken

An Euler system bounds a Selmer group from above in terms of the index of its
bottom class. If the bottom class is torsion, the bound is not weak — it is
absent. Kolyvagin's argument gives `rank E(K) = 1` precisely because the bottom
class is a single non-torsion point; with a torsion bottom class the derived
classes vanish and no inequality survives.

A "rank 2 Euler system" would need a genuinely different object. The known
systems (Heegner, Beilinson–Kato, Beilinson–Flach) each rest on a single
cohomology class propagated up a tower, and the machine that converts them into
Selmer bounds consumes that class. This is why rank 2 is not an incremental
step past rank 1.

### 3.3 Parity is known; the rank is not

The parity side is in much better shape than the rank side. The `p`-parity
theorems (Nekovář; Dokchitser–Dokchitser 2010) give, for every `E/Q` and every
prime `p`,

```
corank_p Sel_{p^∞}(E/Q)  ≡  ord_{s=1} L(E,s)   (mod 2)
```

unconditionally. So if the analytic rank is even we know the Selmer corank is
even. This is compatible with rank 2 and with rank 0 plus a `Sha` of corank 2;
parity cannot distinguish them. Parity is a congruence, and a congruence mod 2
cannot certify a rank of 2.

### 3.4 Sha is not known to be finite for any rank >= 2 curve

Not one example. Finiteness of `Sha` is known exactly where Kolyvagin's argument
runs, which is exactly analytic rank <= 1. For `389a1` — the smallest conductor
with rank 2 — the finiteness of `Sha` is open.

> **Update (2026-09-13).** Per-prime rank-2 results now exist. Castella–Hsieh
> (Forum Math. Sigma 2022) show that a nonzero generalised Kato class gives a
> 2-dimensional `p`-adic Selmer group. Castella (arXiv:2204.09608) treats CM
> curves. arXiv:2609.08431 gives a unit criterion for `Sha[p^∞] = 0` via
> `L_p''` for rank-2 CM curves, outside an explicit set of primes. What
> remains open is `Sha[p] = 0` for **all but finitely many** `p` at once. See
> [`research/01_global_rational_shadow.md`](research/01_global_rational_shadow.md).

### 3.5 The algorithmic gap is the same gap

Descent computes Selmer groups, and for each `n`,

```
0 -> E(Q)/nE(Q) -> Sel^n(E/Q) -> Sha(E/Q)[n] -> 0
```

so a descent gives `rank <= dim Sel^n - dim(torsion part)`, with the slack
being exactly `Sha[n]`. Iterating descents separates the rank from `Sha` if and
only if `Sha` is finite. Hence there is no known algorithm that provably
terminates with the rank of an arbitrary `E/Q`; the termination proof would
require the finiteness of `Sha`, which is what rank >= 2 lacks. **The
computational barrier and the theoretical barrier are the same barrier.**

## 4. What this means for this repository

Applying the above to `docs/FINDINGS.md`:

| curve | analytic rank | status of our rank claim |
|---|---|---|
| 11a1/2/3, 14a1, 15a1, 32a3, 37b1, 571a1, 681b1 | 0 | **certified**: Kolyvagin applies, rank is 0 |
| 37a1, 681a1 | 1 | **certified**: Gross–Zagier + Kolyvagin, rank is 1 |
| 389a1, 571b1, 5077a1 | 2, 2, 3 | **not certified**: see below |

For the last three, what we actually have is:

1. a **lower** bound from the naive point search in `mordellweil.py`,
2. an **upper** bound from descent, when `bsdlab/descent.py` applies — and it
   bounds `rank + dim Sha[φ]`, not the rank,
3. **parity**, consistent in every case, which is genuinely unconditional,
4. the analytic rank, computed numerically, which nothing proves equals the
   algebraic rank at this order of vanishing.

Our `#Sha_an` for these curves is therefore conditional on **two** things, not
one: on strong BSD, and on the algebraic rank equalling the analytic rank. That
second condition is not a formality. It is `rank >= 2` BSD itself.

The `rank_certified` flag in `data/results.json` encodes exactly the boundary in
§1: it is `True` when the analytic rank is 0 or 1, and `False` otherwise. It is
not a statement about our numerical confidence, which is high for all fourteen
curves. It is a statement about which theorems exist.

## 5. The open boundary, as a list

1. Construct, for some `E/Q` of analytic rank 2, two independent points of
   infinite order from the L-function. No method is known.
2. Prove `Sha(E/Q)` finite for a single curve of rank >= 2.
3. Prove `rank E(Q) >= 2` follows from `ord_{s=1} L(E,s) >= 2`, for any curve.
4. Remove the ordinarity / `p >= 5` / residual-irreducibility hypotheses from the
   known `p`-part results, even in rank 0.
5. Give an algorithm that provably computes `rank E(Q)`, for every `E/Q`.

(1)–(3) are the same barrier viewed from three sides. (5) is that barrier
restated as a computation. (4) is a different and more tractable-looking problem
that has nonetheless resisted.

## References

- B. Gross, D. Zagier, *Heegner points and derivatives of L-series*, Invent. Math. 84 (1986).
- V. Kolyvagin, *Euler systems*, in The Grothendieck Festschrift II (1990).
- K. Kato, *p-adic Hodge theory and values of zeta functions of modular forms*, Astérisque 295 (2004).
- T. Dokchitser, V. Dokchitser, *On the Birch–Swinnerton-Dyer quotients modulo squares*, Ann. of Math. 172 (2010).
- C. Skinner, E. Urban, *The Iwasawa Main Conjectures for GL2*, Invent. Math. 195 (2014).
- W. Zhang, *Selmer groups and the indivisibility of Heegner points*, Camb. J. Math. 2 (2014).
- M. Bhargava, C. Skinner, W. Zhang, *A majority of elliptic curves over Q satisfy the Birch and Swinnerton-Dyer conjecture* (2014).
- D. Jetchev, C. Skinner, X. Wan, *The Birch–Swinnerton-Dyer formula for elliptic curves of analytic rank one*, Camb. J. Math. 5 (2017).
