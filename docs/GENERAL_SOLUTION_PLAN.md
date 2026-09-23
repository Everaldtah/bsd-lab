# The road to a general solution of BSD

*Master plan, 2026-09-23. It replaces the sequencing in `PLAN.md` §8 and
keeps that document's decomposition. Read with
[`research/06_route_b_exhaustive.md`](research/06_route_b_exhaustive.md),
which is the evidence behind §3.*

> **Honest framing.** No accepted program for the general case exists. This
> plan does not claim one. It turns the problem into a finite list of
> precisely stated missing theorems (MT1–MT4). It says which are within
> reach, which need new mathematics, and what this lab computes toward each.
> Anything in this repository that looks like a proof of an open case is a
> bug.

## 1. Where things stand

| component | status | reference |
|---|---|---|
| Modularity of `E/Q` | **proved** | Wiles; Breuil–Conrad–Diamond–Taylor |
| Parity of `p^∞`-Selmer rank | **proved** | Nekovář; Dokchitser–Dokchitser |
| Rank equality + Sha finite, `ord L ≤ 1` | **proved** | Gross–Zagier + Kolyvagin |
| `p`-part of exact formula, `ord L ≤ 1` | **partial** (`p ≥ 5` ordinary, big image, …) | Kato; Skinner–Urban; Jetchev–Skinner–Wan; … |
| Weak BSD for *individual* curves with `ord L ≤ 3` | **provable case by case** (descent + BGZ trick) | e.g. `389a` rank 2, `5077a` rank 3 (BGZ 1985) |
| Selmer rank 2 from one `p` | **proved under hypotheses** (`κ_p ≠ 0`) | Castella–Hsieh; Castella |
| Sha finite for *any* rank ≥ 2 curve | **open** | — |
| Any curve proven to have `ord L ≥ 4` | **open** | — |
| Rank equality for all `E` with `ord L ≥ 2` | **open** | — |

The lab's own verified state is in `FINDINGS.md`. In short:

- The full invariant pipeline works on 14 reference curves, with rigorous
  error bounds.
- 559 curves match the LMFDB on local data, and 39,790 curves have been
  enumerated.
- The `p`-adic stack (modular symbols, Kedlaya/E2, MTT, MST heights) is
  anchored to Sage.
- The `𝔅_p = 1` table covers four curves.
- RB1: 43/43 Heegner points on the rank-2 and rank-3 curves are torsion.

## 2. The theorem-level target, decomposed

Full BSD equals the conjunction of four missing theorems (stated in full in
Note 06 §4):

- **MT1 — uniformity in `p`:** `Sha(E)[p] = 0` for almost all `p` when
  `ord L ≥ 2`.
- **MT2 — complex ↔ `p`-adic:** `ord L = ord L_p`, and `L''(E,1) ≠ 0 ⇒
  κ_p(E) ≠ 0`.
- **MT3 — a higher Gross–Zagier over `Q`:** `L^(r)(1)/r! ≈ det⟨P_i,P_j⟩`
  for special-cycle classes.
- **MT4 — Route A:** the `p`-part of the formula at every `p`
  (2, 3, supersingular, additive, reducible).

Dependencies:

```
MT4 ───────────────────────────────────────────────┐
MT1 (one curve) ──► R2: Sha finite for one rank-2 E │
MT1 (uniform) + MT2 ──► rank equality, ord L = 2    ├──► full BSD (R3)
MT3 ──► MT1 + MT2 in every rank; ord L ≥ 4 wall ────┘
```

## 3. Route B verdicts (from Note 06)

| candidate | verdict | why |
|---|---|---|
| Heegner points / twists | **killed** as a point source | No-Go Lemma; RB1 43/43 torsion |
| Stark–Heegner / Darmon points | **killed** for Route B | governed by `L'(E/K,ψ,1)`: same lemma |
| Beilinson–Flach elements | **dormant** | no known link to `L''(E,1)` |
| Kolyvagin systems | **alive** | sees Selmer rank 2 per `p`; needs MT1 |
| Diagonal cycles / generalised Kato | **alive** (strongest) | `κ_p ≠ 0 ⇒ dim Sel = 2`; needs MT1 + MT2 |
| Anticyclotomic / derived heights | **alive** per `p` | Castella–Hsieh, Banwait; needs MT1 |
| Higher GZ (shtukas) | **the right shape** | exists only over function fields → MT3 |

**The strategic consequence.** Route B is not "find rank-2 points". The
points are already known for the reference curves. Route B is: **prove the
Selmer group is no bigger than the points, at every prime simultaneously.**
All the live machinery attacks that one prime at a time.

## 4. The rungs, re-scoped

- **R1** — BSD for all curves of analytic rank ≤ 1. This is MT4 only. It is
  the profession's tractable flank.
- **R2** — `Sha(E/Q)` finite for one rank-2 curve, with `389a1` as the
  target. This is MT1 for one curve. **It is the next historic milestone,
  and the lab is aimed at it.**
- **R2+** — rank equality for all curves of analytic rank 2. This needs
  MT1 uniformly in `E`, plus MT2.
- **R3** — the general case. This needs MT3, or an equivalent new object,
  plus everything above.

## 5. The lab program from here

Each item is either an unconditional deliverable or a sharper test of one
missing theorem. They are ordered by value per unit of work.

- **L7. Kolyvagin classes mod `p` for 389a1 (MT1, Kolyvagin form).** This
  was experiment E2 in `research/EXPERIMENTS.md`.
  - Choose `K` with `rank E^K = 1`.
  - Compute `κ_{ℓ₁ℓ₂}` mod `p` via Brandt modules and Gross points. This
    is how Jetchev–Lauter–Stein did it for small `p`.
  - Do this for every `p` in `[5, 100]`.
  - Deliverable: the first table of `ν_p = 2` witnesses across many `p`.
    That is empirical evidence *for* the uniform statement MT1 needs, in a
    form a proof would have to explain.
  - Kill condition: a `p` with `Sha(E/K)[p] = 0` known, yet every
    `κ_{ℓ₁ℓ₂} ≡ 0` for all `ℓᵢ` below a large bound. That would refute the
    naive anti-concentration model (Note 3).
- **L8. Generalised Kato classes for 389a1 (MT2).**
  - Compute `κ_p(E)` numerically via Darmon–Lauder–Rotger `p`-adic iterated
    integrals, using Lauder's overconvergent projection, for `p = 5, 7, 11`.
  - Check `κ_p ≠ 0` and `loc_p(κ_p) = 0` (it lies in the strict Selmer
    group, per Castella).
  - This builds on the existing `padicl`/`modsym` stack.
- **L2′. Full 2-descent without rational 2-torsion.**
  - This proves `rank(389a1) = rank(571b1) = 2` and `rank(5077a1) ≤ 3`
    in-lab. Weak BSD for those curves then becomes a lab theorem, not a
    citation.
  - It needs cubic-field arithmetic (class groups and units of `Q(θ)`).
- **L4. Scale certification** from 14 to about 500 curves
  (`scratch/certify_scale.py`, staged).
- **E3. Geodesic barycenters (Note 4).** This is the one experiment aimed
  at MT3 on the archimedean side.
- **Housekeeping:**
  - finish the twist census for 5077a1;
  - wire `rootnumbers` into `lseries`;
  - fix the float round-trip in `tools/fetch_lmfdb.py`.

## 6. The theory program (what a mathematician would have to prove)

**T1 (toward MT1, one curve).** Prove, for `389a1` and a fixed `K`, that
the mod-`p` Kolyvagin system has a non-zero class with `ν = 2` for all
`p > p₀`.

- The most promising line of attack is an equidistribution-mod-`p`
  (large-sieve) statement for Heegner points of conductor `ℓ₁ℓ₂` on the
  supersingular locus. See Note 3.
- Note 1's Borel–Cantelli no-go shows that independent per-prime tests
  cannot suffice. The argument must use one global object.

**T2 (toward MT2).** Show that the Castella–Hsieh leading term of the
anticyclotomic `L_p` equals a period of the same motivic extension as
`L''(E,1)`. See Note 2 on motivic transfer.

**T3 (toward MT3).**

- Construct, over `Spec Z`, a replacement for "shtukas with `r` legs" whose
  self-intersection computes `L^(r)`.
- There is no existing candidate. Arithmetic intersection theory on
  products of Shimura varieties (higher-dimensional Gross–Zagier, the
  Kudla program) is the nearest technology, but it is still first-order.
- This is the Wiles-scale new idea the problem needs. We record it without
  pretending to supply it.

**T4 (MT4).**

- Supersingular `±`-main conjectures at every sign.
- `p = 2, 3` Iwasawa theory.
- Additive and ramified-level main conjectures.
- The Eisenstein / reducible case.

These are active subfields with no conceptual wall.

## 7. Risk statement

- R1 is realistic for the profession over years.
- R2 is realistic only if T1 is found. The lab can supply the evidence
  base (L7) that makes T1 a precise conjecture about explicit objects.
- R3 needs T3. Nothing in the literature or in this repository comes close
  to T3.

This plan narrows the general case to named statements and removes dead
branches: two Route-B candidates are killed by theorem, and one is dormant.
It does not remove the need for new mathematics.
