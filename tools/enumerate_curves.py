"""Enumerate elliptic curves over Q from first principles.

Everything here is derived from the curve equations via ``bsdlab``; this tool
never reads any reference table (in particular not data/lmfdb_reference.json).
The pipeline is:

    enumerate_models  -> every integral model in a bounded box
    reduce_to_minimal -> one minimal model per curve (with a j = 0 workaround)
    dedupe            -> one representative per Q-isomorphism class
    sweep             -> conductor + Tamagawa numbers per curve, as JSON rows

Run as::

    python tools/enumerate_curves.py --a4 20 --a6 40 --out data/curves.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

# Allow running as a plain script from any working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bsdlab import reduction
from bsdlab.curve import EllipticCurve

# bsdlab.curve.minimal_model raises ValueError ("cannot convert float NaN to
# integer") on every curve with c4 == 0 (j = 0); it can also raise
# ArithmeticError when Laska-Kraus-Connell scaling strands without an integral
# model. Both are worked around here, not fixed in bsdlab.
_MINIMAL_FAILURES = (ArithmeticError, ValueError, ZeroDivisionError,
                     OverflowError, TypeError)


def enumerate_models(a4_bound: int, a6_bound: int):
    """Yield every model y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6 with
    a1 in {0,1}, a2 in {-1,0,1}, a3 in {0,1} and a4, a6 in [-bound, bound].

    Every curve over Q has such a model, so this box is exhaustive up to the
    bounds. Singular cubics (discriminant 0) are skipped: they are not curves.
    """
    for a1 in (0, 1):
        for a2 in (-1, 0, 1):
            for a3 in (0, 1):
                for a4 in range(-a4_bound, a4_bound + 1):
                    for a6 in range(-a6_bound, a6_bound + 1):
                        E = EllipticCurve(a1, a2, a3, a4, a6)
                        if E.discriminant == 0:
                            continue
                        yield E


def reduce_to_minimal(E: EllipticCurve):
    """Return (minimal model, ok).

    ok is False when the known bsdlab j = 0 bug (c4 == 0) — or any other
    minimal_model failure — forced us to fall back to E's own model, which may
    be non-minimal. The caller records that in the row's "minimal" flag rather
    than crashing the sweep.
    """
    try:
        M = E.minimal_model()
    except _MINIMAL_FAILURES:
        return E, False
    if any(not isinstance(v, int) for v in M.ainvs) or M.discriminant == 0:
        # Inf/NaN contamination from the same bug, in the variant where the
        # call returns garbage instead of raising.
        return E, False
    return M, True


def dedupe(curves) -> list:
    """One representative per Q-isomorphism class.

    Keyed on (c4, c6) of the minimal model: two curves over Q are isomorphic
    exactly when their minimal models share (c4, c6), whereas the a-invariants
    themselves vary with the chosen model.
    """
    # Two Q-curves are isomorphic iff their minimal models have the same (c4, c6).
    seen = {}
    for E in curves:
        M, _ok = reduce_to_minimal(E)
        key = (M.c4, M.c6)
        if key not in seen:
            seen[key] = M
    return list(seen.values())


def _row(E: EllipticCurve) -> dict:
    """The JSON row for one (minimal or fallback) model, all values plain ints."""
    N = reduction.conductor(E)
    tam = reduction.tamagawa_numbers(E)
    j = E.j_invariant
    try:
        minimal = E.is_minimal()
    except _MINIMAL_FAILURES:
        # j = 0 workaround: minimality cannot be certified, and the emitted
        # model is the unreduced fallback.
        minimal = False
    return {
        "ainvs": [int(v) for v in E.ainvs],
        "c4": int(E.c4),
        "c6": int(E.c6),
        "discriminant": int(E.discriminant),
        "j_invariant": "%d/%d" % (j.numerator, j.denominator),
        "conductor": int(N),
        "tamagawa": {str(p): int(c) for p, c in sorted(tam.items())},
        "minimal": bool(minimal),
    }


def sweep(a4_bound: int, a6_bound: int, max_conductor=None) -> list:
    """Deduped minimal curves in the box with local invariants attached.

    Rows whose conductor computation raises are skipped and counted, with a
    summary line per failure reason; they are never silently dropped.
    """
    models = list(enumerate_models(a4_bound, a6_bound))
    print("enumerated %d models" % len(models))
    curves = dedupe(models)
    print("after dedupe: %d curves" % len(curves))

    rows = []
    skipped = Counter()
    for E in curves:
        try:
            row = _row(E)
        except Exception as ex:  # counted below, never silent
            skipped[type(ex).__name__] += 1
            continue
        if max_conductor is not None and row["conductor"] > max_conductor:
            continue
        rows.append(row)

    if skipped:
        total = sum(skipped.values())
        for reason in sorted(skipped):
            print("skipped %d curves: %s: %d"
                  % (total, reason, skipped[reason]))

    # Canonical order so runs are byte-reproducible.
    rows.sort(key=lambda r: (r["conductor"], r["c4"], r["c6"]))
    print("kept %d rows%s" % (
        len(rows),
        " (max conductor %d)" % max_conductor if max_conductor is not None
        else ""))
    return rows


def print_histogram(rows) -> None:
    """A compact conductor histogram summary on stdout."""
    if not rows:
        print("conductor histogram: (no rows)")
        return
    counts = Counter(r["conductor"] for r in rows)
    buckets = [(1, 10), (11, 50), (51, 100), (101, 500), (501, 1000),
               (1001, 5000), (5001, 10000), (10001, float("inf"))]
    print("conductor histogram:")
    bars = [sum(c for N, c in counts.items() if lo <= N <= hi)
            for lo, hi in buckets]
    scale = max(1, max(bars) // 50)
    for (lo, hi), n in zip(buckets, bars):
        if n:
            label = ">10000" if hi == float("inf") else \
                    ("%d-%d" % (lo, hi) if lo > 1 else "1-%d" % hi)
            print("  %-10s %6d %s" % (label, n, "#" * max(1, n // scale)))
    print("  conductor range: %d .. %d"
          % (min(counts), max(counts)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Enumerate elliptic curves over Q from the equations up.")
    ap.add_argument("--a4", type=int, required=True,
                    help="bound on |a4| (inclusive)")
    ap.add_argument("--a6", type=int, required=True,
                    help="bound on |a6| (inclusive)")
    ap.add_argument("--out", type=str, default=None,
                    help="write the JSON dataset here")
    ap.add_argument("--max-conductor", type=int, default=None,
                    help="drop rows whose conductor exceeds this")
    args = ap.parse_args(argv)

    rows = sweep(args.a4, args.a6, args.max_conductor)
    print_histogram(rows)

    if args.out:
        path = Path(args.out)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="\n") as f:
            json.dump(rows, f, indent=2)
            f.write("\n")
        print("wrote %d rows to %s" % (len(rows), path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
