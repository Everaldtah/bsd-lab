"""Check the j = 0 (c4 = 0) minimalisation fix against the LMFDB oracle.

`EllipticCurve.minimal_model()` used to raise on every c4 = 0 input: the
p-adic valuation of 0 is infinity, and the Laska-Kraus-Connell scaling took
``min()`` over ``float('inf') // n`` -- which is NaN in CPython -- poisoning
the exponent search (see `_lkc_scaling` in bsdlab/curve.py).  j = 0 curves
were therefore silently dropped by tools/enumerate_curves.py's failure guard.

This script is the authoritative regression check for that fix: it walks
every curve in data/lmfdb_oracle.json (all LMFDB curves of conductor <= 400;
oracle rows ARE minimal models), keeps the ones with j-invariant 0 (computed
here from the a-invariants via bsdlab -- j = 0 <=> c4 = 0, so no j field is
needed in the oracle), pushes each through minimal_model() and demands the
exact oracle a-invariants back.  Passing through a minimal model must be a
no-op.  It also verifies the construction properties the task requires:
integrality (implicit in the a-invariant type), minimality at every prime,
and preservation of the j-invariant (Q-isomorphism, since j = 0 is rigid).

Usage:
    python tools/check_j0_fix.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bsdlab.curve import EllipticCurve  # noqa: E402

ORACLE = ROOT / 'data' / 'lmfdb_oracle.json'

MAX_PRINTED_FAILURES = 10


def main() -> int:
    blob = json.loads(ORACLE.read_text())
    rows = blob.get('curves')
    if not rows:
        sys.exit(f'{ORACLE} holds no completed oracle (curves key missing or '
                 f'empty; an interrupted fetch leaves a checkpoint only). '
                 f'Run tools/fetch_lmfdb_oracle.py to finish it.')

    total = passed = 0
    failures: list[dict] = []
    for row in rows:
        E = EllipticCurve(*row['ainvs'])
        if E.j_invariant != 0:
            continue  # not a j = 0 (c4 = 0) curve
        total += 1
        label = row.get('lmfdb_label') or '?'
        try:
            M = E.minimal_model()
        except Exception as exc:
            failures.append({'label': label, 'ainvs': row['ainvs'],
                             'why': f'minimal_model() raised {exc.__class__.__name__}: {exc}'})
            continue
        # Primary requirement: the oracle's minimal a-invariants, exactly.
        if list(M.ainvs) != list(row['ainvs']):
            failures.append({'label': label, 'ainvs': row['ainvs'],
                             'why': f'got {list(M.ainvs)}'})
            continue
        # Construction properties: still j = 0 (Q-isomorphic) and minimal.
        why = None
        if M.j_invariant != 0:
            why = f'j-invariant changed to {M.j_invariant}'
        elif not M.is_minimal():
            why = 'returned model is not minimal'
        if why:
            failures.append({'label': label, 'ainvs': row['ainvs'], 'why': why})
            continue
        passed += 1

    print(f'j = 0 curves in oracle (conductor <= {blob.get("meta", {}).get("bound", 400)}): '
          f'{total}')
    print(f'PASS: {passed}   FAIL: {len(failures)}')
    for f in failures[:MAX_PRINTED_FAILURES]:
        print(f"  {f['label']} ainvs={f['ainvs']}: {f['why']}")
    if len(failures) > MAX_PRINTED_FAILURES:
        print(f'  ... and {len(failures) - MAX_PRINTED_FAILURES} more')
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
