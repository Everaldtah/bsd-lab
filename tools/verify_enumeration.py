"""Cross-validate the box-bounded enumeration against the LMFDB oracle.

Joins data/curves.json against data/lmfdb_oracle.json (built by
tools/fetch_lmfdb_oracle.py) on ainvs and checks, for every curve present
on both sides:

  * the conductor agrees;
  * the Tamagawa product agrees (product of the per-prime values in the
    enumeration row vs the oracle's number).

Curves in the oracle but NOT in the enumeration are reported as a count
only: the enumeration is box-bounded, so misses are expected and are not
failures.  The exit code is 0 iff the matched set has zero mismatches.

stdlib only (json/math/sys); nothing from bsdlab is needed -- this script
compares two data files, it does not recompute anything.

Usage:
    python tools/verify_enumeration.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CURVES = ROOT / 'data' / 'curves.json'
ORACLE = ROOT / 'data' / 'lmfdb_oracle.json'

MAX_PRINTED_MISMATCHES = 20


def norm_ainvs(value) -> tuple[int, ...]:
    """Normalise an ainvs representation to a tuple of ints.

    The API and both data files store lists of ints, but ainvs also has a
    string life on the LMFDB web forms ("a1;a2;a3;a4;a6"), so strings with
    ';' or ',' separators are accepted too.
    """
    if isinstance(value, str):
        parts = value.replace(',', ';').split(';')
        return tuple(int(p.strip()) for p in parts)
    if isinstance(value, (list, tuple)):
        return tuple(int(x) for x in value)
    raise ValueError(f'unrecognised ainvs {value!r}')


def enum_tamagawa_product(row: dict) -> int | None:
    tam = row.get('tamagawa') or {}
    vals = [v for v in tam.values() if v is not None]
    return math.prod(vals) if vals else None


def oracle_tamagawa_product(row: dict) -> int | None:
    if row.get('tamagawa_product') is not None:
        return row['tamagawa_product']
    tam = row.get('tamagawa') or {}
    vals = [v for v in tam.values() if v is not None]
    return math.prod(vals) if vals else None


def main() -> int:
    enum_rows = json.loads(CURVES.read_text())
    oracle_blob = json.loads(ORACLE.read_text())
    oracle_rows = oracle_blob['curves']
    bound = oracle_blob.get('meta', {}).get('bound', 400)

    # Index the enumeration by ainvs (keys are unique: checked offline, and
    # a duplicate here is reported rather than silently collapsed).
    enum: dict[tuple[int, ...], dict] = {}
    dupes = 0
    for row in enum_rows:
        key = norm_ainvs(row['ainvs'])
        if key in enum:
            dupes += 1
        enum[key] = row
    if dupes:
        print(f'NOTE: {dupes} duplicate ainvs in curves.json '
              f'(last occurrence kept)')

    matches = 0
    mismatches: list[dict] = []
    unverifiable = 0
    misses: list[dict] = []
    for orc in oracle_rows:
        key = norm_ainvs(orc['ainvs'])
        row = enum.get(key)
        if row is None:
            misses.append(orc)
            continue
        matches += 1
        fields = {}
        if row['conductor'] != orc['conductor']:
            fields['conductor'] = (row['conductor'], orc['conductor'])
        e_prod = enum_tamagawa_product(row)
        o_prod = oracle_tamagawa_product(orc)
        if e_prod is None or o_prod is None:
            # Nothing to compare (e.g. a curve with no localdata rows);
            # surfaced, never silently passed.
            unverifiable += 1
        elif e_prod != o_prod:
            fields['tamagawa_product'] = (e_prod, o_prod,
                                          row.get('tamagawa'),
                                          orc.get('tamagawa'))
        if fields:
            mismatches.append({
                'lmfdb_label': orc.get('lmfdb_label'),
                'lmfdb_iso': orc.get('lmfdb_iso'),
                'ainvs': list(key),
                'fields': fields,
            })

    # Reverse direction, informational only: enumeration rows inside the
    # oracle's conductor range that the oracle does not list.  Non-minimal
    # models are expected there (LMFDB lists minimal models only); minimal
    # ones would be surprising.
    enum_only_min, enum_only_nonmin = [], []
    oracle_keys = {norm_ainvs(o['ainvs']) for o in oracle_rows}
    for row in enum_rows:
        if row['conductor'] > bound or norm_ainvs(row['ainvs']) in oracle_keys:
            continue
        (enum_only_min if row.get('minimal', True) else enum_only_nonmin).append(row)

    print(f'oracle size (conductor <= {bound}): {len(oracle_rows)}')
    print(f'matches (joined on ainvs):         {matches}')
    print(f'mismatches:                        {len(mismatches)}')
    for m in mismatches[:MAX_PRINTED_MISMATCHES]:
        print(f"  {m['lmfdb_label']} (iso {m['lmfdb_iso']}) ainvs={m['ainvs']}")
        for name, detail in m['fields'].items():
            if name == 'conductor':
                print(f'    conductor: enumeration={detail[0]} '
                      f'oracle={detail[1]}')
            else:
                print(f'    {name}: enumeration={detail[0]} oracle={detail[1]} '
                      f'(enum per-prime={detail[2]}, oracle per-prime={detail[3]})')
    if len(mismatches) > MAX_PRINTED_MISMATCHES:
        print(f'  ... and {len(mismatches) - MAX_PRINTED_MISMATCHES} more')
    if unverifiable:
        print(f'unverifiable Tamagawa comparisons (null on one side): {unverifiable}')
    print(f'in oracle but not in enumeration (expected, box-bounded): {len(misses)}')
    if misses:
        sample = ', '.join(m.get('lmfdb_label', '?') for m in misses[:8])
        print(f'  e.g. {sample}')
    print(f'informational: enumeration rows N<={bound} absent from oracle: '
          f'{len(enum_only_min)} minimal, {len(enum_only_nonmin)} non-minimal')
    for row in enum_only_min[:MAX_PRINTED_MISMATCHES]:
        print(f"  MINIMAL row missing from oracle: ainvs={row['ainvs']} "
              f"N={row['conductor']}")

    return 1 if mismatches else 0


if __name__ == '__main__':
    raise SystemExit(main())
