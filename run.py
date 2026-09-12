"""Batch driver: run the full BSD pipeline over a set of curves and cross-check
against the LMFDB reference.

This is the *user-facing* entry point.  It loads reference invariants from
``data/lmfdb_reference.json`` (produced by ``tools/fetch_lmfdb.py``), recomputes
every quantity from scratch via ``bsdlab.bsd.analyse``, and reports how the
recomputation agrees with the reference.  The LMFDB numbers are used only as an
independent oracle; ``bsdlab`` never sees them.

Usage:
    python run.py --labels 11a1 37a1 389a1
    python run.py --max-conductor 200
    python run.py --labels 571a1 681b1 --prec 80
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

import mpmath

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from bsdlab.curve import EllipticCurve            # noqa: E402
from bsdlab.bsd import analyse                     # noqa: E402


def _load_reference(path: Path) -> dict:
    """Map Cremona label -> reference row."""
    rows = json.loads(path.read_text())
    return {r['Clabel']: r for r in rows}


def _nstr(v, n: int = 15) -> str:
    if v is None:
        return '-'
    if isinstance(v, (int,)) and not isinstance(v, bool):
        return str(v)
    if isinstance(v, mpmath.mpf):
        return mpmath.nstr(v, n)
    return str(v)


def _run_one(ref: dict, prec: int, height_bound: int, verbose: bool) -> dict:
    label = ref['Clabel']
    E = EllipticCurve.from_list(ref['ainvs'])
    d = analyse(E, label=label, prec=prec, height_bound=height_bound)

    # Cross-checks against the independent LMFDB reference.
    rank_match = (d.analytic_rank is not None
                  and d.analytic_rank == ref.get('analytic_rank'))
    sha_match = (d.sha_rounded is not None and ref.get('sha') is not None
                 and d.sha_rounded == ref['sha'])

    result = {
        'label': label,
        'ainvs': list(E.ainvs),
        'conductor': d.conductor,
        'conductor_ref': ref.get('conductor'),
        'root_number': d.root_number,
        'analytic_rank': d.analytic_rank,
        'analytic_rank_ref': ref.get('analytic_rank'),
        'rank_match': rank_match,
        'leading_coefficient': (_nstr(d.leading_coefficient)
                                if d.leading_coefficient is not None else None),
        'real_period': _nstr(d.real_period),
        'rank_lower_bound': d.rank_lower_bound,
        'rank_certified': d.rank_certified,
        'regulator': _nstr(d.regulator),
        'regulator_ref': ref.get('regulator'),
        'torsion_order': d.torsion_order,
        'tamagawa': d.tamagawa,
        'tamagawa_product': d.tamagawa_product,
        'sha_analytic': _nstr(d.sha_analytic),
        'sha_rounded': d.sha_rounded,
        'sha_ref': ref.get('sha'),
        'sha_match': sha_match,
        'sha_integrality_error': _nstr(d.sha_integrality_error, 8),
        'sha_is_square': d.sha_is_square,
        'consistent': d.consistent,
        'notes': d.notes,
    }

    if verbose:
        print(f'\n=== {label}  ainvs={list(E.ainvs)} ===')
        print(f'  conductor  {d.conductor}  (ref {ref.get("conductor")})')
        print(f'  root no.   {d.root_number}')
        print(f'  a-rank     {d.analytic_rank}  (ref {ref.get("analytic_rank")})'
              f'  certified={d.rank_certified}')
        print(f'  L^r(1)/r!  {_nstr(d.leading_coefficient)}')
        print(f'  Omega      {_nstr(d.real_period)}')
        print(f'  rank lb    {d.rank_lower_bound}   gens={d.generators}')
        print(f'  Reg        {_nstr(d.regulator)}  (ref {ref.get("regulator")})')
        print(f'  tors       {d.torsion_order}   prod c_p = {d.tamagawa_product}')
        print(f'  #Sha_an    {_nstr(d.sha_analytic)}  -> {d.sha_rounded}'
              f'  (ref {ref.get("sha")})  square={d.sha_is_square}')
        for note in d.notes:
            print(f'  note: {note}')
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--labels', nargs='*', default=[],
                    help='Cremona labels, e.g. 11a1 37a1 571a1')
    ap.add_argument('--max-conductor', type=int, default=0,
                    help='run every reference curve of conductor <= MAX')
    ap.add_argument('--ref', default='data/lmfdb_reference.json')
    ap.add_argument('--out', default='data/results.json')
    ap.add_argument('--prec', type=int, default=60)
    ap.add_argument('--height-bound', type=int, default=200,
                    help='Neron-Tate height bound for the point search')
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args()

    refmap = _load_reference(Path(args.ref))
    if args.labels:
        selected = []
        for lab in args.labels:
            if lab in refmap:
                selected.append(refmap[lab])
            else:
                print(f'warning: {lab} not in reference, skipping',
                      file=sys.stderr)
    elif args.max_conductor:
        selected = [r for r in refmap.values()
                    if r.get('conductor', 10**9) <= args.max_conductor]
        selected.sort(key=lambda r: (r['conductor'], r['Clabel']))
    else:
        ap.error('give --labels or --max-conductor')

    results = []
    for ref in selected:
        try:
            results.append(_run_one(ref, args.prec, args.height_bound,
                                    args.verbose))
        except Exception as exc:            # noqa: BLE001 - report, don't die
            print(f'error on {ref["Clabel"]}: {exc!r}', file=sys.stderr)
            results.append({'label': ref['Clabel'], 'error': repr(exc)})

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=1))

    # Summary table.
    print('\n' + '=' * 100)
    hdr = (f'{"label":<8}{"cond":>6}{"a-rk":>6}{"rk?":>5}{"Sha":>8}'
           f'{"Sha?":>6}{"square":>8}{"ok":>5}  notes')
    print(hdr)
    print('-' * 100)
    n_ok = 0
    for r in results:
        if 'error' in r:
            print(f'{r["label"]:<8}  ERROR {r["error"]}')
            continue
        ok = 'yes' if r['consistent'] else 'no '
        if r['consistent']:
            n_ok += 1
        note = (r['notes'][0][:40] if r['notes'] else '')
        print(f'{r["label"]:<8}{r["conductor"] or 0:>6}'
              f'{r["analytic_rank"] if r["analytic_rank"] is not None else -1:>6}'
              f'{("y" if r["rank_match"] else "n"):>5}'
              f'{r["sha_rounded"] if r["sha_rounded"] is not None else -1:>8}'
              f'{("y" if r["sha_match"] else "n"):>6}'
              f'{str(r["sha_is_square"]):>8}{ok:>5}  {note}')
    print('-' * 100)
    print(f'{n_ok}/{len(results)} curves fully consistent '
          f'(rank + Sha integrality + square). Wrote {args.out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
