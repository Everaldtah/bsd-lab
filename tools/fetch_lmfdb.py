"""Fetch reference BSD invariants from the LMFDB public API.

This is the *reference* side of the cross-validation: bsdlab computes every
quantity from scratch, and the numbers fetched here are what it is compared
against.  Nothing in ``bsdlab/`` imports this module -- the library itself has
no network dependency and no knowledge of the LMFDB.

Endpoint notes, established empirically (the API is undocumented in places):

  * The table is ``ec_curvedata`` and the working query form is exactly
    ``https://www.lmfdb.org/api/ec_curvedata/?<column>=<value>&_format=json``.
  * ``_fields``, ``_limit``, ``_offset`` and range syntax such as
    ``conductor=lt50`` are NOT accepted.  Each one returns an HTML error page,
    and repeated rejected requests trip a reCAPTCHA rate limiter that then
    refuses even valid queries for a while.  Hence: one conductor per request,
    a delay between requests, and a disk cache so a rerun costs nothing.
  * ``lmfdb_label`` takes the dotted label (``37.a1``); the Cremona label is
    returned as ``Clabel`` (``37a1``).
  * Real numbers arrive wrapped as
    ``{"__RealLiteral__": 0, "data": "0.3059...", "prec": 97}``.

Columns available here cover rank, analytic_rank, regulator, sha and torsion.
The real period and the Tamagawa numbers live in tables that the API does not
expose, so those two are cross-checked only indirectly -- through whether
#Sha_an comes out at the LMFDB's sha value.  That is stated in data/README.md
rather than papered over.

Usage:
    python tools/fetch_lmfdb.py --max-conductor 200
    python tools/fetch_lmfdb.py --labels 11.a2 37.a1 389.a1
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

API = 'https://www.lmfdb.org/api/ec_curvedata/'
ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / 'data' / 'cache' / 'lmfdb'
USER_AGENT = 'bsd-lab/0.1 (independent recomputation; https://github.com/Everaldtah/bsd-lab)'

#: Seconds between live requests.  The limiter is undocumented; this is chosen
#: to be slower than any plausible threshold rather than tuned against it.
DELAY = 1.5

#: Columns kept.  Everything else in the table is about isogeny classes, Galois
#: images and Faltings heights, none of which bsdlab computes.
KEEP = (
    'Clabel', 'lmfdb_label', 'ainvs', 'conductor', 'signD',
    'rank', 'analytic_rank', 'regulator', 'sha', 'torsion',
    'torsion_structure', 'num_int_pts', 'jinv', 'cm', 'bad_primes',
)


class RateLimited(RuntimeError):
    """The API returned an HTML page instead of JSON."""


def _unwrap(value):
    """Turn LMFDB's wrapped literals into plain Python.

    Real numbers come back as ``{"__RealLiteral__":0,"data":"0.30599","prec":97}``
    and are kept as the *string* ``data`` so no precision is lost to a float on
    the way in; the comparison code decides what precision to read them at.
    """
    if isinstance(value, dict):
        if '__RealLiteral__' in value:
            return value['data']
        return {k: _unwrap(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_unwrap(v) for v in value]
    return value


def _get(params: dict) -> dict:
    query = urllib.parse.urlencode({**params, '_format': 'json'})
    req = urllib.request.Request(API + '?' + query, headers={'User-Agent': USER_AGENT})
    with urllib.request.urlopen(req, timeout=60) as resp:
        body = resp.read().decode('utf-8', 'replace')
    if not body.lstrip().startswith('{'):
        raise RateLimited(f'non-JSON response for {query!r} (rate limited or bad column)')
    return json.loads(body)


def _get_with_backoff(params: dict, attempts: int = 5) -> dict:
    wait = 20.0
    for i in range(attempts):
        try:
            return _get(params)
        except (RateLimited, urllib.error.URLError, TimeoutError) as exc:
            if i == attempts - 1:
                raise
            print(f'  retry in {wait:.0f}s ({exc})', file=sys.stderr)
            time.sleep(wait)
            wait *= 2
    raise AssertionError('unreachable')


def fetch_conductor(n: int, refresh: bool = False) -> list[dict]:
    """Every curve of conductor ``n``, from cache when possible."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f'N{n}.json'
    if path.exists() and not refresh:
        return json.loads(path.read_text())

    rows: list[dict] = []
    params = {'conductor': n}
    while True:
        payload = _get_with_backoff(params)
        rows.extend(payload.get('data', []))
        nxt = payload.get('next')
        if not nxt:
            break
        # `next` is a full URL; re-derive its query so the same headers apply.
        params = {k: v[0] for k, v in
                  urllib.parse.parse_qs(urllib.parse.urlparse(nxt).query).items()}
        params.pop('_format', None)
        time.sleep(DELAY)

    kept = [{k: _unwrap(r[k]) for k in KEEP if k in r} for r in rows]
    kept.sort(key=lambda r: r.get('Clabel', ''))
    path.write_text(json.dumps(kept, indent=1, sort_keys=True))
    return kept


def fetch_label(lmfdb_label: str, refresh: bool = False) -> dict | None:
    """One curve by dotted LMFDB label, e.g. ``37.a1``."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f'L{lmfdb_label}.json'
    if path.exists() and not refresh:
        return json.loads(path.read_text())
    payload = _get_with_backoff({'lmfdb_label': lmfdb_label})
    rows = payload.get('data', [])
    if not rows:
        return None
    kept = {k: _unwrap(rows[0][k]) for k in KEEP if k in rows[0]}
    path.write_text(json.dumps(kept, indent=1, sort_keys=True))
    return kept


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--max-conductor', type=int, default=0,
                    help='fetch every curve of conductor 1..MAX')
    ap.add_argument('--min-conductor', type=int, default=1)
    ap.add_argument('--labels', nargs='*', default=[],
                    help='dotted LMFDB labels, e.g. 37.a1')
    ap.add_argument('--refresh', action='store_true', help='ignore the disk cache')
    ap.add_argument('--out', default='data/lmfdb_reference.json')
    args = ap.parse_args()

    collected: list[dict] = []

    for label in args.labels:
        row = fetch_label(label, refresh=args.refresh)
        print(f'{label}: {"ok" if row else "NOT FOUND"}')
        if row:
            collected.append(row)
        time.sleep(DELAY)

    if args.max_conductor:
        for n in range(args.min_conductor, args.max_conductor + 1):
            cached = (CACHE / f'N{n}.json').exists() and not args.refresh
            rows = fetch_conductor(n, refresh=args.refresh)
            collected.extend(rows)
            if rows:
                print(f'N={n}: {len(rows)} curves{" (cached)" if cached else ""}')
            if not cached:
                time.sleep(DELAY)

    if collected:
        out = ROOT / args.out
        out.parent.mkdir(parents=True, exist_ok=True)
        collected.sort(key=lambda r: (r.get('conductor', 0), r.get('Clabel', '')))
        out.write_text(json.dumps(collected, indent=1, sort_keys=True))
        print(f'\nwrote {len(collected)} curves to {args.out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
