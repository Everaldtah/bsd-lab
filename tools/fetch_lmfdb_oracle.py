"""Fetch the LMFDB cross-validation oracle for every curve of conductor <= 400.

Unlike fetch_lmfdb.py (the small 14-curve reference set feeding the BSD
pipeline), this pulls *all* curves in a conductor range so that
tools/verify_enumeration.py can join the box-bounded enumeration in
data/curves.json against a complete LMFDB snapshot.

Strategy, established empirically (2026-09, extending the notes in
fetch_lmfdb.py):

  * The /api/ endpoint refuses every range syntax tried: ``conductor__lte``
    answers 404, a ``py:{'$lte': 400}`` operator value answers 404, and any
    unknown parameter (e.g. ``count``) is read as a field name and answers
    404.  There is no ranged batch at ANY bound, so falling back to a
    smaller bound would not help.
  * The /api/ endpoint DOES accept a filterless query with ``_sort`` and
    ``_fields``, paginated via the server's own ``next`` links (100 rows per
    page).  Walking ``ec_localdata`` in ``_sort=conductor`` order therefore
    delivers every Tamagawa row for small conductors in the first ~20 pages,
    terminating as soon as a page starts past the bound.
  * The *web* search route (``/EllipticCurve/Q/``) accepts dash ranges
    (``conductor=1-400``) that the API refuses, and its official Download
    button issues ``?download=1&query={'conductor': {'$gte': 1, '$lte':
    400}}&Submit=csv`` -- ALL matching rows in ONE response, with the
    default displayed columns already including lmfdb_label, lmfdb_iso,
    conductor, rank, torsion_structure and ainvs.  The curve list comes
    from there (single call); Tamagawa numbers come from the API walk
    (they live in ec_localdata, one row per bad prime, joined here on
    lmfdb_label).
  * The rate limiter (reCAPTCHA interstitial returning HTML where JSON or
    CSV was expected) is survived, not avoided: back off 160s then double.
    A 404 is a semantic refusal and is raised immediately instead.
  * ainvs arrives as "[0, -1, 1, -7820, -263580]" in the CSV; it is parsed
    to the list-of-ints shape data/curves.json uses.

Call budget: 1 CSV download + ~1 API call per 100 Tamagawa rows
(~20-25 pages for bound 400) + a handful of probes -- O(tens), never
thousands.  The output file doubles as its own resume checkpoint (raw
Tamagawa rows and the pending ``next`` offset live in meta until the run
completes), so an interrupted run loses at most the page in flight.

Usage:
    python tools/fetch_lmfdb_oracle.py                 # conductors 1..400
    python tools/fetch_lmfdb_oracle.py --max-conductor 100
"""

from __future__ import annotations

import argparse
import csv
import http.client
import io
import json
import math
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fetch_lmfdb  # noqa: E402  (reuse _unwrap / RateLimited / USER_AGENT)

EC_SEARCH = 'https://www.lmfdb.org/EllipticCurve/Q/'
LOCAL_API = 'https://www.lmfdb.org/api/ec_localdata/'
ROOT = Path(__file__).resolve().parent.parent

#: Politeness delay between HTTP calls (brief asks for ~0.5s; the extra
#: margin is for the undocumented rate limiter).
DELAY = 2.0

#: Columns requested from the API walk.  ``conductor`` is only needed as
#: the walk's stop signal.
LOCAL_FIELDS = 'lmfdb_label,prime,tamagawa_number,conductor'

#: The walk needs ~2200 Tamagawa rows / 100 per page; the cap exists so a
#: server-side sort change cannot turn the walk into an endless loop.
MAX_PAGES = 200


class Fetcher:
    """Same request shape as fetch_lmfdb._get, parameterised by URL."""

    def __init__(self, calls: int = 0):
        self.calls = calls

    def get_text(self, url: str) -> str:
        req = urllib.request.Request(url, headers={'User-Agent': fetch_lmfdb.USER_AGENT})
        with urllib.request.urlopen(req, timeout=120) as resp:
            body = resp.read().decode('utf-8', 'replace')
        self.calls += 1
        return body

    def get_with_backoff(self, url: str, attempts: int = 6, expect: str = 'json') -> str:
        # A rejected request seems to extend the limiter's penalty window,
        # so the first retry must already be longer than a full cooldown
        # (~100-150s, observed) rather than a short network-style backoff.
        # The format check lives INSIDE the loop: the limiter serves HTML
        # where JSON/CSV was expected, and that must retry like any other
        # rejection rather than escape and kill the run.
        wait = 160.0
        for i in range(attempts):
            try:
                body = self.get_text(url)
                first = body.splitlines()[0] if body else ''
                if expect == 'json' and not body.lstrip().startswith('{'):
                    raise fetch_lmfdb.RateLimited('non-JSON response')
                if expect == 'csv' and 'HYPERLINK' not in first:
                    raise fetch_lmfdb.RateLimited('download returned HTML, not CSV')
                return body
            except urllib.error.HTTPError:
                # A 404 is a semantic refusal, not congestion; re-sending
                # the same query would only hammer the limiter.
                raise
            except (fetch_lmfdb.RateLimited, urllib.error.URLError, TimeoutError,
                    ConnectionError, http.client.HTTPException) as exc:
                if i == attempts - 1:
                    raise
                print(f'  retry in {wait:.0f}s ({exc})', file=sys.stderr, flush=True)
                time.sleep(wait)
                wait *= 2
        raise AssertionError('unreachable')


def download_curves(fetcher: Fetcher, lo: int, hi: int) -> list[dict]:
    """Every curve with ``lo <= conductor <= hi``, via the web Download button."""
    query = json.dumps({'conductor': {'$gte': lo, '$lte': hi}},
                       separators=(', ', ': ')).replace('"', "'")
    params = urllib.parse.urlencode({'download': '1', 'query': query,
                                     'conductor': f'{lo}-{hi}',
                                     'count': '50', 'Submit': 'csv'})
    body = fetcher.get_with_backoff(EC_SEARCH + '?' + params, expect='csv')
    rows = list(csv.reader(io.StringIO(body)))
    header, data = rows[0], rows[1:]
    # Header cells look like =HYPERLINK(".../ec.q.conductor", "conductor").
    names = [re.search(r', "([^"]+)"\)$', c).group(1) if re.search(r', "([^"]+)"\)$', c) else c
             for c in header]
    idx = {name: i for i, name in enumerate(names)}
    for need in ('lmfdb_label', 'conductor', 'rank', 'torsion_structure', 'ainvs'):
        if need not in idx:
            raise RuntimeError(f'CSV download lacks the {need!r} column; '
                               f'columns were {names}')

    out: list[dict] = []
    for r in data:
        tors = json.loads(r[idx['torsion_structure']]) if r[idx['torsion_structure']].strip() else []
        out.append({
            'lmfdb_label': r[idx['lmfdb_label']],
            'lmfdb_iso': r[idx['lmfdb_iso']] if 'lmfdb_iso' in idx else None,
            'conductor': int(r[idx['conductor']]),
            'rank': int(r[idx['rank']]),
            'torsion': math.prod(tors) if tors else 1,
            'ainvs': json.loads(r[idx['ainvs']]),
        })
    return out


def walk_localdata(fetcher: Fetcher, bound: int, checkpoint: dict | None) -> list[dict]:
    """All ec_localdata rows with conductor <= bound, via the API sort-walk.

    ``checkpoint`` (from a previous interrupted run) carries already-fetched
    rows plus the ``next`` URL to resume from.
    """
    rows: list[dict] = list(checkpoint['rows']) if checkpoint else []
    url: str | None = checkpoint['next_url'] if checkpoint else (
        LOCAL_API + '?' + urllib.parse.urlencode(
            {'_format': 'json', '_sort': 'conductor', '_fields': LOCAL_FIELDS}))

    for page in range(checkpoint['pages'] if checkpoint else 0, MAX_PAGES):
        if url is None:
            return rows
        body = fetcher.get_with_backoff(url, expect='json')
        payload = json.loads(body)
        data = payload.get('data', [])
        if not data:
            return rows  # empty page = end (its `next` is a bogus _offset=0 self-link)
        if int(data[0]['conductor']) > bound:
            return rows  # _sort=conductor: nothing later can be inside the bound
        for r in data:
            if int(r['conductor']) <= bound:
                rows.append({k: fetch_lmfdb._unwrap(r[k]) for k in
                             ('lmfdb_label', 'prime', 'tamagawa_number')})
        nxt = payload.get('next')
        next_offset = int(urllib.parse.parse_qs(
            urllib.parse.urlparse(nxt).query).get('_offset', ['0'])[0]) if nxt else 0
        done = (not nxt) or next_offset <= 0
        # Only a validated next URL goes into the checkpoint: an empty page
        # carries a bogus self-referential next (_offset=0) that serves HTML
        # if ever followed after a crash-resume.
        yield_checkpoint = {
            'rows': rows,
            'next_url': None if done else urllib.parse.urljoin('https://www.lmfdb.org', nxt),
            'pages': page + 1,
        }
        yield yield_checkpoint
        if done:
            return rows
        url = yield_checkpoint['next_url']
        time.sleep(DELAY)
    raise RuntimeError(f'page cap of {MAX_PAGES} exceeded during ec_localdata walk')


def build_oracle(curves: list[dict], local_rows: list[dict]) -> tuple[list[dict], list[str]]:
    """Attach per-prime Tamagawa numbers to the curve rows."""
    tam: dict[str, dict[str, int]] = {}
    warnings: list[str] = []
    for row in local_rows:
        label, p = row['lmfdb_label'], row['prime']
        if str(p) in tam.get(label, {}):
            warnings.append(f'duplicate localdata row for {label} at p={p}')
        tam.setdefault(label, {})[str(p)] = row['tamagawa_number']

    out: list[dict] = []
    for c in curves:
        label = c['lmfdb_label']
        n = c['conductor']
        t = tam.pop(label, None)
        if t is None:
            warnings.append(f'no ec_localdata rows for {label}; Tamagawa left null')
        for p in (t or {}):
            if n % int(p):
                warnings.append(f'localdata prime {p} does not divide conductor {n} ({label})')
        out.append({
            'ainvs': c['ainvs'],
            'conductor': n,
            'rank': c['rank'],
            'torsion': c['torsion'],
            'tamagawa': t,
            'tamagawa_product': math.prod(t.values()) if t else None,
            'lmfdb_label': label,
            'lmfdb_iso': c['lmfdb_iso'],
        })
    for label in sorted(tam):
        warnings.append(f'ec_localdata label {label} has no ec_curvedata row')
    return out, warnings


def write_oracle(path: Path, blob: dict) -> None:
    """Atomic rewrite so a crash mid-write cannot corrupt the checkpoint."""
    tmp = path.parent / (path.name + '.tmp')
    tmp.write_text(json.dumps(blob, indent=1, sort_keys=True))
    os.replace(tmp, path)


API_NOTES = (
    'The LMFDB /api/ endpoint refuses all range syntax (conductor__lte and '
    'py:{\'$lte\': ...} operator values both answer 404, as does any unknown '
    'parameter such as count), so no ranged batch exists at any bound.  The '
    'curve list therefore comes from ONE call to the web search route\'s '
    'official Download button (?download=1&query={\'conductor\': {\'$gte\': 1, '
    '\'$lte\': 400}}&Submit=csv), which accepts the dash range the API '
    'rejects and returns every row at once with ainvs/conductor/rank/ '
    'torsion.  Tamagawa numbers are not in that table: they were walked '
    'from /api/ec_localdata with a filterless _sort=conductor query, '
    'following the server\'s next links (~100 rows/page) and stopping at '
    'the first page whose rows all exceed the bound, joined to the curves '
    'on lmfdb_label.'
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--max-conductor', type=int, default=400)
    ap.add_argument('--min-conductor', type=int, default=1)
    ap.add_argument('--out', default='data/lmfdb_oracle.json')
    ap.add_argument('--no-resume', action='store_true',
                    help='start over even if a partial oracle exists')
    args = ap.parse_args()

    out_path = ROOT / args.out
    started = datetime.now(timezone.utc).isoformat()
    blob = None
    if out_path.exists() and not args.no_resume:
        blob = json.loads(out_path.read_text())
        meta = blob.get('meta', {})
        if meta.get('bound') != args.max_conductor:
            print('existing oracle has a different bound; starting over', flush=True)
            blob = None
        elif 'curves' in blob and blob.get('curves'):
            print(f'oracle already complete ({len(blob["curves"])} curves); nothing to do',
                  flush=True)
            return 0

    fetcher = Fetcher(calls=(blob['meta']['http_calls'] if blob else 0))

    # Phase 1: the curve list (1 call, re-issued on resume -- it is cheap
    # and idempotent, unlike caching a second raw copy inside the oracle).
    print(f'downloading curves with {args.min_conductor} <= N <= {args.max_conductor} ...',
          flush=True)
    curves = download_curves(fetcher, args.min_conductor, args.max_conductor)
    if curves:
        print(f'  {len(curves)} curves '
              f'(N = {curves[0]["conductor"]}..{curves[-1]["conductor"]})', flush=True)
    else:
        print('  0 curves (no elliptic curve has a conductor in that range)', flush=True)

    # Phase 2: the Tamagawa walk, checkpointed into the oracle file after
    # every page so an interruption loses at most one page.
    checkpoint = blob['meta'].get('localdata') if blob else None
    if checkpoint:
        print(f'resuming localdata walk: {len(checkpoint["rows"])} rows, '
              f'page {checkpoint["pages"]}', flush=True)
    state = {'rows': [], 'pages': 0}  # stays empty if the walk yields nothing
    for state in walk_localdata(fetcher, args.max_conductor, checkpoint):
        write_oracle(out_path, {'meta': {
            'bound': args.max_conductor,
            'started': started,
            'updated': datetime.now(timezone.utc).isoformat(),
            'http_calls': fetcher.calls,
            'localdata': state,
            'api_notes': API_NOTES,
        }})
        print(f'  localdata page {state["pages"]}: {len(state["rows"])} rows so far',
              flush=True)
    local_rows = state['rows']

    # Phase 3: join and final write (raw walk state dropped from meta).
    oracle, warnings = build_oracle(curves, local_rows)
    for w in warnings:
        print(f'  WARNING: {w}', file=sys.stderr, flush=True)
    write_oracle(out_path, {
        'meta': {
            'bound': args.max_conductor,
            'started': started,
            'updated': datetime.now(timezone.utc).isoformat(),
            'http_calls': fetcher.calls,
            'api_notes': API_NOTES,
        },
        'curves': sorted(oracle, key=lambda r: (r['conductor'], r['lmfdb_label'])),
    })
    print(f'\nwrote {len(oracle)} curves '
          f'(conductors {args.min_conductor}..{args.max_conductor}) to {args.out}; '
          f'{fetcher.calls} HTTP calls total', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
