"""Experiment E4/E1: the p-adic BSD quotient

    B_p(E) = [T^r] L_p(E,T) * |tors|^2 * log_p(1+p)^r
             / ( (1 - 1/alpha)^2 * Reg_p * prod c_v )

for 37a1 (r=1), 389a1, 571b1 (r=2), 5077a1 (r=3) at good ordinary p.
p-adic BSD (MTT) predicts B_p = #Sha = 1 for all four.

STATUS OF THIS EXPERIMENT (docs/research/01_global_rational_shadow.md
§8.1): it independently reproduces known numerics (Stein-Wuthrich 2013 Thm
12.3 for 389a1) with a Sage-free implementation.  It cannot test the
uniformity in p that GRS adds, and a mismatch would most likely be a bug.

Precision: c_r is taken at Riemann-sum level n; its digits are the ones
agreeing with level n-1 (heuristic stability, NOT Sage's proven bound).
Every other factor is exact to its stated precision.

Usage: python tests_research/e4_beta_table.py [pmax] [budget]
Writes data/e4_beta_table.json incrementally.
"""
import json
import os
import sys
import time
from fractions import Fraction as F
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bsdlab.curve import EllipticCurve
from bsdlab.padicl import Measure, PlusSymbols, lp_series, unit_root
from bsdlab.pheights import _plog_unit, padic_regulator
from bsdlab.reduction import conductor, tamagawa_numbers

CURVES = {
    "37a1": ([0, 0, 1, -1, 0], [(0, 0)]),
    "389a1": ([0, 1, 1, -2, 0], [(-1, 1), (0, 0)]),
    "571b1": ([0, 1, 1, -4, 2], [(1, 0), (0, 1)]),
    "5077a1": ([0, 0, 1, -7, 6], [(0, -3), (-1, -4), (1, -1)]),
}
OUT = Path(os.environ.get("E4_OUT", str(ROOT / "data" / "e4_beta_table.json")))


def vp(x, p, cap):
    if x == 0:
        return cap
    v = 0
    while x % p == 0 and v < cap:
        x //= p
        v += 1
    return v


def primes(lo, hi):
    return [q for q in range(lo, hi + 1)
            if q > 1 and all(q % d for d in range(2, int(q ** 0.5) + 1))]


def choose_level(p, budget, r):
    n = 2
    while (p - 1) * p ** n <= budget:
        n += 1
    return max(n, 2)


def beta_cell(E, lab, gens, r, p, P, budget, K=14):
    """One (curve, p) cell -> dict."""
    t0 = time.time()
    R = p ** K
    n = choose_level(p, budget, r)
    M = Measure(P, E, p, K=K)
    hi = lp_series(M, n, r)[r]
    lo = lp_series(M, n - 1, r)[r]
    agree = vp((hi - lo) % R, p, K)          # absolute precision of c_r
    vc = vp(hi % p ** agree, p, agree)
    if vc >= agree:
        return {"p": p, "level": n, "status": "c_r not resolved",
                "c_r_abs_prec": agree, "secs": round(time.time() - t0, 1)}
    u_reg, e_reg = padic_regulator(E, gens, p, K)
    v_reg_u = vp(u_reg, p, K)
    if v_reg_u >= K:
        return {"p": p, "level": n, "status": "Reg_p not resolved",
                "secs": round(time.time() - t0, 1)}
    alpha, _ = unit_root(E, p, K)
    eps = (1 - pow(alpha, -1, R)) ** 2 % R
    v_eps = vp(eps, p, K)
    lg = _plog_unit(1 + p, p, K)               # valuation 1
    v_lg = vp(lg, p, K)
    tam = 1
    for c in tamagawa_numbers(E).values():
        tam *= c
    # valuations and relative precisions
    rel = min(agree - vc, K - v_reg_u, K - v_eps, K - v_lg)
    val = vc + r * v_lg - v_eps - (v_reg_u - e_reg) - vp(tam, p, K)
    mod = p ** rel
    unit = lambda x, v: (x // p ** v) % mod
    num = unit(hi, vc) * pow(unit(lg, v_lg), r, mod)
    den = unit(eps, v_eps) * unit(u_reg, v_reg_u) * (tam // p ** vp(tam, p, K))
    b_unit = num * pow(den % mod, -1, mod) % mod
    ok = (val == 0 and b_unit % mod == 1 % mod)
    return {"p": p, "a_p": E.ap(p), "level": n, "valuation": val,
            "rel_prec": rel, "B_p_unit_mod_p^rel": b_unit,
            "c_r_val": vc, "Reg_p_val": v_reg_u - e_reg, "eps_val": v_eps,
            "B_p_equals_1": ok, "status": "ok",
            "secs": round(time.time() - t0, 1)}


def main():
    pmax = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    budget = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    only = sys.argv[3].split(",") if len(sys.argv) > 3 else list(CURVES)
    table = json.load(open(OUT)) if OUT.exists() else {}
    for lab in only:
        ai, gl = CURVES[lab]
        E = EllipticCurve(*ai)
        N = conductor(E)
        gens = [(F(a), F(b)) for a, b in gl]
        r = len(gens)
        t = time.time()
        P = PlusSymbols(E, N)
        print("%s: symbols ready (N=%d, dim %d) %.0fs" % (lab, N, P.M.dim,
                                                        time.time() - t), flush=True)
        rows = table.setdefault(lab, {})
        for p in primes(5, pmax):
            if str(p) in rows or N % p == 0 or E.ap(p) % p == 0:
                continue
            try:
                cell = beta_cell(E, lab, gens, r, p, P, budget)
            except Exception as ex:  # noqa: BLE001
                cell = {"p": p, "status": "error: %s" % ex}
            rows[str(p)] = cell
            json.dump(table, open(OUT, "w"), indent=1)
            print(lab, json.dumps(cell), flush=True)


if __name__ == "__main__":
    main()
