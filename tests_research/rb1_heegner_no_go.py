"""Experiment RB1: the first-derivative no-go, measured.

Route B candidate (a) (PLAN.md 5.2(a)) and every Heegner-type point source
share one mechanism: a point y_K in E(K), K imaginary quadratic, whose
canonical height is governed by the FIRST derivative (Gross-Zagier)

    h(y_K)  ~  L'(E/K, 1),        L(E/K, s) = L(E, s) L(E^K, s).

If ord_{s=1} L(E,s) >= 2 then L'(E/K,1) = 0 for every K, so y_K must be
torsion for every K -- a Heegner point can never see rank >= 2.  This
script checks that prediction numerically, independently of Sage/PARI:

    z_K = sum over Cl(K) of phi(tau_A),   phi(tau) = sum_n a_n/n q^n,

reduced modulo the Neron lattice (Manin constant 1 for these optimal
curves).  z_K in (1/m) Lambda for small m  <=>  y_K torsion.

Calibration: 37a1 (rank 1) must give z_K = k * z(P) + torsion with k != 0.
Targets: 389a1 (rank 2, w=+1) and 5077a1 (rank 3, w=-1) must give torsion.

Usage: python tests_research/rb1_heegner_no_go.py [dps] [Dmax]
Writes data/rb1_heegner.json.
"""
import json
import sys
import time
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path

import mpmath
from mpmath import mp, mpf, mpc, exp, pi, sqrt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bsdlab.curve import EllipticCurve
from bsdlab.lseries import an_coefficients
from bsdlab.periods import period_lattice, elliptic_logarithm

CURVES = {
    "37a1": ([0, 0, 1, -1, 0], 37, 1, [(0, 0)]),
    "389a1": ([0, 1, 1, -2, 0], 389, 2, [(-1, 1), (0, 0)]),
    "5077a1": ([0, 0, 1, -7, 6], 5077, 3, [(0, -3), (-1, -4), (1, -1)]),
}


def fundamental(D):
    if D % 4 == 1:
        m = -D
        return all(m % (p * p) for p in range(3, isqrt(m) + 1, 2)) and m % 4 != 0
    if D % 4 == 0:
        d = D // 4
        if d % 4 not in (2, 3):
            return False
        m = -d
        return all(m % (p * p) for p in range(2, isqrt(m) + 1))
    return False


def reduce_form(a, b, c):
    while True:
        if c < a or (c == a and b < 0):
            a, b, c = c, -b, a
            continue
        if b > a or b <= -a:
            k = (a - b) // (2 * a)
            b, c = b + 2 * k * a, a * k * k + b * k + c
            continue
        return (a, b, c)


def class_number(D):
    h = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a) == 0:
                c = (b * b - D) // (4 * a)
                if c >= a and gcd(gcd(a, abs(b)), c) == 1 and not (b < 0 and (a == c)):
                    h += 1
        a += 1
    return h


def heegner_taus(D, N):
    """One tau per class of Cl(D): forms [N a', B, C], B == beta mod 2N,
    chosen with a' small (maximises Im tau).  None if D fails Heegner hyp."""
    betas = [b for b in range(2 * N) if (b * b - D) % (4 * N) == 0]
    if not betas or gcd(D, N) != 1:          # N must split, not ramify, in K
        return None
    beta = betas[0]
    h = class_number(D)
    seen, taus = set(), []
    ap = 1
    while len(taus) < h and ap < 400:
        A = N * ap
        for B in range(-A + (beta + A) % (2 * N), A + 1, 2 * N):
            if (B * B - D) % (4 * A):
                continue
            C = (B * B - D) // (4 * A)
            if gcd(gcd(A, abs(B)), C) != 1:
                continue
            cls = reduce_form(A, B, C)
            if cls in seen:
                continue
            seen.add(cls)
            taus.append((A, B, C))
        ap += 1
    return taus if len(taus) == h else None


def phi(an, tau, nmax):
    q = exp(2j * pi * tau)
    s, qn = mpc(0), mpc(1)
    for n in range(1, nmax + 1):
        qn *= q
        if an[n]:
            s += mpf(an[n]) / n * qn
    return s


def best_denominator(xy, mmax=60, tol=None):
    """Smallest m <= mmax with m*x, m*y within tol of integers."""
    for m in range(1, mmax + 1):
        err = max(abs(m * t - mpmath.nint(m * t)) for t in xy)
        if err < tol:
            return m, err
    return None, min(abs(t - mpmath.nint(t)) for t in xy)


def run(label, dps, Dmax):
    ai, N, r, gens = CURVES[label]
    E = EllipticCurve(*ai)
    mp.dps = dps + 10
    w1, w2 = period_lattice(E, dps + 10)
    W2 = w2.imag                         # rectangular: disc > 0 for all three
    rows = []
    tol = mpf(10) ** (-(dps // 2))
    Ds = [D for D in range(-3, -Dmax - 1, -1) if fundamental(D) and D not in (-3, -4)]
    jobs = []
    for D in Ds:
        taus = heegner_taus(D, N)
        if taus:
            jobs.append((D, taus))
    Imin = min(sqrt(-D) / (2 * A) for D, t in jobs for (A, B, C) in t)
    nmax = int(dps * 2.31 / (2 * float(pi) * float(Imin))) + 20
    t0 = time.time()
    an = an_coefficients(E, nmax, N)
    print("%s: %d discriminants, a_n to %d (%.0fs)" % (label, len(jobs), nmax,
                                                       time.time() - t0), flush=True)
    gz = [elliptic_logarithm(E, (F(x), F(y)), dps + 10) for x, y in gens]
    for D, taus in jobs:
        z = mpc(0)
        for A, B, C in taus:
            tau = mpc(-B, sqrt(-D)) / (2 * A)
            nm = int(dps * 2.31 / (2 * float(pi) * float(tau.imag))) + 20
            z += phi(an, tau, nm)
        x, y = z.real / w1, z.imag / W2
        m, err = best_denominator([x, y], tol=tol)
        row = {"D": D, "h": len(taus), "z_re": mpmath.nstr(z.real, 20),
               "z_im": mpmath.nstr(z.imag, 20),
               "torsion_denominator": m, "residual": mpmath.nstr(err, 3)}
        if m is None and r == 1:
            # rank-1 calibration: z = k z(P) + (1/2-lattice torsion) ?
            zP = gz[0]
            for k in range(-40, 41):
                if k == 0:
                    continue
                d = z - k * zP
                mm, e2 = best_denominator([d.real / w1, d.imag / W2], mmax=4, tol=tol)
                if mm:
                    row["k_times_generator"] = k
                    row["k_residual"] = mpmath.nstr(e2, 3)
                    break
        rows.append(row)
        print(label, json.dumps(row), flush=True)
    return {"N": N, "rank": r, "dps": dps, "nmax": nmax, "rows": rows,
            "all_torsion": all(rw["torsion_denominator"] for rw in rows)}


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    Dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    only = sys.argv[3].split(",") if len(sys.argv) > 3 else list(CURVES)
    out = ROOT / "data" / "rb1_heegner.json"
    res = json.load(open(out)) if out.exists() else {}
    for lab in only:
        res[lab] = run(lab, dps, Dmax)
        json.dump(res, open(out, "w"), indent=1)
    if "37a1" in only:
        # Gross-Zagier cross-check on the control: exact L(E^D,1) from modular
        # symbols; GZ predicts L_alg = 0 <=> y_K torsion, and L_alg ~ k^2.
        from bsdlab.modsym import algebraic_twisted_l
        E = EllipticCurve(*CURVES["37a1"][0])
        chk = {}
        for rw in res["37a1"]["rows"]:
            v = algebraic_twisted_l(E, rw["D"], 37)
            tors = rw["torsion_denominator"] is not None
            chk[str(rw["D"])] = {"L_alg": str(v), "k": rw.get("k_times_generator"),
                                 "heegner_torsion": tors,
                                 "consistent": (v == 0) == tors}
        res["37a1_twist_L_values"] = chk
        json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: v["all_torsion"] for k, v in res.items() if "rows" in v}))


if __name__ == "__main__":
    main()
