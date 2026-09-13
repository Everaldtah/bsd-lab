"""Gates for bsdlab.padicl (W4p3: MTT p-adic L-series).

P1  distribution relation (theorem; = Hecke T_p + alpha^2 - a_p alpha + p = 0),
    37a1 and 389a1 at p = 5, 7, levels n = 1..3.
P2  external anchor, Sage doctest (padic_lseries.measure, 37a, p=5):
    Sage's measure is in w1 units, ours in Neron units Omega_E = 2 w1, so
    Sage = 2 * ours:  mu(1 + 25 Z_5) = 2 + 3*5 + 4*5^3 + 2*5^4 + 3*5^5 +
    3*5^6 + 4*5^7 + 4*5^8 + O(5^9).
P3  external anchor, Sage doctest (padic_lseries class docstring, 389a,
    p=5): series(2) = O(5^4) + O(5)T + 4T^2 + 2T^3 + 3T^4 (mod 5 each);
    series(3) = O(5^5) + O(5^2)T + (4+4*5)T^2 + (2+4*5)T^3 + (3+O(5^2))T^4.
P4  interpolation (theorem, MTT): L_p(0) = (1 - 1/alpha)^2 L(E,1)/Omega_E
    for 11a1 [L/Omega = 1/5] at p = 7, 13; checked at level 5 mod p^4.
P5  p-adic BSD order of vanishing (prediction, reported): 37a1 at p=5,7
    has c0 = 0 and c1 != 0; 389a1 at p=5,7 has c0 = c1 = 0, c2 != 0.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bsdlab.curve import EllipticCurve
from bsdlab.padicl import Measure, PlusSymbols, lp_series, unit_root


def digits(x, p, k):
    out = []
    for _ in range(k):
        out.append(x % p)
        x //= p
    return out


def vp(x, p, cap):
    v = 0
    while x % p == 0 and v < cap:
        x //= p
        v += 1
    return v


def main() -> int:
    fails = 0
    E37 = EllipticCurve(0, 0, 1, -1, 0)
    E389 = EllipticCurve(0, 1, 1, -2, 0)
    E11 = EllipticCurve(0, -1, 1, -10, -20)
    syms = {"37a1": PlusSymbols(E37, 37), "389a1": PlusSymbols(E389, 389),
            "11a1": PlusSymbols(E11, 11)}
    curves = {"37a1": E37, "389a1": E389, "11a1": E11}

    for lab in ("37a1", "389a1"):
        for p in (5, 7):
            M = Measure(syms[lab], curves[lab], p, K=10)
            R = p ** 10
            bad = sum(1 for n in (1, 2, 3) for a in range(1, p ** n)
                      if a % p and M.mu(a, n) !=
                      sum(M.mu(a + j * p ** n, n + 1) for j in range(p)) % R)
            fails += bad > 0
            print("P1 %-5s p=%d distribution failures %d %s"
                  % (lab, p, bad, "PASS" if not bad else "FAIL"))

    M = Measure(syms["37a1"], E37, 5, K=9)
    sage = [2, 3, 0, 4, 2, 3, 3, 4, 4]
    got = digits(2 * M.mu(1, 2) % 5 ** 9, 5, 9)
    fails += got != sage
    print("P2 37a1 2*mu(1+25Z_5) %s vs Sage %s %s"
          % (got, sage, "PASS" if got == sage else "FAIL"))

    M = Measure(syms["389a1"], E389, 5, K=12)
    c2 = lp_series(M, 2, 4)
    c3 = lp_series(M, 3, 4)
    ok2 = (vp(c2[0], 5, 4) >= 4 and c2[1] % 5 == 0 and
           [c2[i] % 5 for i in (2, 3, 4)] == [4, 2, 3])
    ok3 = (vp(c3[0], 5, 5) >= 5 and c3[1] % 25 == 0 and
           [c3[i] % 25 for i in (2, 3, 4)] == [24, 22, 3])
    fails += (not ok2) + (not ok3)
    print("P3 389a1 p=5 series(2) %s, series(3) %s"
          % ("PASS" if ok2 else "FAIL", "PASS" if ok3 else "FAIL"))

    for p in (7, 13):
        K = 10
        R = p ** K
        M = Measure(syms["11a1"], E11, p, K=K)
        c0 = lp_series(M, 5 if p == 7 else 4, 1)[0]
        a, _ = unit_root(E11, p, K)
        e = (1 - pow(a, -1, R)) ** 2 % R
        want = e * pow(5, -1, R) % R
        ok = (c0 - want) % p ** 4 == 0
        fails += not ok
        print("P4 11a1 p=%d L_p(0) == (1-1/alpha)^2/5 mod p^4: %s"
              % (p, "PASS" if ok else "FAIL"))

    for lab, r in (("37a1", 1), ("389a1", 2)):
        for p in (5, 7):
            M = Measure(syms[lab], curves[lab], p, K=10)
            cs = lp_series(M, 5 if p == 5 else 4, 3)
            z = [vp(c, p, 10) for c in cs]
            ok = all(zz >= 3 for zz in z[:r]) and z[r] < 3
            print("P5 %-5s p=%d valuations c0..c3 %s -> ord_T = %d predicted: %s"
                  % (lab, p, z, r, "consistent" if ok else "INCONSISTENT"))
            fails += not ok
    print("failures", fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
