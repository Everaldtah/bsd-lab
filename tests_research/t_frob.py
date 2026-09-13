"""Gates for bsdlab.frobenius (W4p1a, orchestrator implementation).

G1  charpoly: trace F == a_p, det F == p (mod p^K) -- proves the matrix.
G2  truncation stability: F at N and N+3 agree mod p^K.
G3  Serre/Katz congruence, evaluated on (E, omega):
        E_{p+1}(E, omega) == E2(E, omega) * E_{p-1}(E, omega)   (mod p)
    with E4 = c4, E6 = -c6.  (As q-expansions E_{p+1} == E2 mod p; as
    functions on curves the Hasse invariant A = E_{p-1} appears, since
    E_{p-1} == 1 only on q-expansions.  The first version of this gate
    omitted A and failed 20/29 correct values.)
G4  External anchor: Sage/MAGMA doctest (sage/schemes/elliptic_curves/
    padics.py, padic_E2): E2(37a1 = [-1,1/4], p=5) =
    2 + 4*5 + 2*5^3 + 5^4 + 3*5^5 + 2*5^6 + 5^8 + 3*5^9 + 4*5^10 + 2*5^11
    + 2*5^12 + 2*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 4*5^18 + 2*5^19 + O(5^20).
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bsdlab.curve import EllipticCurve
from bsdlab.frobenius import (e2_padic, eisenstein_value_mod_p,
                              frobenius_matrix)

CURVES = {"11a1": [0, -1, 1, -10, -20], "37a1": [0, 0, 1, -1, 0],
          "389a1": [0, 1, 1, -2, 0], "571b1": [0, 1, 1, -4, 2],
          "5077a1": [0, 0, 1, -7, 6], "32a3": [0, 0, 0, 4, 0]}
PRIMES = [5, 7, 11, 13, 17, 19]
K = 8


def main() -> int:
    fails = 0
    for lab, ai in CURVES.items():
        E = EllipticCurve.from_list(ai)
        for p in PRIMES:
            if E.discriminant % p == 0:
                continue
            ap = E.ap(p)
            t = time.time()
            F = frobenius_matrix(E, p, K)
            R = p ** K
            tr = (F[0][0] + F[1][1]) % R
            det = (F[0][0] * F[1][1] - F[0][1] * F[1][0]) % R
            g1 = tr == ap % R and det == p % R
            F2 = frobenius_matrix(E, p, K, N=K + 6)
            g2 = F2 == F
            line = "%-6s p=%-2d a_p=%-3d G1 %s G2 %s" % (
                lab, p, ap, "PASS" if g1 else "FAIL", "PASS" if g2 else "FAIL")
            fails += (not g1) + (not g2)
            if ap % p:
                e2 = e2_padic(E, p, K - 2)
                ek = eisenstein_value_mod_p(E, p + 1, p)
                hasse = eisenstein_value_mod_p(E, p - 1, p)
                g3 = (e2 * hasse - ek) % p == 0
                fails += not g3
                line += " G3 %s (E2 mod p=%d, E_%d=%d) E2 mod p^%d=%d" % (
                    "PASS" if g3 else "FAIL", e2 % p, p + 1, ek, K - 2, e2)
            else:
                line += " (supersingular: no E2)"
            print(line + "  %.1fs" % (time.time() - t), flush=True)
    digits = {0: 2, 1: 4, 3: 2, 4: 1, 5: 3, 6: 2, 8: 1, 9: 3, 10: 4, 11: 2,
              12: 2, 14: 2, 15: 3, 16: 3, 17: 3, 18: 4, 19: 2}
    want = sum(d * 5 ** e for e, d in digits.items())
    got = e2_padic(EllipticCurve.from_list(CURVES["37a1"]), 5, 20)
    g4 = got == want
    fails += not g4
    print("G4 E2(37a1, 5) mod 5^20 vs Sage/MAGMA: %s (got %d, want %d)"
          % ("PASS" if g4 else "FAIL", got, want))
    print("failures", fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
