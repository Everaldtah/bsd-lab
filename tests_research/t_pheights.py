"""Gates for bsdlab.pheights (W4p1b). External anchors are Sage doctests
(sage/schemes/elliptic_curves/padics.py), whose normalisation makes MTT
p-adic BSD hold.

H1  h_5(37a1, (0,0)) = 5 + 5^2 + 5^3 + 3*5^6 + 4*5^7 + 5^9 + O(5^10).
H2  Reg_5(5077a1) = 5 + 5^2 + 4*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 4*5^7 + 2*5^8
    + 5^9 + O(5^10)   (rank 3, anomalous: 5 | #E(F_5) = 10).
    The basis must be a Z-basis: (0,-3), (-1,-4), (1,-1) has classical
    regulator 0.417143558758384 = the tabulated value; the often-quoted
    (-3,0), (-2,3), (-1,3) spans an index-2 sublattice (regulator 4x).
H3  quadratic form: h(2P) = 4 h(P), h(3P) = 9 h(P) for 37a1 at p = 5, 7, 11.
H4  torsion: h = 0 on the 5-torsion point (5,5) of 11a1 at p = 7.
"""
import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bsdlab.curve import EllipticCurve
from bsdlab.pheights import PadicHeight, padic_regulator


def value_digits(u, e, p, k):
    """base-p digits of u / p^e starting at p^0 (requires v_p(u) >= e)."""
    if u % p ** e:
        return None
    x = u // p ** e
    out = []
    for _ in range(k):
        out.append(x % p)
        x //= p
    return out


def main() -> int:
    fails = 0
    E37 = EllipticCurve(0, 0, 1, -1, 0)
    H = PadicHeight(E37, 5, 12)
    got = value_digits(*H((F(0), F(0))), 5, 10)
    want = [0, 1, 1, 1, 0, 0, 3, 4, 0, 1]
    fails += got != want
    print("H1 h_5(37a1 gen) %s vs Sage %s %s" % (got, want,
                                                  "PASS" if got == want else "FAIL"))

    E5077 = EllipticCurve(0, 0, 1, -7, 6)
    basis = [(F(0), F(-3)), (F(-1), F(-4)), (F(1), F(-1))]
    got = value_digits(*padic_regulator(E5077, basis, 5, 16), 5, 10)
    want = [0, 1, 1, 4, 2, 2, 2, 4, 2, 1]
    fails += got != want
    print("H2 Reg_5(5077a1) %s vs Sage %s %s" % (got, want,
                                                 "PASS" if got == want else "FAIL"))

    P = (F(0), F(0))
    for p in (5, 7, 11):
        H = PadicHeight(E37, p, 10)
        R = p ** 8
        h1 = H(P)[0]
        h2 = H(E37.multiply(P, 2))[0]
        h3 = H(E37.multiply(P, 3))[0]
        ok = (h2 - 4 * h1) % R == 0 and (h3 - 9 * h1) % R == 0 and h1 % R
        fails += not ok
        print("H3 37a1 p=%d h(2P)=4h(P), h(3P)=9h(P) mod p^8: %s"
              % (p, "PASS" if ok else "FAIL"))

    E11 = EllipticCurve(0, -1, 1, -10, -20)
    u, e = PadicHeight(E11, 7, 8)((F(5), F(5)))
    ok = u % 7 ** 8 == 0
    fails += not ok
    print("H4 11a1 torsion (5,5) h_7 = 0: %s" % ("PASS" if ok else "FAIL"))
    print("failures", fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
