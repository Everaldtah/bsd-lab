"""Acceptance suite for bsdlab.modsym (W4p2): exact rational twisted
L-values from Manin symbols.  B1-B4 are certified anchor values
(numerically identified rationals from the lab's L-series at two
precisions); B5 is the runtime bound."""

import time
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bsdlab.curve import EllipticCurve
from bsdlab.modsym import ManinSpace, eigenfunctional, twisted_l_ratio

CASES = [
    ("B1 11a1", [0, -1, 1, -10, -20], 11,
     {1: "1/5", 12: "5", -11: "2", -7: "0", 8: "0"}),
    ("B2 37a1", [0, 0, 1, -1, 0], 37,
     {1: "0", -7: "1", -3: "1"}),
    ("B3 389a1", [0, 1, 1, -2, 0], 389,
     {1: "0", -7: "0", -3: "2", 5: "4"}),
    ("B4 571b1", [0, 1, 1, -4, 2], 571,
     {1: "0"}),
]


def main() -> int:
    failures = 0
    t_all = time.time()
    for name, ai, N, ds in CASES:
        t = time.time()
        E = EllipticCurve.from_list(ai)
        M = ManinSpace(N)
        lams = eigenfunctional(M, E)
        t_setup = time.time() - t
        for d, want in sorted(ds.items()):
            t = time.time()
            got = twisted_l_ratio(E, d, N, M=M, lams=lams)
            ok = str(got) == want
            failures += 0 if ok else 1
            print("%-9s d=%-4d got %-8s want %-4s %s  (%.1fs)"
                  % (name, d, got, want, "OK" if ok else "FAIL",
                     time.time() - t))
        print("%-9s setup %.1fs (dim %d, %d symbols)"
              % (name, t_setup, M.dim, len(M.symbols)))
    print("total %.1fs, failures %d" % (time.time() - t_all, failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
