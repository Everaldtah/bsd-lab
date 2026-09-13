"""Out-of-sample check of bsdlab.modsym: exact rational vs 40-digit numeric
L(E_d,1)/Omega(E_d) on twists NOT in the worker's anchor set."""
import sys, time
from fractions import Fraction
from mpmath import mp, nstr
from bsdlab.curve import EllipticCurve
from bsdlab import twists, modsym
from bsdlab.reduction import conductor
from bsdlab.lseries import l_value
from bsdlab.periods import real_period

mp.dps = 40
CURVES = {"11a1": [0,-1,1,-10,-20], "37a1": [0,0,1,-1,0], "389a1": [0,1,1,-2,0]}
ANCHOR = {"11a1": set(), "37a1": set(), "389a1": set()}  # post-fix: check anchors too
D = [d for d in twists.fundamental_discriminants(60) if abs(d) <= 60]
bad = 0; n = 0
for lab, ai in CURVES.items():
    E = EllipticCurve(*ai); N = conductor(E)
    M = modsym.ManinSpace(N) if 'ManinSpace' in dir(modsym) else None
    lams = modsym.eigenfunctional(M, E)
    for d in D:
        if d in ANCHOR[lab] or N % abs(d) == 0 or __import__('math').gcd(d, N) != 1:
            continue
        Ed = twists.twist(E, d).minimal_model(); Nd = conductor(Ed)
        if Nd > 60000:
            continue
        exact = modsym.twisted_l_ratio(E, d, N, M=M, lams=lams)
        num = l_value(Ed, Nd, prec=40) / real_period(Ed, prec=40)
        ok = abs(num - mp.mpf(exact.numerator)/exact.denominator) < mp.mpf(10)**-20
        n += 1; bad += (not ok)
        print(f"{lab:6s} d={d:4d} Nd={Nd:6d} modsym={str(exact):8s} numeric={nstr(num,12):>16s} {'OK' if ok else 'MISMATCH'}", flush=True)
print(f"checked {n}, mismatches {bad}")
