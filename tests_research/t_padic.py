"""A1: padic acceptance tests. Run: python scratch/t_padic.py"""
import random
import sys

sys.path.insert(0, ".")
from bsdlab.padic import Padic, psqrt, plog, pexp


def test_algebra():
    random.seed(1)
    for p in (5, 7):
        for _ in range(200):
            K = 20
            a = Padic.from_int(p, K, random.randrange(1, p ** 6))
            b = Padic.from_int(p, K, random.randrange(1, p ** 6))
            lhs = (a + b) * (a - b)
            rhs = a * a - b * b
            assert lhs == rhs, (p, a, b, lhs, rhs)
    print("A1a (a+b)(a-b) == a^2-b^2 mod p^20, p=5,7: PASS")


def test_psqrt():
    random.seed(2)
    for p in (5, 7, 13):
        for _ in range(100):
            K = 12
            n = random.randrange(1, p ** 5)
            x = Padic.from_int(p, K, n)
            try:
                r = psqrt(x)
            except ValueError as err:
                if not x.is_unit():
                    continue  # non-unit: correctly refused
                # unit refused => must be a genuine non-residue mod p
                assert "non-residue" in str(err)
                assert pow(n % p, (p - 1) // 2, p) == p - 1
                continue
            assert r * r == x, (p, x, r)
    print("A1b psqrt(x)^2 == x (non-residues correctly refused): PASS")


def test_plog_pexp():
    random.seed(3)
    for p in (5, 7):
        for _ in range(100):
            K = 10
            # x in 1 + pZ_p (convergence region of pexp)
            x = Padic.from_int(p, K, 1 + p * random.randrange(0, p ** 4))
            y = x - Padic.from_int(p, K, 1)  # v_p(y) >= 1
            e = pexp(y)
            assert plog(e) == y, (p, x, y, e)  # log(exp(y)) = y
            assert e == pexp(plog(e))
    print("A1c plog(pexp(x)) == x on 1+pZ_p, p=5,7: PASS")


if __name__ == "__main__":
    test_algebra()
    test_psqrt()
    test_plog_pexp()
