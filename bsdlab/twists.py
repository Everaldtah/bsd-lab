"""Quadratic twists E_d and the frontier-curve twist census (PLAN 5.2(a)).

For an integral model with a1 = a3 = 0, E: y^2 = x^3 + A x^2 + B x + C,
the quadratic twist by a fundamental discriminant d is
E_d: y^2 = x^3 + A d x^2 + B d^2 x + C d^3  (map (x,y) -> (x/d, y/(d sqrt d))).
Any E is Q-isomorphic to such a model: scaling x = X/4, Y = 8y + 4a1 x + 4a3
puts it at [0, b2, 0, 8 b4, 16 b6].  The twist is then fed through the lab's
own exact pipeline (conductor, root number, L-series) -- no shortcuts.
"""

from __future__ import annotations

from sympy.functions.combinatorial.numbers import jacobi_symbol

from bsdlab.curve import EllipticCurve


def kronecker(a: int, n: int) -> int:
    """The Kronecker symbol (a/n), n != 0 (n even handled by the 2-rules)."""
    if n == -1:
        return -1 if a < 0 else 1
    if n == 2:
        if a % 2 == 0:
            return 0
        return 1 if a % 8 in (1, 7) else -1
    sign = 1
    if n < 0:
        sign = kronecker(a, -1)
        n = -n
    return sign * jacobi_symbol(a % n, n)


def fundamental_discriminants(bound: int) -> list:
    """All fundamental discriminants d with 0 < |d| <= bound.

    d is fundamental iff d is squarefree with d == 1 mod 4, or d = 4m with
    m squarefree, m == 2 or 3 mod 4 (d = 4m never has m == 1 mod 4: that m
    would make d = 4m non-fundamental since m == 1 mod 4 => d/4 integral
    discriminant).  Ordered by |d|, then sign.
    """
    out = []
    for k in range(1, bound + 1):
        for d in (-k, k):
            if d % 4 == 1 and _is_squarefree(d):
                out.append(d)
            elif d % 16 in (8, 12):
                m = d // 4
                if _is_squarefree(m):
                    out.append(d)
    return out


def _is_squarefree(n: int) -> bool:
    n = abs(n)
    p = 2
    while p * p <= n:
        if n % (p * p) == 0:
            return False
        p += 1
    return n != 0


def a13_zero_model(E: EllipticCurve) -> EllipticCurve:
    """A Q-isomorphic integral model [0, b2, 0, 8*b4, 16*b6] (a1 = a3 = 0).

    From Y = 8y + 4a1 x + 4a3, X = 4x one gets Y^2 = X^3 + b2 X^2 + 8b4 X
    + 16b6 exactly; j is preserved, so this is checked against E.
    """
    M = EllipticCurve.from_list([0, E.b2, 0, 8 * E.b4, 16 * E.b6])
    if M.j_invariant != E.j_invariant:
        raise ArithmeticError("a13-zero model changed j for %r" % (E,))
    return M


def twist(E: EllipticCurve, d: int) -> EllipticCurve:
    """The quadratic twist E_d as a bsdlab EllipticCurve (nonminimal model).

    y^2 = x^3 + A x^2 + B x + C  ->  y^2 = x^3 + A d x^2 + B d^2 x + C d^3.
    """
    A, B, C = E.a2, E.a4, E.a6
    return EllipticCurve.from_list(
        [0, A * d, 0, B * d * d, C * d * d * d])
