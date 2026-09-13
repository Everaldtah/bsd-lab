"""Local root numbers w_p at primes of bad reduction, from Kodaira data.

Table-route computation of the local root numbers, replacing the
NotImplementedError fallback in ``lseries._local_root_number`` for additive
primes.  The tables implemented here are:

- multiplicative (split / nonsplit I_n):  w_p = -1 / +1  (classical);
- additive p >= 5:  Halberstadt's classification, C. R. Acad. Sci. Paris
  326 (1998), 1047-1052 (via Rohrlich, Compositio 87 (1993), Prop. 2 and
  Rizzo, Compositio 136 (2003), Table I); expressible by Kodaira type;
- additive p = 3:  Kraus, Manuscripta Math. 69 (1990), Table II as
  reproduced by Rizzo (2003), Table II -- implemented for the decoded
  rows only;
- additive p = 2:  honest gap.  Only the derivable cases are implemented
  (good reduction; multiplicative split via c6 mod 8, Rizzo Lemma 7;
  additive potentially multiplicative via Connell's formula / Rizzo
  Lemma 9).  Every other additive p = 2 case returns None with reason
  ``P2_TABLE_MISSING`` -- never a guess.

The reduction data (Kodaira type, reduction type) comes from the lab's own
Tate's algorithm, ``bsdlab.reduction.local_data``, so this module adds no
new reduction logic and no external CAS dependency.
"""

P2_TABLE_MISSING = "additive p=2 table entry missing"

# Status of a local computation: ('ok', w) or ('none', reason).
from typing import Optional, Tuple

Status = Tuple[str, object]


def _vp(x: int, p: int) -> int:
    """v_p(x) for a nonzero integer x (v_p(0) is never needed: callers use
    float('inf') semantics via the 0 special case)."""
    if x == 0:
        return 1 << 30  # stands for +infinity; valuation arithmetic below
        # only ever compares and adds nonneg integers far below 2^30
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def _unit(x: int, p: int) -> int:
    """Unit part x / p^{v_p(x)} as an integer (0 stays 0)."""
    if x == 0:
        return 0
    while x % p == 0:
        x //= p
    return x


def _legendre(a: int, p: int) -> int:
    """Legendre symbol (a/p), p an odd prime, a coprime to p -> +1/-1."""
    a %= p
    if a == 0:
        raise ArithmeticError("Legendre symbol with p | a")
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def local_root_number(E, p) -> Status:
    """Local root number w_p of E at p, as ('ok', +1/-1) or ('none', reason).

    E is an bsdlab.curve.EllipticCurve; p a prime.  The computation is
    exact integer arithmetic on the minimal-model invariants c4, c6 and
    the Kodaira type from bsdlab.reduction.  The multiplicative case is
    cross-checked against the split/nonsplit verdict of Tate's algorithm;
    a disagreement raises ArithmeticError rather than returning a value.
    """
    from bsdlab.reduction import local_data
    ld = local_data(E, p)
    if ld.reduction == "good":
        # w_p = +1 at good reduction (Rohrlich; Rizzo Fact 3(2)).
        return ("ok", 1)
    if ld.reduction in ("split multiplicative", "nonsplit multiplicative"):
        # Classical: w_p = -1 iff split (Rizzo Fact 3(3)).  The split
        # test is made independently of Tate's algorithm:
        #   p odd  : split iff -c6 is a square mod p (Rizzo Lemma 7);
        #   p = 2  : split iff c6 = 1 mod 8            (Rizzo Lemma 7).
        c6m = _model_for(E, p).c6
        if p == 2:
            split = _unit(c6m, 2) % 8 == 1
        else:
            split = _legendre(-_unit(c6m, p), p) == 1
        tate_split = ld.reduction == "split multiplicative"
        if split != tate_split:
            raise ArithmeticError(
                "split test disagreement at p=%d: invariants say %s, "
                "Tate says %s" % (p, split, tate_split))
        return ("ok", -1 if split else 1)
    # additive reduction
    if p >= 5:
        return _additive_p_ge5(ld, p)
    if p == 3:
        return _additive_p3(_model_for(E, p), ld)
    return _additive_p2(_model_for(E, p), ld)


def _model_for(E, p):
    """The model bsdlab.reduction.local_data runs Tate's algorithm on,
    so that c4/c6 valuations here refer to the same model."""
    try:
        return E.minimal_model()
    except (ArithmeticError, ValueError):
        return E


def _additive_p_ge5(ld, p) -> Status:
    """Additive reduction at p >= 5: Halberstadt's classification.

    Halberstadt, C. R. Acad. Sci. Paris 326 (1998) 1047-1052 (as stated
    in Rohrlich, Compositio 87 (1993), Prop. 2, and tabulated in Rizzo,
    Compositio 136 (2003), Table I -- two independent published sources,
    which agree line by line).  For p >= 5 the entry depends only on the
    Kodaira type:

      II  : w_p = (-1/p)          [Table I: (>=1,1,2)]
      III : w_p = (-2/p)          [Table I: (1,>=2,3)]
      IV  : w_p = (-3/p)          [Table I: (>=2,2,4)]
      I0* : w_p = (-1/p)          [Table I: (2,>=3,6) and (>=2,3,6)]
      In* : w_p = (-1/p)          [additive potentially multiplicative,
                                   (2,3,>=7); also Fact 3(4): w_p = (-1/p)]
      IV* : w_p = (-3/p)          [Table I: (>=3,4,8)]
      III*: w_p = (-2/p)          [Table I: (3,>=5,9)]
      II* : w_p = (-1/p)          [Table I: (>=4,5,10)]
    """
    k = ld.kodaira
    if k in ("II", "II*", "I0*") or (k.endswith("*") and k.startswith("I") and
                                     k != "I0*" and k[1:-1].isdigit()):
        return ("ok", _legendre(-1, p))       # II, II*, I0*, In*
    if k in ("III", "III*"):
        return ("ok", _legendre(-2, p))       # III, III*
    if k in ("IV", "IV*"):
        return ("ok", _legendre(-3, p))       # IV, IV*
    raise ArithmeticError("unknown Kodaira symbol %r at p=%d" % (k, p))


def _triplet(Emin, p):
    """Rizzo's (a,b,c): smallest nonnegative triplet with
    (a,b,c) = (v_p(c4), v_p(c6), v_p(Delta)) mod (4,6,12), plus the raw
    valuations and unit parts, all on the minimal model Emin."""
    c4, c6, D = Emin.c4, Emin.c6, Emin.discriminant
    v4, v6, vD = _vp(c4, p), _vp(c6, p), _vp(D, p)
    a, b, c = v4 % 4, v6 % 6, vD % 12
    u4, u6 = _unit(c4, p), _unit(c6, p)
    # c_{n,e} = c_n / p^o  with  o = n*floor(v_p(c_n)/n) + e   (Rizzo p.2)
    o4 = 4 * (v4 // 4) + 2
    o6 = 6 * (v6 // 6) + 4
    return (a, b, c, v4, v6, vD, u4, u6,
            c4 // p ** min(o4, v4) if c4 else 0,
            c6 // p ** min(o6, v6) if c6 else 0)


def _additive_p3(Emin, ld) -> Status:
    """Additive reduction at p = 3: Kraus's table (Manuscripta Math. 69
    (1990)), as reproduced in Rizzo, Compositio 136 (2003), Table II.

    Only rows that could be decoded with certainty from the published
    table are implemented; each case cites its triplet (a,b,c).  Any
    other triplet/Kodaira combination returns ('none', reason) rather
    than a guessed value.
    """
    a, b, c, v4, v6, vD, u4, u6 = _triplet(Emin, 3)[:8]
    k = ld.kodaira
    # (2,3,>=6) I_{c-6}* and (3,>=6,6) I0*:  both rows give W3 = -1
    # (potentially multiplicative, Lemma 8(2), and the I0* row).
    if k == "I0*" or (k.endswith("*") and k[1:-1].isdigit()):
        return ("ok", -1)
    if k in ("III", "III*"):
        # Rows (>=2,3,3) [cond. c6'^2+2^6 = -3 c4'^2 (9)] and (2,>=5,3):
        # W3 = +1 unconditionally on both decoded III/III* rows.
        return ("ok", 1)
    if k == "II":
        # v(Delta)=3: row (>=2,3,3), W3 = +1 iff c6' = 4,7,8 mod 9.
        if c == 3 and b == 3 and a >= 2:
            return ("ok", 1 if u6 % 9 in (4, 7, 8) else -1)
        # v(Delta)=4: row (2,3,4): W3 = +1.
        if c == 4:
            return ("ok", 1)
        # v(Delta)=5: row (>=3,4,5): W3 = +1 iff c6' = -2 mod 3.
        if c == 5 and a >= 3:
            return ("ok", 1 if u6 % 3 == 1 else -1)
        return ("none", "p=3 additive II row not decoded (triplet %r)"
                % ((a, b, c),))
    return ("none", "p=3 additive table entry not decoded (triplet %r, %s)"
            % ((a, b, c), k))


def _additive_p2(Emin, ld) -> Status:
    """Additive reduction at p = 2: only what is PROPERLY derivable.

    Implemented cases (Halberstadt, CRAS 326 (1998); Connell via Rizzo
    Fact 3(4), Lemma 9; Rizzo Table III):
      - additive potentially multiplicative (v2(Delta) > 3*v2(c4)):
        W2 = +1 iff c6' = 3 mod 4 (Connell's formula), where the split
        sub-case c6' = 1 mod 4 with triplet (0,0,>=7) forces W2 = -1;
      - potentially good, Kodaira I3*, triplet (2,b,0) (12 | v2(Delta),
        v2(c4) = 2 mod 4, c_{4,2} = 1 mod 4): Rizzo Table III rows
          b = 4:  W2 = +1 iff c_{4,2} + 4 c_{6,4} = 9,13 (16)
          b >= 5 (incl. c6 = 0): W2 = +1 iff ... = 5,9 (16).
    Everything else: ('none', P2_TABLE_MISSING).
    """
    c4, c6, D = Emin.c4, Emin.c6, Emin.discriminant
    v4, v6, vD = _vp(c4, 2), _vp(c6, 2), _vp(D, 2)
    u4, u6 = _unit(c4, 2), _unit(c6, 2)
    if vD > 3 * v4:
        # Additive potentially multiplicative (Rizzo Fact 6 / Lemma 9(2)):
        # triplet (0,0,>=7) forces c6' = 1 mod 4 and W2 = -1; triplet
        # (2,3,>=7): W2 = +1 iff c6' = 3 mod 4 (Fact 3(4)).
        if v6 % 6 == 0:
            return ("ok", -1)   # (0,0,>=7) with c6' = 1 mod 4
        return ("ok", 1 if u6 % 4 == 3 else -1)
    if ld.kodaira == "I3*" and vD % 12 == 0 and v4 % 4 == 2 and u4 % 4 == 1:
        o6 = 6 * (v6 // 6) + 4          # c_{6,4} = c6/2^o6 (0 stays 0)
        c64 = c6 // 2 ** min(o6, v6) if c6 else 0
        lhs = (u4 + 4 * c64) % 16
        if c6 == 0 or v6 % 6 == 5:      # row (2,>=5,0): c6=0 counts as v=inf
            return ("ok", 1 if lhs in (5, 9) else -1)
        if v6 % 6 == 4:                 # row (2,4,0)
            return ("ok", 1 if lhs in (9, 13) else -1)
    return ("none", P2_TABLE_MISSING)


def global_root_number_two_route(E) -> dict:
    """Global root number by two independent routes.

    Returns {'table_route':  {'sign': +1/-1} or {'sign': None, 'reason': s},
             'numeric_route': the Fricke numeric sign (lseries), or None,
             'agree':         True/False, or None when either route missing}

    The table route is w = -prod_{p | N} w_p with every w_p from
    ``local_root_number``; a single missing local entry (e.g. an
    undecoded additive p = 2 table row) blanks the product with its
    reason -- it is never filled in by guesswork.
    """
    from bsdlab.reduction import conductor
    from bsdlab.lseries import root_number_numeric
    N = conductor(E)
    w, missing = -1, None
    for p in E.bad_primes():
        status, val = local_root_number(E, p)
        if status != "ok":
            missing = "p=%d: %s" % (p, val)
            break
        w *= val
    table = {"sign": None, "reason": missing} if missing else {"sign": w}
    try:
        numeric = root_number_numeric(E, N)
    except (ArithmeticError, ImportError):
        numeric = None
    agree = None
    if table["sign"] is not None and numeric is not None:
        agree = (table["sign"] == numeric)
        if not agree:
            raise ArithmeticError(
                "table route %d disagrees with numeric route %d for %r (N=%d)"
                % (table["sign"], numeric, E, N))
    return {"table_route": table, "numeric_route": numeric, "agree": agree}
