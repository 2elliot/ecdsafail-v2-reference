"""The rank-one start's wire and cell changes, checked exhaustively at toy sizes, with mutants that must fail.

  units    rank_one.inc (with a control and without) and dec on every input of 1 to 7 bits and
           both control values, every other wire 0 after
  residue  RankOne.residue on every U in [1, 2^w) and D in (-2^w, 2^w) for w = 2..6: D mod U, the
           sign and the quotient of |D|, every other wire 0 after
  fold     RankOne.fold maps X to (X mod q, floor(X/q)) on every X < 2^w * q with every other wire 0,
           at primes q = 2^n - f whose f splits as 2^t + c with zero wires between c's product and
           t (the wires the fold now lends), and at primes with no split
  M        DialogMultiplier's body maps (x, y) to (x, x*y mod q) with every other wire 0, on every
           x in [1, q) and y in [0, q), at primes whose pass takes the start (61: J = 3, 127: J = 6,
           251: J = 4, 509: J = 7)
  mutants  an incrementer without its lowest bit; one that clears its carries after the sum; a
           cell that swaps on its J wire; a J wire that holds b; and, for each free call in fold
           and convert, the same circuit with the first freed wire flipped just before, at q = 61
           and 127 and in each fold: a wire given back clean must matter, so every one has to fail

usage: start_wires_check.py [q ...]"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import rank_one
from qec.arith import Mod
from qec.builder import Circuit
from qec.point_add import DialogMultiplier


def run(C, init, lanes):
    """Bit-sliced run of a flat circuit: init maps wire -> lane mask."""
    full = (1 << lanes) - 1
    s = [0] * C.num_wires
    for w, v in init.items():
        s[w] = v
    p = C.ops
    for i in range(0, len(p), 4):
        k, a, b, c = p[i:i + 4]
        if k == 0:
            s[a] ^= full
        elif k == 1:
            s[b] ^= s[a]
        elif k == 2:
            s[c] ^= s[a] & s[b]
        else:
            s[a], s[b] = s[b], s[a]
    return s


def pack(values, wires):
    out = {w: 0 for w in wires}
    for lane, v in enumerate(values):
        for i, w in enumerate(wires):
            if v >> i & 1:
                out[w] |= 1 << lane
    return out


def differs(C, state, want):
    return any(state[w] != want.get(w, 0) for w in range(C.num_wires))


def unit_bad(k):
    """Lanes on which inc, dec or the zero test is wrong at k bits."""
    bad = 0
    vals = [(x, c) for x in range(1 << k) for c in (0, 1)]
    for name, f in (("inc", lambda x, c: (x + c) % (1 << k)), ("dec", lambda x, c: (x - c) % (1 << k))):
        C = Circuit()
        X, (ctl,) = C.alloc(k), C.alloc(1)
        getattr(rank_one, name)(C, X, ctl)
        init = pack([x for x, _ in vals], X) | pack([c for _, c in vals], [ctl])
        want = pack([f(x, c) for x, c in vals], X) | pack([c for _, c in vals], [ctl])
        bad += differs(C, run(C, init, len(vals)), want)
    C = Circuit()
    X = C.alloc(k)
    rank_one.inc(C, X)
    xs = list(range(1 << k))
    bad += differs(C, run(C, pack(xs, X), len(xs)), pack([(x + 1) % (1 << k) for x in xs], X))
    return bad


def residue_bad(w):
    """Whether RankOne.residue is wrong at w bits: on every U in [1, 2^w) and every D
    in (-2^w, 2^w), D mod U lands on D's low wires, the sign stays, Q takes floor(|D| / U)."""
    ro = rank_one.RankOne.__new__(rank_one.RankOne)
    ro.w = w
    C = Circuit()
    U, D, Q = C.alloc(w), C.alloc(w + 1), C.alloc(w)
    ro.residue(C, U, D, Q)
    vals = [(u, d) for u in range(1, 1 << w) for d in range(-(1 << w) + 1, 1 << w)]
    init = pack([u for u, _ in vals], U) | pack([d % (2 << w) for _, d in vals], D)
    want = pack([u for u, _ in vals], U) | pack([d % u | (d < 0) << w for u, d in vals], D) | pack(
        [abs(d) // u for u, d in vals], Q)
    return differs(C, run(C, init, len(vals)), want)


def m_bad(q, flip=None):
    """Whether M is wrong on some (x, y) at q. flip = i flips the first wire of the i-th nonempty
    free call made from fold or convert. Returns (bad, number of such calls)."""
    C = Circuit()
    mod = Mod(q)
    X, Y = C.alloc(mod.n), C.alloc(mod.n)
    seen = [0]
    free = C.free

    def spy(wires):
        if wires and sys._getframe(1).f_code.co_name in ("fold", "convert"):
            if seen[0] == flip:
                C.x(wires[0])
            seen[0] += 1
        free(wires)

    C.free = spy
    m = DialogMultiplier(C, mod)
    assert m._rank_j(True, False), f"q = {q} does not take the rank-one start"
    xin, yin = list(X), list(Y)
    # the body leaves its registers on the wires its lists name after
    m._body(X, Y)
    pairs = [(x, y) for x in range(1, q) for y in range(q)]
    init = pack([x for x, _ in pairs], xin) | pack([y for _, y in pairs], yin)
    want = pack([x for x, _ in pairs], X) | pack([x * y % q for x, y in pairs], Y)
    return differs(C, run(C, init, len(pairs)), want), seen[0]


def fold_cases():
    """(q, w, gap) for small primes whose fold has zero wires under t, and two with no split."""
    out, plain = [], []
    for n in range(9, 15):
        for q in range((1 << n) - 1, (1 << n) - (1 << n - 3), -2):
            if any(q % d == 0 for d in range(3, int(q ** 0.5) + 1, 2)):
                continue
            f = (1 << n) - q
            for w in (2, 3, 4):
                if not ((2 ** w + 2) * f <= 2 ** n and ((2 ** w) * f).bit_length() <= n):
                    continue
                cut = rank_one.split(f, w + 1)
                if cut is None:
                    if len(plain) < 2:
                        plain.append((q, w, 0))
                    continue
                gap = cut[1] - (cut[0] * ((1 << w + 1) - 1)).bit_length()
                if gap > 0 and cut[0] and len([c for c in out if c[2] == gap]) < 2 and len(out) < 6:
                    out.append((q, w, gap))
    return out + plain


def fold_bad(q, w, flip=None):
    C = Circuit()
    ro = rank_one.RankOne(q, 3, w)
    n = ro.n
    X = C.alloc(n + w)
    seen = [0]
    free = C.free

    def spy(wires):
        if wires and sys._getframe(1).f_code.co_name == "fold":
            if seen[0] == flip:
                C.x(wires[0])
            seen[0] += 1
        free(wires)

    C.free = spy
    ro.fold(C, X)
    xs = list(range(q << w))
    want = pack([x % q | (x // q) << n for x in xs], X)
    return differs(C, run(C, pack(xs, X), len(xs)), want), seen[0]


def mutant_inc_low(C, x, ctl=None):
    if ctl is None:
        C.x(x[0])
        return
    k = len(x)
    cs = [ctl] + C.alloc(k - 1)
    for i in range(1, k):
        C.ccx(cs[i - 1], x[i - 1], cs[i])
    for i in reversed(range(1, k)):
        C.cx(cs[i], x[i])
        C.ccx(cs[i - 1], x[i - 1], cs[i])
    C.free(cs[1:])


def mutant_inc_order(C, x, ctl=None):
    if ctl is None:
        (ctl,) = C.alloc(1)
    k = len(x)
    cs = [ctl] + C.alloc(k - 1)
    for i in range(1, k):
        C.ccx(cs[i - 1], x[i - 1], cs[i])
    for i in range(1, k):
        C.cx(cs[i], x[i])
    C.cx(ctl, x[0])
    for i in reversed(range(1, k)):
        C.ccx(cs[i - 1], x[i - 1], cs[i])
    C.free(cs[1:])


def mutant_kind(C, g0, b):
    (k,) = C.alloc(1)
    C.cx(b, k)
    return k




def main():
    primes = [int(a) for a in sys.argv[1:]] or [61, 127, 251, 509]
    ok = True

    def report(label, good, expect_fail=False):
        nonlocal ok
        passed = (not good) if expect_fail else good
        ok &= passed
        print(f"{'ok  ' if passed else 'FAIL'} {label}", flush=True)

    report("units: inc and dec on 1..7 bits", not any(unit_bad(k) for k in range(1, 8)))
    report("residue on every U and D at w = 2..6", not any(residue_bad(w) for w in range(2, 7)))
    for q in primes:
        bad, frees = m_bad(q)
        report(f"M exact on every (x, y) at q = {q} ({frees} free calls in fold and convert)", not bad)
    cases = fold_cases()
    assert sum(1 for c in cases if c[2]) >= 3
    for q, w, gap in cases:
        bad, frees = fold_bad(q, w)
        report(f"fold exact on every X < 2^{w} * q at q = {q} ({gap} lent wires)", not bad)
        alive = [i for i in range(frees) if not fold_bad(q, w, flip=i)[0]]
        report(f"mutants, a dirty wire at each of {frees} free calls in that fold: all fail"
               + (f" (alive: {alive})" if alive else ""), bool(alive), expect_fail=True)
    for name, attr, owner, mutant in (("incrementer without its low bit", "inc", rank_one, mutant_inc_low),
                                      ("incrementer clearing after the sum", "inc", rank_one, mutant_inc_order),
                                      ("a swap on the J wire", "_kind_swap", rank_one, lambda C, b, k: k),
                                      ("a J wire that is b", "_kind", rank_one, mutant_kind)):
        good = owner.__dict__[attr]
        setattr(owner, attr, mutant)
        try:
            unit = attr.startswith("_kind") or any(unit_bad(k) for k in range(2, 8)) or any(residue_bad(w) for w in (2, 3))
            whole = any(m_bad(q)[0] for q in primes[:2])
        finally:
            setattr(owner, attr, good)
        if not attr.startswith("_kind"):
            report(f"mutant, {name}: units fail", not unit, expect_fail=True)
        report(f"mutant, {name}: M fails", not whole, expect_fail=True)
    for q in primes[:2]:
        _, frees = m_bad(q)
        alive = [i for i in range(frees) if not m_bad(q, flip=i)[0]]
        report(f"mutants, a dirty wire at each of {frees} free calls in fold and convert, q = {q}: all fail"
               + (f" (alive: {alive})" if alive else ""), bool(alive), expect_fail=True)
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
