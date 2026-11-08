# Run from the repository root: python3 -I scripts/gt_blocked_check.py .   (about 6 minutes)
"""Exhaustive check of arith.gt_into_blocked against its contract, on every input.

Contract (as used in TwoBitMod.tick): t ^= ctrl*[a > b]; a, b, ctrl, the dirty wires (any value)
and the zeros wires are restored; zeros may hold anything when ctrl = 0 but read 0 when ctrl = 1;
fresh wires end at 0.
"""
import itertools, sys
sys.path.insert(0, sys.argv[1])
from qec.builder import Circuit
from qec.arith import gt_into_blocked

def compositions(w, m):
    for cuts in itertools.combinations(range(1, w), m - 1):
        b = (0,) + cuts + (w,)
        yield [b[i + 1] - b[i] for i in range(m)]

def run(w, sizes, n_dirty, n_zero):
    C = Circuit()
    a, b = C.alloc(w), C.alloc(w)
    (t, ctrl) = C.alloc(2)
    dirty, zeros = C.alloc(n_dirty), C.alloc(n_zero)
    first_fresh = C.num_wires
    gt_into_blocked(C, a, b, t, ctrl, sizes, zeros, dirty)
    nw = C.num_wires
    ops = C.ops
    # lanes: every (a, b, ctrl, t, dirty, zeros) with zeros = 0 when ctrl = 1
    lanes = []
    for av in range(1 << w):
        for bv in range(1 << w):
            for c in (0, 1):
                for tv in (0, 1):
                    for dv in range(1 << n_dirty):
                        for zv in (range(1 << n_zero) if c == 0 else (0,)):
                            lanes.append((av, bv, c, tv, dv, zv))
    L = len(lanes); ALL = (1 << L) - 1
    s = [0] * nw
    def setbits(wires, vals):
        for i, wire in enumerate(wires):
            m = 0
            for l, v in enumerate(vals):
                if v >> i & 1: m |= 1 << l
            s[wire] = m
    setbits(a, [x[0] for x in lanes]); setbits(b, [x[1] for x in lanes])
    setbits([ctrl], [x[2] for x in lanes]); setbits([t], [x[3] for x in lanes])
    setbits(dirty, [x[4] for x in lanes]); setbits(zeros, [x[5] for x in lanes])
    init = list(s)
    for k in range(0, len(ops), 4):
        o, x, y, z = ops[k:k + 4]
        if o == 0: s[x] ^= ALL
        elif o == 1: s[y] ^= s[x]
        elif o == 2: s[z] ^= s[x] & s[y]
        else: s[x], s[y] = s[y], s[x]
    want_t = 0
    for l, (av, bv, c, tv, dv, zv) in enumerate(lanes):
        if tv ^ (c & (av > bv)): want_t |= 1 << l
    bad = []
    if s[t] != want_t: bad.append("t")
    for wire in [*a, *b, ctrl, *dirty, *zeros]:
        if s[wire] != init[wire]: bad.append(f"wire {wire} not restored")
    for wire in range(first_fresh, nw):
        if s[wire]: bad.append(f"fresh {wire} dirty")
    ccx = ops[0::4].count(2)
    return bad, L, ccx

total = 0
for w in range(1, 7):
    for m in range(1, w + 1):
        for sizes in compositions(w, m):
            for n_dirty in range(0, min(3, m - 1) + 1):
                n_zero = min(2, 2 * m - 1 - n_dirty)
                for nz in range(0, n_zero + 1):
                    bad, L, ccx = run(w, sizes, n_dirty, nz)
                    total += 1
                    if bad:
                        print("FAIL", w, sizes, n_dirty, nz, bad[:3]); sys.exit(1)
print("all", total, "configurations pass")
