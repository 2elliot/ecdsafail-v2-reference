# Run from the repository root: python3 -I scripts/crotate_check.py .
"""Exhaustive check of arith.crotate against the bubble chain it replaces.
For every m in 0..10, every k in 1..6, every control and data value, with copies starting at 0:
crotate must give exactly what the old chain `for i in range(m - 1): cswap(ctrl, a[i], a[i+1])` gives
(state equality on every wire, so the control is unchanged and every copy is back to zero), and with the
control set it must be the rotation a[i] <- a[(i + 1) % m]."""
import sys
sys.path.insert(0, sys.argv[1])
from qec.builder import Circuit
from qec.arith import crotate, cswap

def sim(ops, state):
    s = list(state)
    for k in range(0, len(ops), 4):
        o, a, b, t = ops[k:k + 4]
        if o == 0: s[a] ^= 1
        elif o == 1: s[b] ^= s[a]
        elif o == 2: s[t] ^= s[a] & s[b]
        else: s[a], s[b] = s[b], s[a]
    return s

checked = 0
for m in range(0, 11):
    for k in range(1, 7):
        C = Circuit(); ctrls = C.alloc(k); a = C.alloc(m)
        crotate(C, ctrls, a)
        R = Circuit(); rc = R.alloc(k); ra = R.alloc(m)
        for i in range(m - 1):
            cswap(R, rc[0], [ra[i]], [ra[i + 1]])
        ccx_new, ccx_old = C.prims[0::4].count(2), R.prims[0::4].count(2)
        assert ccx_new == ccx_old == max(m - 1, 0), (m, k, ccx_new, ccx_old)
        for c in (0, 1):
            for v in range(1 << m):
                for junk in ((0,) * (k - 1),):
                    st = [c, *junk] + [v >> i & 1 for i in range(m)]
                    got = sim(C.prims, st)
                    want = sim(R.prims, st)
                    assert got == want, (m, k, c, v, junk, got, want)
                    assert got[:k] == [c] + [0] * (k - 1)
                    if c and m:
                        exp = [v >> ((i + 1) % m) & 1 for i in range(m)]
                        assert got[k:] == exp
                    checked += 1
print("crotate exhaustive: OK", checked, "cases (m 0..10, k 1..6, control 0/1, clean copies)")
