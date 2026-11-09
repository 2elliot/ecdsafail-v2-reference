"""Exhaustive check of point_add.eq_const_into and and_tree in tree form (point_add.EQ_TREE).

    python3.12 -I -B scripts/eq_tree_check.py <tree> <nmax> [broken]

For n = 1..nmax bits, every value < 2^n, with and without a control, on all inputs (bits, ctrl,
target): target must end target ^ ctrl*[bits == value], bits and ctrl restored, every other wire 0.
and_tree: for k = 1..nmax wires, the root must hold the AND and clear_and must leave every temp 0.
With "broken", and_tree drops its last CCX; the check must fail.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
from qec import point_add, sim
from qec.builder import Circuit

NMAX = int(sys.argv[2])
point_add.EQ_TREE = True
if len(sys.argv) > 3 and sys.argv[3] == "broken":
    orig = point_add.and_tree

    def broken(C, wires):
        root, history = orig(C, wires)
        if history:
            a, b, t = history[-1]
            C.ccx(a, b, t)
        return root, history
    point_add.and_tree = broken


def lane_mask(i, N):
    return sum(1 << lane for lane in range(1 << N) if lane >> i & 1)


bad = cases = 0
for n in range(1, NMAX + 1):
    for value in range(1 << n):
        for ctrl in (True, False):
            C = Circuit()
            c, bits, (t,) = C.alloc(1), C.alloc(n), C.alloc(1)
            point_add.eq_const_into(C, c[0] if ctrl else None, bits, value, t)
            ins = c + bits + [t]
            N = len(ins)
            init = {w: lane_mask(i, N) for i, w in enumerate(ins)}
            st = sim.run(C.ops, C.num_wires, init, 1 << N)
            for lane in range(1 << N):
                cv, tv = lane & 1, lane >> (n + 1) & 1
                bv = lane >> 1 & ((1 << n) - 1)
                want = tv ^ ((cv if ctrl else 1) & (bv == value))
                bad += (st[t] >> lane & 1) != want
            bad += any(st[w] != init[w] for w in c + bits)
            bad += any(st[w] for w in range(C.num_wires) if w not in init)
            cases += 1 << N
    C = Circuit()
    ws = C.alloc(n)
    root, history = point_add.and_tree(C, ws)
    (out,) = C.alloc(1)
    C.cx(root, out)
    point_add.clear_and(C, history)
    init = {w: lane_mask(i, n) for i, w in enumerate(ws)}
    st = sim.run(C.ops, C.num_wires, init, 1 << n)
    bad += st[out] != 1 << ((1 << n) - 1)
    bad += any(st[w] for w in range(C.num_wires) if w not in init and w != out)
print("cases", cases, "bad", bad)
print("OK" if bad == 0 else "FAIL")
sys.exit(0 if bad == 0 else 1)
