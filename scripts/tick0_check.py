"""Exhaustive equivalence of the first Euclid step in tree form (point_add.TICK0_TREE) with the chain form.

    python3.12 -I -B scripts/tick0_check.py <tree> <nmax> [broken]

For every odd q with 4 <= bit length n <= nmax, both forms of _euclid_tick0 run on U = q, every
n-bit V, E = 0, a clean cell and clean scratch (the domain the pass gives it), with 1, 3 and 31 copies
of g0. The tree form must leave U, V, E and the cell as the chain form does, and every scratch wire 0.
At secp256k1 the tree block must use no wire outside its registers, or the pass would skip it.
With "broken", the tree form's rotation reverses one pair too few; the check must fail.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
from qec import point_add, sim
from qec.arith import Mod
from qec.builder import Circuit
from qec.curve import SECP256K1

NMAX = int(sys.argv[2])
if len(sys.argv) > 3 and sys.argv[3] == "broken":
    orig = point_add._par_rotate
    point_add._par_rotate = lambda C, two, V, cp: orig(C, two, V[:-1], cp)


def run(q, tree, copies):
    n = q.bit_length()
    point_add.TICK0_TREE, point_add.TICK0_COPIES, point_add.TICK0_SCRATCH = tree, copies, 4 * n + 2 * copies + 8
    C = Circuit()
    m = point_add.DialogMultiplier(C, Mod(q))
    blk = m._euclid_tick0()
    k = m._tick0_scratch()
    U, V, E, cell, scr = C.alloc(n), C.alloc(n), C.alloc(1), C.alloc(2), C.alloc(k)
    regs = [U, V, E, cell, scr]
    C.apply(blk, regs)
    lanes = 1 << n
    init = {w: (lanes and (1 << lanes) - 1) if q >> i & 1 else 0 for i, w in enumerate(U)}
    for i, w in enumerate(V):
        init[w] = sum(1 << v for v in range(lanes) if v >> i & 1)
    st = sim.run(C.ops, C.num_wires, init, lanes)
    out = [st[w] for w in U + V + E + cell]
    live = set(U + V + E + cell)
    dirty = sum(st[w] != 0 for w in range(C.num_wires) if w not in live)
    return out, dirty


bad = cases = 0
for q in range(9, 1 << NMAX, 2):
    ref, d0 = run(q, False, 1)
    bad += d0
    for copies in (1, 3, 31):
        got, d1 = run(q, True, copies)
        bad += (got != ref) + d1
        cases += 1 << q.bit_length()
point_add.TICK0_TREE, point_add.TICK0_COPIES, point_add.TICK0_SCRATCH = True, 31, 450
C = Circuit(True)
m = point_add.DialogMultiplier(C, Mod(SECP256K1.q))
t = m._euclid_tick0()
contained = t.wires <= set(t.flat_in) and set(t.flat_out) == set(t.flat_in)
bad += not contained
print("cases", cases, "bad", bad, "secp256k1 block self-contained", contained, "scratch", m._tick0_scratch(),
      "tape at step 0", m.tape_external - 2)
print("OK" if bad == 0 else "FAIL")
sys.exit(0 if bad == 0 else 1)
