"""Exhaustive check of arith.add_lookahead, chain and tree forms, as the square's rows use it.

    python3.12 -I -B scripts/row_add_check.py <tree> <wmax> [broken]

For w = 1..wmax and temps in {1, 2, w, 2w}, on all 2^(2w + 1) inputs (a, b, cout): b must end
(a + b) mod 2^w, cout must end cout ^ carry out, a must be restored, and the temps and every other
wire must end 0. With "broken", the tree form's adjacent_xor runs its CXs in reverse order; the
check must fail.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
from qec import arith, sim
from qec.builder import Circuit

WMAX = int(sys.argv[2])
if len(sys.argv) > 3 and sys.argv[3] == "broken":
    orig = arith.adjacent_xor
    def broken(C, a):
        g = arith._Gates(); orig(g, a)
        for o in reversed(g.ops):
            C.cx(*o)
    arith.adjacent_xor = broken


def lane_mask(i, N):
    L = 1 << N
    period = 1 << (i + 1)
    block = ((1 << (1 << i)) - 1) << (1 << i)
    return ((1 << L) - 1) // ((1 << period) - 1) * block


def unpack(st, wires, lanes):
    vals = [0] * lanes
    for j, w in enumerate(wires):
        x = st[w]
        for lane in range(lanes):
            if x >> lane & 1:
                vals[lane] |= 1 << j
    return vals


ok = True
for tree in (False, True):
    for w in range(1, WMAX + 1):
        for nt in sorted({1, 2, w, 2 * w}):
            C = Circuit()
            a, b, (cout,), T = C.alloc(w), C.alloc(w), C.alloc(1), C.alloc(nt)
            arith.add_lookahead(C, a, b, cout, T, tree=tree)
            free = a + b + [cout]
            N = len(free)
            lanes = 1 << N
            init = {x: lane_mask(i, N) for i, x in enumerate(free)}
            st = sim.run(C.ops, C.num_wires, init, lanes)
            A, B, Co = unpack(init, a, lanes), unpack(init, b, lanes), unpack(init, [cout], lanes)
            A2, B2, Co2 = unpack(st, a, lanes), unpack(st, b, lanes), unpack(st, [cout], lanes)
            good = A2 == A and all(B2[i] == (A[i] + B[i]) % (1 << w) and Co2[i] == Co[i] ^ ((A[i] + B[i]) >> w)
                                   for i in range(lanes))
            good &= all(not st[x] for x in range(C.num_wires) if x not in init)
            ok &= good
            print(f"tree={tree} w={w} temps={nt} lanes=2^{N} ok={good}", flush=True)
print("row add check:", "OK" if ok else "FAIL")
sys.exit(0 if ok else 1)
