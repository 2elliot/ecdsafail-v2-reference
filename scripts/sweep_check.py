"""Exhaustive checks for lookahead_tree at small blocks, and for the Euclid step built with them.

    python3.12 -I -B scripts/sweep_check.py <tree> <nmax> <wmax> [broken]

1. lookahead_tree(U, V, temps, r) for r = 1, 2, 3, n = 2..nmax, on all 2^(2n - 1) inputs (U[0] and
   U[j], V[j] for j >= 1): U must end f_j = x_j ^ p_j f_(j-1) from f_0 = U[0], V unchanged, temps 0.
2. The step _euclid_step(w, w - 2) with RIPPLE_MAX_W = 0 (the tree difference) for w = 9..wmax, on the
   step domain (cell 0, U odd, every U[1:], V, E): the tree's sweep sizes against the old (3, 4, 5, 6, 8),
   equal on U, V, E, cell, ring and scratch all 0; and the same step on a ring of 2w wires.
With "broken", lookahead_tree drops its last fix CCX; the checks must fail.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
from qec import arith, point_add, sim
from qec.arith import Mod
from qec.builder import Circuit

NMAX, WMAX = int(sys.argv[2]), int(sys.argv[3])
if len(sys.argv) > 4 and sys.argv[4] == "broken":
    orig_tree = arith._lookahead_tree
    class Drop:
        def __init__(self, C): self.C, self.last = C, None
        def cx(self, *a): self.flush(); self.C.cx(*a)
        def ccx(self, *a): self.flush(); self.last = a
        def flush(self):
            if self.last: self.C.ccx(*self.last); self.last = None
    def broken(C, U, V, q, r):
        d = Drop(C); orig_tree(d, U, V, q, r)  # the final CCX is never flushed
    arith._lookahead_tree = broken


def lane_mask(i, N):
    L = 1 << N
    period = 1 << (i + 1)
    block = ((1 << (1 << i)) - 1) << (1 << i)
    return ((1 << L) - 1) // ((1 << period) - 1) * block


ok = True
for r in (1, 2, 3):
    for n in range(2, NMAX + 1):
        C = Circuit()
        U, V, T = C.alloc(n), C.alloc(n), C.alloc(arith.tree_temps(n, r))
        arith.lookahead_tree(C, U, V, T, r)
        free = U + V[1:]
        N = len(free)
        init = {x: lane_mask(i, N) for i, x in enumerate(free)}
        st = sim.run(C.ops, C.num_wires, init, 1 << N)
        full = (1 << (1 << N)) - 1
        f = init[U[0]]
        good = st[U[0]] == f
        for j in range(1, n):
            f = init[U[j]] ^ (init[V[j]] & f)
            good &= st[U[j]] == f and st[V[j]] == init[V[j]]
        good &= not any(st[t] for t in T) and all(not st[x] for x in range(C.num_wires) if x not in init)
        ok &= good
        print(f"tree r={r} n={n:2d} lanes=2^{N} temps={len(T)} ok={good}", flush=True)

OLD = (3, 4, 5, 6, 8)


def step_run(w, sizes, k=None):
    saved = arith.SWEEP_BLOCKS
    arith.SWEEP_BLOCKS = sizes
    point_add.RIPPLE_MAX_W = 0
    q = next(p for p in (251, 509, 1021, 2039, 4093, 8191, 16381, 32749) if p.bit_length() >= w)
    C = Circuit()
    m = point_add.DialogMultiplier(C, Mod(q))
    k = w - 2 if k is None else k
    U, V, E, cell, ring = C.alloc(w), C.alloc(w), C.alloc(1), C.alloc(2), C.alloc(k)
    C.apply(m._euclid_step(w, k), [U, V, E, cell, ring])
    arith.SWEEP_BLOCKS = saved
    free = U[1:] + V + E
    N = len(free)
    init = {x: lane_mask(i, N) for i, x in enumerate(free)}
    init[U[0]] = (1 << (1 << N)) - 1
    st = sim.run(C.ops, C.num_wires, init, 1 << N)
    io = U + V + E + cell
    rest = [x for x in range(C.num_wires) if x not in set(io)]
    return [st[x] for x in io], not any(st[x] for x in rest), C.ops[0::4].count(2), N


for w in range(9, WMAX + 1):
    new, new_clean, new_ccx, N = step_run(w, arith.SWEEP_BLOCKS)
    old, old_clean, old_ccx, _ = step_run(w, OLD)
    good = new == old and new_clean and old_clean
    ok &= good
    print(f"step w={w:2d} lanes=2^{N} ccx new/old={new_ccx}/{old_ccx} equal_and_clean={good}", flush=True)
    wide, wide_clean, wide_ccx, _ = step_run(w, arith.SWEEP_BLOCKS, 2 * w)
    good = wide == old and wide_clean
    ok &= good
    print(f"step w={w:2d} ring={2 * w} (more than the tree needs) ccx={wide_ccx} equal_and_clean={good}", flush=True)
print("sweep check:", "OK" if ok else "FAIL")
sys.exit(0 if ok else 1)
