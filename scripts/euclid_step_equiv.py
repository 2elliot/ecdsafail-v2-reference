"""Exhaustive gate-level check of the Euclid step with the narrow-width ripple (port of
2elliot depth/v2's scripts/euclid_step_equiv.py to the v2-depth / bez-jack step shape).

    python3.12 -I -B scripts/euclid_step_equiv.py <tree> <wmax_all> <wmax_domain>

For each w, the step is _euclid_step(w, k) on U[:w], V, E, cell and a ring of k = w - 2 clean
wires (the ring PR8 hands every narrow step), plus k = 0 and k = 3 for coverage. Three builds:
  NEW  the tree as is (ripple for w <= RIPPLE_MAX_W),
  OLD  RIPPLE_MAX_W = 0, the variant PR8 uses for that k (u6 for k >= 7, u1 for k >= 3, else chain),
  REF  RIPPLE_MAX_W = 0 with k = 0, the chain fallback.
All inputs (w <= wmax_all): every assignment of U, V, E, cell (ring and scratch start 0). NEW must
equal REF on U, V, E, cell, the ring wires and the multiset of final scratch masks, on every lane.
Domain (w <= wmax_domain): cell = 0 and U odd, every assignment of the rest. NEW must equal OLD
on U, V, E, cell, and every ring and scratch wire must end 0 in both.
Exit 0 only if every comparison holds.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
from qec import point_add, sim  # noqa: E402
from qec.arith import Mod  # noqa: E402
from qec.builder import Circuit  # noqa: E402

WALL, WDOM = int(sys.argv[2]), int(sys.argv[3])
NEW_MAX = point_add.RIPPLE_MAX_W


def lane_mask(i, N):
    L = 1 << N
    period = 1 << (i + 1)
    block = ((1 << (1 << i)) - 1) << (1 << i)
    reps = ((1 << L) - 1) // ((1 << period) - 1)
    return reps * block


def build(w, k, ripple_max):
    point_add.RIPPLE_MAX_W = ripple_max
    q = next(p for p in (251, 509, 1021, 2039, 4093, 8191, 16381, 32749) if p.bit_length() >= w)
    C = Circuit()
    m = point_add.DialogMultiplier(C, Mod(q))
    U, V, E, cell, ring = C.alloc(w), C.alloc(w), C.alloc(1), C.alloc(2), C.alloc(k)
    blk = m._euclid_step(w, k) if k else m._euclid_step(w)
    C.apply(blk, [U, V, E, cell, ring] if k else [U, V, E, cell])
    point_add.RIPPLE_MAX_W = NEW_MAX
    return C, U, V, E, cell, ring


def run(w, k, ripple_max, domain):
    C, U, V, E, cell, ring = build(w, k, ripple_max)
    io = U + V + E + cell
    if domain:
        free = U[1:] + V + E
        N = len(free)
        full = (1 << (1 << N)) - 1
        init = {x: lane_mask(i, N) for i, x in enumerate(free)}
        init[U[0]] = full
    else:
        N = len(io)
        init = {x: lane_mask(i, N) for i, x in enumerate(io)}
    st = sim.run(C.ops, C.num_wires, init, 1 << N)
    scratch = [x for x in range(C.num_wires) if x not in set(io) | set(ring)]
    return ([st[x] for x in io], [st[x] for x in ring], sorted(st[x] for x in scratch if st[x]),
            C.ops[0::4].count(2), N)


ok = True
for w in range(2, max(WALL, WDOM) + 1):
    for k in sorted({0, 3, w - 2}):
        if k < 0:
            continue
        if w <= WALL:
            new = run(w, k, NEW_MAX, False)
            ref = run(w, 0, 0, False)
            same = new[0] == ref[0] and all(r == 0 for r in new[1]) and new[2] == ref[2]
            ok &= same
            print(f"all    w={w:2d} k={k:2d} lanes=2^{new[4]} ccx new/ref={new[3]}/{ref[3]} equal_to_chain={same}")
        if w <= WDOM:
            new = run(w, k, NEW_MAX, True)
            old = run(w, k, 0, True)
            same = new[0] == old[0] and not any(new[1]) and not new[2] and not any(old[1]) and not old[2]
            ok &= same
            print(f"domain w={w:2d} k={k:2d} lanes=2^{new[4]} ccx new/old={new[3]}/{old[3]} equal_to_pr8_and_clean={same}")
print("euclid step equivalence:", "OK" if ok else "FAIL")
sys.exit(0 if ok else 1)
