# Run from the repository root: python3 -I scripts/euclid_step_equiv.py . 8   (compare its hashes with another generator's)
"""Gate-level fingerprint of one Euclid step on EVERY input state (reachable or not).
    python3 -I step_equiv.py <generator dir> <wmax>
For w = 2..wmax builds a flat circuit holding only _euclid_step(w) applied to U[:w], V, E, cell, simulates all
2^(2w+3) assignments of those wires (scratch starts zero), and prints, per w: the SHA-256 of the output wires'
lane masks in output-layout order, whether every scratch wire ends zero in every lane, and the CCX count.
Two generators compute the same step function iff their hashes match."""
import hashlib, sys
sys.path.insert(0, sys.argv[1])
from qec import sim
from qec.arith import Mod
from qec.builder import Circuit
from qec.point_add import DialogMultiplier

def lane_mask(i, N):
    L = 1 << N
    period = 1 << (i + 1)
    block = ((1 << (1 << i)) - 1) << (1 << i)      # ones in the upper half of each period
    reps = ((1 << L) - 1) // ((1 << period) - 1)
    return reps * block

for w in range(2, int(sys.argv[2]) + 1):
    q = next(p for p in (251, 509, 1021, 2039, 4093) if p.bit_length() >= w)
    C = Circuit()
    m = DialogMultiplier(C, Mod(q))
    U, V, E, cell = C.alloc(w), C.alloc(w), C.alloc(1), C.alloc(2)
    ins = U + V + E + cell
    C.apply(m._euclid_step(w), [U, V, E, cell])
    outs = U + V + E + cell
    N = len(ins)
    init = {wire: lane_mask(i, N) for i, wire in enumerate(ins)}
    st = sim.run(C.ops, C.num_wires, init, 1 << N)
    h = hashlib.sha256()
    for wire in outs:
        h.update(st[wire].to_bytes(((1 << N) + 7) // 8, "little"))
    scratch = [x for x in range(C.num_wires) if x not in set(outs)]
    # The first scratch wire is `two` in both versions (the baseline's only one, ctl[0] here); fold
    # its final mask into the hash so the comparison covers it on every lane, reachable or not.
    h.update(st[scratch[0]].to_bytes(((1 << N) + 7) // 8, "little"))
    rest_zero = all(st[x] == 0 for x in scratch[1:])
    # On the step's domain (cell starts 0, U odd) every scratch wire must end zero.
    dom = ((1 << (1 << N)) - 1) & ~init[cell[0]] & ~init[cell[1]] & init[U[0]]
    clean_dom = all(st[x] & dom == 0 for x in scratch)
    print(f"w={w} lanes=2^{N} scratch={len(scratch)} others_zero_everywhere={rest_zero} "
          f"clean_on_domain={clean_dom} ccx={C.ops[0::4].count(2)} out+two={h.hexdigest()[:16]}")
