"""Time one component in isolation under the board's tick model (all wires free at tick 0).

    python3 -I -B iso.py <tree> <component> [n]

Components: add_la (Mike's add_or_shift_la), add_tb (the tree's TwoBitMod._add_or_shift),
diff (the tree's DialogMultiplier._shared_difference at w = n; ring of 193 if the tree takes one).
Prints CCX, raw ticks, and ticks after the tree's depth_pass if it has one.
"""
import sys
from pathlib import Path

tree = Path(sys.argv[1]).resolve()
comp = sys.argv[2]
n = int(sys.argv[3]) if len(sys.argv) > 3 else 256
sys.path.insert(0, str(tree))
from qec.builder import Circuit, CCX, CX, SWAP, X  # noqa: E402
from qec.curve import SECP256K1  # noqa: E402


def ticks(ops):
    r = {}
    ccx = 0
    for i in range(0, len(ops), 4):
        k, a, b, c = ops[i:i + 4]
        if k == X:
            r[a] = r.get(a, 0) + 1
        elif k == CX:
            t = max(r.get(a, 0), r.get(b, 0)) + 1
            r[a] = r[b] = t
        elif k == SWAP:
            t = max(r.get(a, 0), r.get(b, 0)) + 3
            r[a] = r[b] = t
        else:
            ccx += 1
            ra, rb, rt = r.get(a, 0), r.get(b, 0), r.get(c, 0)
            rt += 1
            x = max(rb, rt) + 1; rb = rt = x
            rt += 1
            x = max(ra, rt) + 1; ra = rt = x
            rt += 1
            x = max(rb, rt) + 1; rb = rt = x
            rt += 1
            x = max(ra, rt) + 1; ra = rt = x
            rb += 1; rt += 1
            rt += 1
            x = max(ra, rb) + 1; ra = rb = x
            ra += 1; rb += 1
            x = max(ra, rb) + 1; ra = rb = x
            r[a], r[b], r[c] = ra, rb, rt
    return ccx, max(r.values(), default=0)


def build(fn, sizes):
    C = Circuit(True)
    regs = [C.alloc(k) for k in sizes]
    blk = C.record(lambda *R: fn(C, *R), regs)
    return C, blk


if comp in ("add_la", "add_tb"):
    from qec import two_bit
    mod = two_bit.TwoBitMod(SECP256K1.q) if n == 256 else None
    if comp == "add_la":
        fn = lambda C, J, a, word: two_bit.add_or_shift_la(C, J[0], a, word, n)
    else:
        if hasattr(two_bit, "LOOKAHEAD_ADD"):
            two_bit.LOOKAHEAD_ADD = False
        fn = lambda C, J, a, word: mod._add_or_shift(C, J[0], a, word)
    sizes = [1, n, n + 2]
elif comp == "diff":
    from qec import point_add
    from qec.arith import Mod
    dm = point_add.DialogMultiplier.__new__(point_add.DialogMultiplier)
    ring = int(sys.argv[4]) if len(sys.argv) > 4 else 193
    import inspect
    takes_ring = "ring_wires" in inspect.signature(point_add.DialogMultiplier._shared_difference).parameters
    def fn(C, U, V, g, R):
        dm.C = C
        if takes_ring:
            dm._shared_difference(U, V, g[0], g[1], R)
        else:
            dm._shared_difference(U, V, g[0], g[1])
    sizes = [n, n, 2, ring]
else:
    raise SystemExit("unknown component")

C, blk = build(fn, sizes)
ccx, raw = ticks(blk.ops)
out = {"component": comp, "n": n, "ccx": ccx, "raw_ticks": raw}
try:
    from qec import depth_pass
    order, ids = [blk], {blk: 0}
    bodies = depth_pass.optimize(order, ids)
    ops, _ = bodies[blk]
    out["depth_pass_ticks"] = ticks(ops)[1]
    inv = []
    for i in range(len(ops) - 4, -1, -4):
        inv.extend(ops[i:i + 4])
    out["depth_pass_inverse_ticks"] = ticks(inv)[1]
except ImportError:
    pass
print(out)
