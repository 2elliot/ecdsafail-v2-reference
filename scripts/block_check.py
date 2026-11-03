"""The multiplier with the tape code, on inputs chosen for it, in both directions.

Runs M on every x at primes where two or three blocks fit, on the merged path, where the code
reads the B1-swept tape last cell first, and off it. At secp256k1
it runs structured x. Each case runs M forward from a clean start and
M's reverse from the clean expected output. A control drops one CCX of the code and must fail.
Usage: block_check.py <census inputs file>
"""
import json
import random
import sys
from array import array
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim, tape_block
from qec.arith import Mod
from qec.builder import Circuit
from qec.point_add import DialogMultiplier

SECP = 2**256 - 2**32 - 977


def inverse(ops):
    rev = ops[::-1]
    out = array("i", ops)
    out[0::4], out[1::4], out[2::4], out[3::4] = rev[3::4], rev[2::4], rev[1::4], rev[0::4]
    return out


def pack(values, wires):
    """Lane masks per wire in time linear in the lanes; sim.pack is quadratic at 2**20 lanes."""
    return {w: int("".join("1" if v >> i & 1 else "0" for v in reversed(values)) or "0", 2) for i, w in enumerate(wires)}


def bad_lanes(state, want, num_wires):
    """Count lanes where any wire differs from want, which is zero outside its keys."""
    bad = 0
    for w in range(num_wires):
        bad |= state[w] ^ want.get(w, 0)
    return bad.bit_count()


def run(q, cases):
    mod = Mod(q)
    C = Circuit()
    X, Y = C.alloc(mod.n), C.alloc(mod.n)
    X0, Y0 = list(X), list(Y)
    m = DialogMultiplier(C, mod)
    m.apply(X, Y)
    lanes = len(cases)
    xs, ys = [c[0] for c in cases], [c[1] for c in cases]
    want = [a * b % q if a else b for a, b in cases]
    init = pack(xs, X0) | pack(ys, Y0)
    out = pack(xs, X) | pack(want, Y)
    forward = bad_lanes(sim.run(C.ops, C.num_wires, init, lanes), out, C.num_wires)
    backward = bad_lanes(sim.run(inverse(C.ops), C.num_wires, out, lanes), init, C.num_wires)
    return {"q": hex(q), "merged": m.narrow.merged, "T": m.iters, "blocks": m.blocks, "u_lent": m.u_lent,
            "lanes": lanes, "wires": C.num_wires, "forward failures": forward, "inverse failures": backward}


def primes(lo, hi):
    return [p for p in range(lo, hi) if p % 2 and all(p % d for d in range(3, int(p**.5) + 1, 2))]


def merged_below(bits, count):
    out = []
    for p in reversed(primes(1 << (bits - 1), 1 << bits)):
        if DialogMultiplier(Circuit(), Mod(p)).narrow.merged:
            out.append(p)
            if len(out) == count:
                return out


def off_path(bits, r, count):
    pool = [p for p in r.sample(primes(1 << (bits - 1), 1 << bits), 40) if not DialogMultiplier(Circuit(), Mod(p)).narrow.merged]
    return pool[:count]


def main():
    r = random.Random(20261102)
    out = []
    # Two blocks at 20 and 21 bits, including merged and fallback paths.
    cases = [1048573, 700001, 2097143]
    for q in cases:
        out.append(run(q, [(x, r.randrange(q)) for x in range(q)]))
    q = 1073741789  # 30-bit prime: three blocks; sampled below.
    xs = sorted({0, 1, 2, q - 1, q - 2} | {r.randrange(q) for _ in range(200000)})
    out.append(run(q, [(x, r.randrange(q)) for x in xs]))
    xs = [int(l, 16) for l in Path(sys.argv[1]).read_text().split()]
    xs = r.sample(xs, 8000) + [0, 1, 2, SECP - 1, SECP - 2]
    out.append(run(SECP, [(x, r.randrange(SECP)) for x in xs]))
    # Control: one CCX fewer in the code must be seen.
    prims = list(tape_block.PRIMS)
    k = next(i for i, p in enumerate(prims) if p[0] == "ccx")
    tape_block.PRIMS = tuple(prims[:k] + prims[k + 1:])
    q = out[0]["q"]
    control = run(int(q, 16), [(x, r.randrange(int(q, 16))) for x in range(int(q, 16))])
    out.append({"control": control, "seen": control["forward failures"] > 0})
    print(json.dumps(out, indent=0))
    ok = all(o.get("forward failures", 0) == 0 and o.get("inverse failures", 0) == 0 for o in out[:-1])
    return 0 if ok and out[-1]["seen"] else 1


if __name__ == "__main__":
    sys.exit(main())
