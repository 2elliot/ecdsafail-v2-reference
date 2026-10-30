"""Exhaustive or sampled basis check of the controlled point-addition circuit.

For every lane (c, R): after the circuit the accumulator must hold R + c*A,
the control must be unchanged, and every other wire must be zero."""
import argparse
import io
import json
import random
import sys
import time
from array import array
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import flatten, sim
from qec.builder import CCX, CX, SWAP, X
from qec.curve import SECP256K1, SECP256K1_G, toy_secp_curves
from qec.emit import write_circuit
from qec.point_add import ControlledPointAdd, encode


def structural(ops, num_wires):
    """Every wire in range and every multi-wire gate on distinct wires."""
    it = iter(ops)
    for op, a, b, c in zip(it, it, it, it):
        assert 0 <= a < num_wires and 0 <= b < num_wires and 0 <= c < num_wires
        if op == CCX:
            assert a != b and a != c and b != c
        elif op != X:
            assert a != b
    kinds = ops[0::4]
    return {name: kinds.count(k) for k, name in ((X, "X"), (CX, "CX"), (CCX, "CCX"), (SWAP, "SWAP"))}


def check(gen, ops, cases):
    """cases: list of (c, R). Returns list of failure descriptions."""
    C, lay, n, curve, A = gen.C, gen.layout, gen.n, gen.curve, gen.A
    acc = lay["x"] + lay["y"]
    lanes = len(cases)
    init = sim.pack([cr[0] for cr in cases], lay["control"])
    init.update(sim.pack([encode(R, n) for _, R in cases], acc))
    state = sim.run(ops, C.num_wires, init, lanes)
    out = sim.unpack(state, acc, lanes)
    ctl = sim.unpack(state, lay["control"], lanes)
    live = set(acc) | set(lay["control"])
    dirty = 0
    for w in range(C.num_wires):
        if w not in live:
            dirty |= state[w]
    fails = []
    for j, (c, R) in enumerate(cases):
        want = encode(curve.add(R, A) if c else R, n)
        if out[j] != want or ctl[j] != c or dirty >> j & 1:
            fails.append({"c": c, "R": None if R is None else [hex(R[0]), hex(R[1])], "got": hex(out[j]), "want": hex(want), "dirty": bool(dirty >> j & 1)})
    return fails


def build_ops(gen, hier):
    if not hier:
        C = gen.build()
        return C, C.ops
    C = gen.build(hierarchical=True)
    buf = io.StringIO()
    write_circuit(buf, C, gen.layout)
    codes = {"X": X, "CX": CX, "CCX": CCX, "SWAP": SWAP}
    ops = array("i")
    for chunk in flatten.operations_json(flatten.load(io.StringIO(buf.getvalue()))):
        for op in json.loads("[" + chunk.strip("[],") + "]"):
            wires = [w for k, w in op.items() if k != "op"]
            ops.extend([codes[op["op"]], *wires, *wires[-1:] * (3 - len(wires))])
    return C, ops


def special_points(curve, A):
    twoA = curve.add(A, A)
    return [None, A, curve.neg(A), curve.neg(twoA), twoA, curve.add(twoA, A)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--toy-bits", type=int, nargs="*", default=[])
    ap.add_argument("--toy-per-size", type=int, default=2)
    ap.add_argument("--secp", action="store_true")
    ap.add_argument("--secp-random", type=int, default=64)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--hier", action="store_true", help="simulate the flattened hierarchical emission")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    receipts = []

    for bits in args.toy_bits:
        for curve in toy_secp_curves(bits)[: args.toy_per_size]:
            pts = curve.points()
            for A in [P for P in pts[1:] if P[1] != 0][:2]:
                t0 = time.time()
                gen = ControlledPointAdd(curve, A)
                C, ops = build_ops(gen, args.hier)
                counts = structural(ops, C.num_wires)
                cases = [(c, R) for c in (0, 1) for R in pts]
                fails = check(gen, ops, cases)
                r = {"mode": "EXHAUSTIVE_PROMISE_DOMAIN", "q": curve.q, "curve_order": len(pts), "addend": list(A), "n": gen.n,
                     "wires": C.num_wires, "ops": len(ops) // 4, "counts": counts, "lanes": len(cases), "failures": len(fails),
                     "first_failures": fails[:3], "seconds": round(time.time() - t0, 2)}
                print(json.dumps(r), flush=True)
                receipts.append(r)

    if args.secp:
        t0 = time.time()
        curve, A = SECP256K1, SECP256K1_G
        gen = ControlledPointAdd(curve, A)
        C, ops = build_ops(gen, args.hier)
        counts = structural(ops, C.num_wires)
        t_build = time.time() - t0
        rand = [curve.mul(rng.randrange(1, 2**256), A) for _ in range(args.secp_random)]
        cases = [(c, R) for c in (0, 1) for R in special_points(curve, A) + rand]
        t1 = time.time()
        fails = check(gen, ops, cases)
        r = {"mode": "SAMPLED_256", "q": hex(curve.q), "addend": "secp256k1 G", "n": gen.n, "wires": C.num_wires, "ops": len(ops) // 4,
             "counts": counts, "lanes": len(cases), "special_cases": ["inf", "A", "-A", "-2A", "2A", "3A"], "failures": len(fails),
             "first_failures": fails[:3], "build_seconds": round(t_build, 1), "sim_seconds": round(time.time() - t1, 1)}
        print(json.dumps(r), flush=True)
        receipts.append(r)

    if args.out:
        args.out.write_text(json.dumps(receipts, indent=1))
    sys.exit(1 if any(r["failures"] for r in receipts) else 0)


if __name__ == "__main__":
    main()
