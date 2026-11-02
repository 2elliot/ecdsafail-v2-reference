"""Exhaustive or sampled basis check of the controlled point-addition circuit.

For every lane (c, R): after the circuit the accumulator must hold R + c*A,
the control must be unchanged, and every other wire must be zero."""
import argparse
import hashlib
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
from qec.emit import _op, write_circuit
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


def check(gen, ops, cases, num_wires=None):
    """cases: list of (c, R). Returns list of failure descriptions. num_wires
    defaults to the built circuit's; pass a file's when ops came from it."""
    lay, n, curve, A = gen.layout, gen.n, gen.curve, gen.A
    num_wires = num_wires or gen.C.num_wires
    acc = lay["x"] + lay["y"]
    lanes = len(cases)
    init = sim.pack([cr[0] for cr in cases], lay["control"])
    init.update(sim.pack([encode(R, n) for _, R in cases], acc))
    state = sim.run(ops, num_wires, init, lanes)
    out = sim.unpack(state, acc, lanes)
    ctl = sim.unpack(state, lay["control"], lanes)
    live = set(acc) | set(lay["control"])
    dirty = 0
    for w in range(num_wires):
        if w not in live:
            dirty |= state[w]
    fails = []
    for j, (c, R) in enumerate(cases):
        want = encode(curve.add(R, A) if c else R, n)
        if out[j] != want or ctl[j] != c or dirty >> j & 1:
            fails.append({"c": c, "R": None if R is None else [hex(R[0]), hex(R[1])], "got": hex(out[j]), "want": hex(want), "dirty": bool(dirty >> j & 1)})
    return fails


def expand(doc):
    """The flat ops of a circuit that flatten.load parsed, in the builder's
    (op, a, b, c) layout. It is flatten.operations_json's macro expansion,
    memoized per definition: a call substitutes its wire map into the callee,
    and an inverse call runs the callee reversed with nested flags flipped."""
    code = {flatten._templates[name]: k for name, k in (("X", X), ("CX", CX), ("CCX", CCX), ("SWAP", SWAP))}

    def segments(body):
        segs, run = [], None
        for ins in body:
            if type(ins) is flatten.Call:
                segs.append(ins)
                run = None
                continue
            if run is None:
                run = (array("b"), array("i"), array("i"), array("i"))
                segs.append(run)
            w = ins.wires
            run[0].append(code[ins.template])
            run[1].append(w[0])
            run[2].append(w[1] if len(w) > 1 else w[0])
            run[3].append(w[-1])
        return segs

    forward = [segments(d["operations"]) for d in doc["definitions"]]
    backward = {}

    def body(callee, inverse):
        if not inverse:
            return forward[callee]
        if callee not in backward:
            backward[callee] = [s if type(s) is flatten.Call else tuple(a[::-1] for a in s) for s in reversed(forward[callee])]
        return backward[callee]

    kinds, slots = array("b"), (array("i"), array("i"), array("i"))

    def emit(segs, wires, inverse):
        pick = wires.__getitem__
        for s in segs:
            if type(s) is flatten.Call:
                flip = s.inverse != inverse
                emit(body(s.callee, flip), [wires[w] for w in s.wires], flip)
            else:
                kinds.extend(s[0])
                for out, local in zip(slots, s[1:]):
                    out.extend(map(pick, local))

    emit(segments(doc["operations"]), list(range(doc["num_wires"])), False)
    ops = array("i", bytes(16 * len(kinds)))
    ops[0::4] = array("i", kinds)
    for k, out in enumerate(slots, 1):
        ops[k::4] = out
    return ops


def render(ops):
    it = iter(ops)
    return "[" + ",".join(_op(*t) for t in zip(it, it, it, it)) + "]"


def build_ops(gen, hier):
    """The circuit and its flat ops. For hier, the ops come from the emitted JSON
    text, parsed and validated by qec/flatten.py and expanded, and the third
    value is that parsed document and the text's sha256."""
    if not hier:
        C = gen.build()
        return C, C.ops, None
    C = gen.build(hierarchical=True)
    buf = io.StringIO()
    write_circuit(buf, C, gen.layout)
    text = buf.getvalue()
    doc = flatten.load(io.StringIO(text))
    return C, expand(doc), (doc, hashlib.sha256(text.encode()).hexdigest())


def special_points(curve, A):
    twoA = curve.add(A, A)
    return [None, A, curve.neg(A), curve.neg(twoA), twoA, curve.add(twoA, A)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--toy-bits", type=int, nargs="*", default=[])
    ap.add_argument("--toy-per-size", type=int, default=2)
    ap.add_argument("--secp", action="store_true")
    ap.add_argument("--secp-random", type=int, default=64)
    ap.add_argument("--secp-addend", type=lambda s: int(s, 0), nargs="+", default=[1], metavar="K", help="check A = [K]G for each K (default 1)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--hier", action="store_true", help="simulate the emitted hierarchical JSON, expanded by qec/flatten.py's rules")
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
                C, ops, emitted = build_ops(gen, args.hier)
                counts = structural(ops, C.num_wires)
                cases = [(c, R) for c in (0, 1) for R in pts]
                fails = check(gen, ops, cases)
                r = {"mode": "EXHAUSTIVE_PROMISE_DOMAIN", "q": curve.q, "curve_order": len(pts), "addend": list(A), "n": gen.n,
                     "wires": C.num_wires, "ops": len(ops) // 4, "counts": counts, "lanes": len(cases), "failures": len(fails),
                     "first_failures": fails[:3], "seconds": round(time.time() - t0, 2)}
                if emitted:
                    r["expansion_matches_flatten"] = render(ops) == "".join(flatten.operations_json(emitted[0]))
                print(json.dumps(r), flush=True)
                receipts.append(r)

    for k in args.secp_addend if args.secp else ():
        t0 = time.time()
        curve = SECP256K1
        A = SECP256K1_G if k == 1 else curve.mul(k, SECP256K1_G)
        gen = ControlledPointAdd(curve, A)
        C, ops, emitted = build_ops(gen, args.hier)
        sha = emitted and emitted[1]
        del emitted
        counts = structural(ops, C.num_wires)
        t_build = time.time() - t0
        rand = [curve.mul(rng.randrange(1, 2**256), A) for _ in range(args.secp_random)]
        cases = [(c, R) for c in (0, 1) for R in special_points(curve, A) + rand]
        t1 = time.time()
        fails = check(gen, ops, cases)
        r = {"mode": "SAMPLED_256", "q": hex(curve.q), "addend": "secp256k1 G" if k == 1 else f"[{hex(k)}]G", "n": gen.n, "wires": C.num_wires,
             "ops": len(ops) // 4, "counts": counts, "lanes": len(cases), "special_cases": ["inf", "A", "-A", "-2A", "2A", "3A"], "failures": len(fails),
             "first_failures": fails[:3], "build_seconds": round(t_build, 1), "sim_seconds": round(time.time() - t1, 1)}
        if sha:
            r.update(k=hex(k), json_sha256=sha)
        print(json.dumps(r), flush=True)
        receipts.append(r)
        del C, ops

    if args.out:
        args.out.write_text(json.dumps(receipts, indent=1))
    sys.exit(1 if any(r["failures"] or not r.get("expansion_matches_flatten", True) for r in receipts) else 0)


if __name__ == "__main__":
    main()
