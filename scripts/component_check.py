"""Full-width component checks on boundary inputs that random curve points
almost never reach: values near q and 2^n, x = 0, and long Euclid runs."""
import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.arith import Mod
from qec.builder import Circuit
from qec.model import iteration_bound
from qec.point_add import DialogMultiplier


def euclid_iters(q, x):
    u, v, k = q, x, 0
    while v:
        if v & 1 and u > v:
            u, v = v, u
        if v & 1:
            v -= u
        v >>= 1
        k += 1
    return k


def run_case(build, inputs, expect):
    """build(C) -> (input regs, output regs). inputs: list of tuples."""
    C = Circuit()
    ins, outs = build(C)
    lanes = len(inputs)
    init = {}
    for k, reg in enumerate(ins):
        init.update(sim.pack([v[k] for v in inputs], reg))
    state = sim.run(C.ops, C.num_wires, init, lanes)
    got = [sim.unpack(state, reg, lanes) for reg in outs]
    live = {w for r in outs for w in r}
    dirty = 0
    for w in range(C.num_wires):
        if w not in live:
            dirty |= state[w]
    fails = [i for i, v in enumerate(inputs) if tuple(g[i] for g in got) != expect(*v) or dirty >> i & 1]
    return {"ops": len(C), "wires": C.num_wires, "lanes": lanes, "failures": len(fails), "first": [[hex(x) for x in inputs[i]] for i in fails[:2]]}


def main():
    q = 2**256 - 2**32 - 977
    n = q.bit_length()
    mod = Mod(q)
    rng = random.Random(5)
    edge = [0, 1, 2, 3, q - 1, q - 2, (q - 1) // 2, (q + 1) // 2, 2**255, 2**255 - 1, 2**32 + 977, q - 2**32, 2**256 - 2**33 - 1954 - q % 7]
    edge = [e % q for e in edge]
    pool = edge + [rng.randrange(q) for _ in range(40)]
    cands = [rng.randrange(1, q) for _ in range(20000)] + [1, 2, q - 1] + [2**k for k in range(256)]
    long_x = sorted(cands, key=lambda x: -euclid_iters(q, x))[:24]
    out = {"q": hex(q), "iteration_bound": iteration_bound(q), "longest_euclid_run_in_inputs": euclid_iters(q, long_x[0])}

    def b_modadd(C):
        a, b = C.alloc(n), C.alloc(n)
        mod.modadd(C, a, b)
        return [a, b], [a, b]
    pairs = [(a, b) for a in edge for b in edge] + [(rng.randrange(q), rng.randrange(q)) for _ in range(60)]
    out["modadd"] = run_case(b_modadd, pairs, lambda a, b: (a, (a + b) % q))

    def b_dbl(C):
        s = C.alloc(n)
        s_in = list(s)
        mod.moddbl(C, s)
        return [s_in], [s]
    out["moddbl"] = run_case(b_dbl, [(v,) for v in pool], lambda s: (2 * s % q,))

    def b_neg(C):
        (c,) = C.alloc(1)
        s = C.alloc(n)
        mod.cneg(C, c, s)
        return [[c], s], [[c], s]
    out["cneg"] = run_case(b_neg, [(c, v) for c in (0, 1) for v in pool if v], lambda c, s: (c, (q - s) % q if c else s))

    def b_mul(C):
        x, y = C.alloc(n), C.alloc(n)
        x_in, y_in = list(x), list(y)
        DialogMultiplier(C, mod).apply(x, y)
        return [x_in, y_in], [x, y]
    xs = [0, 1, 2, q - 1, q - 2, 2**255] + long_x + [rng.randrange(1, q) for _ in range(16)]
    ys = [0, 1, q - 1, 2**255 % q, rng.randrange(q)]
    mul_in = [(x, ys[i % len(ys)]) for i, x in enumerate(xs)] + [(0, y) for y in ys] + [(x, q - 1) for x in long_x[:6]]
    t = time.time()
    out["dialog_mul"] = run_case(b_mul, mul_in, lambda x, y: (x, x * y % q if x else y))
    out["dialog_mul"]["seconds"] = round(time.time() - t, 1)
    print(json.dumps(out, indent=1))
    sys.exit(1 if any(isinstance(v, dict) and v["failures"] for v in out.values()) else 0)


if __name__ == "__main__":
    main()
