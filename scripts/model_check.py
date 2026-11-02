"""Classical and bit-sliced gate checks of the dialog multiplier and pass.
Small primes cover every x; y is exhaustive below 300 and sampled above."""
import io
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

sys.path.insert(0, str(Path(__file__).resolve().parent))

from check import expand, structural
from qec import flatten, sim
from qec.arith import Mod, load_const
from qec.builder import Circuit
from qec.emit import write_circuit
from qec.model import euclid, inplace_mul, iteration_bound, width_schedule
from qec.point_add import DialogMultiplier


def gate_mul(q: int, ys, inverse: bool = False, iters: int | None = None) -> int:
    """Count bad output or dirty-scratch lanes; malformed builds raise."""
    C = Circuit(hierarchical=True)
    mod = Mod(q)
    X, Y = C.alloc(mod.n), C.alloc(mod.n)
    inputs = (list(X), list(Y))
    m = DialogMultiplier(C, mod)
    if iters is not None:
        m.iters = iters
    m.apply(X, Y, inverse=inverse)
    assert len(set(X + Y)) == 2 * mod.n, "aliased output register"
    buf = io.StringIO()
    write_circuit(buf, C, {"control": [], "x": X, "y": Y})
    ops = expand(flatten.load(io.StringIO(buf.getvalue())))
    structural(ops, C.num_wires)
    pairs = [(x, y) for x in range(q) for y in ys]
    xs, values = zip(*pairs)
    init = sim.pack(xs, inputs[0]) | sim.pack(values, inputs[1])
    state = sim.run(ops, C.num_wires, init, len(pairs))
    expected = [y if x == 0 else y * (pow(x, -1, q) if inverse else x) % q
                for x, y in pairs]
    want = sim.pack(xs, X) | sim.pack(expected, Y)
    bad = 0
    for w, actual in enumerate(state):
        bad |= actual ^ want.get(w, 0)
    return bad.bit_count()


def gate_pass(q: int, iters: int | None = None) -> int:
    """Count remapped input lanes whose gate pass does not end at (1, 0)."""
    C = Circuit()
    mod = Mod(q)
    m = DialogMultiplier(C, mod)
    if iters is not None:
        m.iters = iters
    U, V = C.alloc(mod.n), C.alloc(mod.n)
    inputs = list(V)
    G = C.alloc(2 * m.iters)
    load_const(C, U, q)
    for i in range(m.iters):
        C.apply(m._euclid_step(m.widths[i]), [U, V, G[2 * i:2 * i + 2]])
    structural(C.prims, C.num_wires)
    state = sim.run(C.prims, C.num_wires, sim.pack([x or 1 for x in range(q)], inputs), q)
    bad = state[U[0]] ^ ((1 << q) - 1)
    for w in U[1:] + V:
        bad |= state[w]
    return bad.bit_count()


def trace(q, x):
    u, v, sizes = q, x, []
    while v:
        sizes.append(max(u, v).bit_length())
        if v & 1 and u > v:
            u, v = v, u
        if v & 1:
            v -= u
        v >>= 1
    return sizes


rows = []
for q in [3, 5, 7, 11, 13, 31, 43, 61, 127, 251, 509, 1021, 4093]:
    B, W = iteration_bound(q), width_schedule(q, iteration_bound(q))
    worst = 0
    width_ok = True
    for x in range(q):
        t = trace(q, x)
        worst = max(worst, len(t))
        width_ok &= all(s <= W[i] for i, s in enumerate(t))
    ys = range(q) if q < 300 else random.Random(q).sample(range(q), 8)
    mul_ok = all(inplace_mul(q, x, y) == (x * y % q if x else y) for x in range(q) for y in ys)
    pass_end = all(euclid(q, x or 1, B)[:2] == (1, 0) for x in range(q))
    remap_control = [x for x in range(q) if euclid(q, x, B)[:2] != (1, 0)] == [0]
    short_rejected = any(euclid(q, x or 1, worst - 1)[1] != 0 for x in range(q))
    short_rejected &= gate_mul(q, ys, iters=worst - 1) > 0
    pass_gates = gate_pass(q) == 0
    forward_bad = gate_mul(q, ys)
    inverse_bad = gate_mul(q, ys, inverse=True)
    rows.append({"q": q, "worst_iterations": worst, "bound": B,
                 "width_schedule_ok": width_ok, "multiplier_ok": mul_ok,
                 "pass_end_ok": pass_end, "remap_control_ok": remap_control,
                 "pass_short_rejected": short_rejected, "pass_gates_ok": pass_gates,
                 "multiplier_gates_ok": forward_bad == inverse_bad == 0,
                 "forward_bad": forward_bad, "inverse_bad": inverse_bad})
    print(json.dumps(rows[-1]))
q = 2**256 - 2**32 - 977
rng = random.Random(1)
pairs = [(rng.randrange(q), rng.randrange(q)) for _ in range(2000)]
ok = all(inplace_mul(q, x, y) == (x * y % q if x else y) for x, y in pairs)
pass_end = all(euclid(q, x or 1, iteration_bound(q))[:2] == (1, 0)
               for x in [0, 1, q - 1] + [x for x, _ in pairs])
print(json.dumps({"q": "secp256k1", "bound": iteration_bound(q), "random_pairs": 2000,
                  "multiplier_ok": ok, "pass_end_ok": pass_end}))
keys = ("width_schedule_ok", "multiplier_ok", "pass_end_ok", "remap_control_ok",
        "pass_short_rejected", "pass_gates_ok", "multiplier_gates_ok")
sys.exit(0 if ok and pass_end and all(r["worst_iterations"] <= r["bound"]
         and all(r[k] for k in keys) for r in rows) else 1)
