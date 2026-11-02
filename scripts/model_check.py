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
from qec.model import euclid, inplace_mul, iteration_bound, offset_pass, sweep_tape, width_schedule
from qec.point_add import DialogMultiplier


def gate_mul(q: int, ys, inverse: bool = False, iters: int | None = None) -> int:
    """Check every x against the supplied y values."""
    return gate_pairs(q, [(x, y) for x in range(q) for y in ys], inverse, iters)


def gate_pairs(q: int, pairs, inverse: bool = False, iters: int | None = None) -> int:
    """Count bad output or dirty-scratch lanes; malformed builds raise."""
    pairs = list(pairs)
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


def _pass_bad(q, xs, iters):
    C = Circuit()
    mod = Mod(q)
    m = DialogMultiplier(C, mod)
    if iters is not None:
        m.iters = iters
    U, V = C.alloc(mod.n), C.alloc(mod.n)
    inputs = list(V)
    G = C.alloc(2 * m.iters)
    load_const(C, U, q)
    m._euclid_pass(U, V, G)
    assert len(set(U + V + G)) == len(U + V + G), "aliased pass register"
    structural(C.prims, C.num_wires)
    state = sim.run(C.prims, C.num_wires, sim.pack(xs, inputs), len(xs))
    # Keep the reference model independent of offset_pass and sweep_tape.
    tapes = [sum((odd + 2 * swap) << (2 * i)
                 for i, (odd, swap) in enumerate(euclid(q, x, m.iters)[2])) for x in xs]
    want = sim.pack([q if x == 0 else 1 for x in xs], U) | sim.pack(tapes, G)
    bad = 0
    for w, actual in enumerate(state):
        bad |= actual ^ want.get(w, 0)
    return bad.bit_count()


def gate_pass(q: int, iters: int | None = None) -> int:
    """Check terminal words, exact tape and all scratch on remapped inputs."""
    return _pass_bad(q, [x or 1 for x in range(q)], iters)


def gate_pass_zero(q: int) -> int:
    """A generic zero caller keeps q and writes zero tape and scratch."""
    return _pass_bad(q, [0], None)


def width_edge_inputs(q: int) -> list[tuple[int, int]]:
    """Schedule-edge powers and neighbours, boundaries, and seed-2026 samples."""
    rng = random.Random(2026)
    xs = {x for k in range(q.bit_length())
          for x in (2**k, 2**k - 1, 2**k + 1) if 0 <= x < q}
    xs |= {0, 1, 3, q - 1, q - 2, (q - 1) // 2, (q + 1) // 2}
    xs |= {rng.randrange(q) for _ in range(64)}
    return [(x, rng.choice([0, 1, q - 1, rng.randrange(q)])) for x in sorted(xs)]


def offset_checks(q, xs):
    """Compare held states and swept cells with the independent reference dialog."""
    T = iteration_bound(q)
    offset_ok = tape_ok = True
    for x in xs:
        states, cells = offset_pass(q, x, T)
        if len(states) != T + 1 or len(cells) != T:
            return False, False
        u, v = q, x
        reference = euclid(q, x, T)[2]
        for i, (U, V, e) in enumerate(states):
            offset_ok &= U == u and V + e == v and V >= 0
            offset_ok &= e == (cells[i - 1][1] if i else 0)
            if i < T:
                odd, swap = reference[i]
                if swap:
                    u, v = v, u
                v = (v - odd * u) // 2
        tape_ok &= sweep_tape(cells) == reference and states[-1][1:] == (0, 0)
    return offset_ok, tape_ok


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
    pass_zero = gate_pass_zero(q) == 0
    forward_bad = gate_mul(q, ys)
    inverse_bad = gate_mul(q, ys, inverse=True)
    pass_offset, pass_tape = offset_checks(q, range(1, q))
    rows.append({"q": q, "worst_iterations": worst, "bound": B,
                 "width_schedule_ok": width_ok, "multiplier_ok": mul_ok,
                 "pass_end_ok": pass_end, "remap_control_ok": remap_control,
                 "pass_short_rejected": short_rejected, "pass_gates_ok": pass_gates,
                 "multiplier_gates_ok": forward_bad == inverse_bad == 0,
                 "forward_bad": forward_bad, "inverse_bad": inverse_bad,
                 "pass_offset_ok": pass_offset, "pass_tape_ok": pass_tape, "pass_zero_ok": pass_zero})
    print(json.dumps(rows[-1]))
q = 2**256 - 2**32 - 977
rng = random.Random(1)
pairs = [(rng.randrange(q), rng.randrange(q)) for _ in range(2000)]
ok = all(inplace_mul(q, x, y) == (x * y % q if x else y) for x, y in pairs)
pass_end = all(euclid(q, x or 1, iteration_bound(q))[:2] == (1, 0)
               for x in [0, 1, q - 1] + [x for x, _ in pairs])
offset_ok = all(offset_checks(q, [1, q - 1] + [x for x, _ in pairs]))
edge_pairs = width_edge_inputs(q)
width_edge_ok = gate_pairs(q, edge_pairs) == 0 and gate_pairs(q, edge_pairs, inverse=True) == 0
print(json.dumps({"q": "secp256k1", "bound": iteration_bound(q), "random_pairs": 2000,
                  "multiplier_ok": ok, "pass_end_ok": pass_end, "offset_ok": offset_ok,
                  "width_edge_ok": width_edge_ok, "width_edge_lanes": len(edge_pairs)}))
keys = ("width_schedule_ok", "multiplier_ok", "pass_end_ok", "remap_control_ok",
        "pass_short_rejected", "pass_gates_ok", "multiplier_gates_ok",
        "pass_offset_ok", "pass_tape_ok", "pass_zero_ok")
sys.exit(0 if ok and pass_end and offset_ok and width_edge_ok and all(r["worst_iterations"] <= r["bound"]
         and all(r[k] for k in keys) for r in rows) else 1)
