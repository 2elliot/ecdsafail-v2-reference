"""Full-width component checks on boundary inputs that random curve points
almost never reach: values near q and 2^n, x = 0, and long Euclid runs. The
controlled and constant adders run with ctrl 0 and 1 on the kernel's constants
for G and on sparse and dense edge constants."""
import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.arith import Mod, add_const, gt_const_into
from qec.builder import Circuit
from qec.curve import SECP256K1_G
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
    f = 2**256 - q
    edge = [0, 1, 2, 3, q - 1, q - 2, (q - 1) // 2, (q + 1) // 2, 2**255, 2**255 - 1, f, q - f, 2**256 - 1, q - 2**32, 2**256 - 2**33 - 1954 - q % 7]
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

    def b_cmodadd(C):
        c, a, b = C.alloc(1), C.alloc(n), C.alloc(n)
        mod.modadd(C, a, b, ctrl=c[0])
        return [c, a, b], [c, a, b]
    out["modadd_ctrl"] = run_case(b_cmodadd, [(c, a, b) for c in (0, 1) for a, b in pairs], lambda c, a, b: (c, a, (b + c * a) % q))

    a_G, b_G = SECP256K1_G
    consts = {"q-a_G": q - a_G, "q-b_G": q - b_G, "3a_G": 3 * a_G % q, "a_G": a_G, "1": 1, "f": f, "2^255": 2**255,
              "q-1": q - 1, "q-2": q - 2, "2^255-1": 2**255 - 1, "2^256-1 mod q": (2**256 - 1) % q}
    for name, v in consts.items():
        near = sorted({x % q for x in (q - v - 1, q - v, q - v + 1, v - 1, v, v + 1)} | set(pool))

        def b_cmc(C, v=v):
            c, b = C.alloc(1), C.alloc(n)
            mod.cmodadd_const(C, c[0], v, b)
            return [c, b], [c, b]
        out[f"cmodadd_const_ctrl[{name}]"] = run_case(b_cmc, [(c, b) for c in (0, 1) for b in near], lambda c, b, v=v: (c, (b + c * v) % q))

        near = sorted({x % 2**n for x in (v - 1, v, v + 1, 2**n - v - 1, 2**n - v, 2**n - v + 1, 2**n - 1)} | set(pool))

        def b_gt(C, v=v):
            c, b, t = C.alloc(1), C.alloc(n), C.alloc(1)
            gt_const_into(C, v, b, t[0], ctrl=c[0])
            return [c, b, t], [c, b, t]
        out[f"gt_const_into_ctrl[{name}]"] = run_case(b_gt, [(c, b, t) for c in (0, 1) for b in near for t in (0, 1)],
                                                     lambda c, b, t, v=v: (c, b, t ^ (c & (v > b))))

        def b_gt_u(C, v=v):
            b, t = C.alloc(n), C.alloc(1)
            gt_const_into(C, v, b, t[0])
            return [b, t], [b, t]
        out[f"gt_const_into[{name}]"] = run_case(b_gt_u, [(b, t) for b in near for t in (0, 1)], lambda b, t, v=v: (b, t ^ (v > b)))

        def b_addc(C, v=v):
            c, b, o = C.alloc(1), C.alloc(n), C.alloc(1)
            add_const(C, b, v, ctrl=c[0], cout=o[0])
            return [c, b, o], [c, b, o]
        out[f"add_const_ctrl_cout[{name}]"] = run_case(b_addc, [(c, b, o) for c in (0, 1) for b in near for o in (0, 1)],
                                                      lambda c, b, o, v=v: (c, (b + c * v) % 2**n, o ^ (c & (b + v >= 2**n))))

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


def _bz_pack(values, wires):
    """Byte-backed transpose avoids quadratic big-int growth at exhaustive sizes."""
    masks = [bytearray((len(values) + 7) // 8) for _ in wires]
    for lane, value in enumerate(values):
        byte, bit = lane >> 3, 1 << (lane & 7)
        for i, mask in enumerate(masks):
            if value >> i & 1:
                mask[byte] |= bit
    return {w: int.from_bytes(mask, 'little') for w, mask in zip(wires, masks)}


def _bz_state(C, ins, outs, inputs, expected, ops=None):
    init, want = {}, {}
    for i, reg in enumerate(ins):
        init.update(_bz_pack([v[i] for v in inputs], reg))
    for i, reg in enumerate(outs):
        want.update(_bz_pack([v[i] for v in expected], reg))
    state = sim.run(C.ops if ops is None else ops, C.num_wires, init, len(inputs))
    wrong = dirty = 0
    for wire, value in enumerate(state):
        wrong |= value ^ want.get(wire, 0)
        if wire not in want:
            dirty |= value
    return wrong, dirty, init, want


def _bz_build(q, m=None, inverse=False, step=False):
    C = Circuit()
    n = q.bit_length()
    regs = [C.alloc(2 if step else 1), C.alloc(n), C.alloc(n)]
    ins = [r[:] for r in regs]
    mod = Mod(q)
    if step:
        block = DialogMultiplier(C, mod)._bezout()
    else:
        block = C.record(lambda g, a, s: mod.dbl_add(C, g[0], a, s, m), regs)
    C.apply(block, regs, inverse=inverse)
    return C, ins, regs


def _bz_exhaust(q, m=None, inverse=False, step=False):
    C, ins, outs = _bz_build(q, m, inverse, step)
    batch, wants = [], []
    failures = lanes = 0
    for g in range(4 if step else 2):
        for r in range(q):
            for s in range(q):
                batch.append((g, r, s))
                v = ((s - g * r) * pow(2, -1, q) if inverse else 2 * s + (g & 1) * r) % q
                wants.append((g, v, r) if step and g & 2 else (g, r, v))
                if len(batch) == 8192:
                    wrong, _, _, _ = _bz_state(C, ins, outs, batch, wants)
                    failures += wrong.bit_count(); lanes += len(batch)
                    batch, wants = [], []
    if batch:
        wrong, _, _, _ = _bz_state(C, ins, outs, batch, wants)
        failures += wrong.bit_count(); lanes += len(batch)
    return {'lanes': lanes, 'failures': failures, 'ops': len(C), 'wires': C.num_wires}


def _bz_census(q):
    from array import array
    from qec.builder import CCX, CX
    C, ins, outs = _bz_build(q)
    inputs = [(g, r, s) for g in (0, 1) for r in range(q) for s in range(q)]
    expected = [(g, r, (2 * s + g * r) % q) for g, r, s in inputs]
    wrong, _, init, want = _bz_state(C, ins, outs, inputs, expected)
    assert wrong == 0
    caught, survivors = [], []
    for pos in range(0, len(C.ops), 4):
        if C.ops[pos] != CCX:
            continue
        ops = array('i', C.ops)
        _, a, _, c = ops[pos:pos+4]
        ops[pos:pos+4] = array('i', (CX, a, c, c))
        state = sim.run(ops, C.num_wires, init, len(inputs))
        bad = any(value != want.get(w, 0) for w, value in enumerate(state))
        (caught if bad else survivors).append(pos // 4)
    return {'lanes': len(inputs), 'failures': int(not caught), 'caught': len(caught),
            'survivors': survivors, 'survivor_reason': 'equivalent on every canonical input, exhaustively checked including scratch'}


def _bz_order_control(q):
    from array import array
    from qec.builder import CCX
    C, ins, outs = _bz_build(q)
    inputs = [(g, r, s) for g in (0, 1) for r in range(q) for s in range(q)]
    expected = [(g, r, (2 * s + g * r) % q) for g, r, s in inputs]
    wrong, _, init, want = _bz_state(C, ins, outs, inputs, expected)
    assert wrong == 0
    # The terminal prefix gate is the last compute CCX targeting a scratch
    # prefix and controlled by the top high bit. Identify its inverse by wires.
    high_top = ins[2][-2]
    candidates = []
    live = {w for reg in ins for w in reg}
    ops = [tuple(C.ops[p:p+4]) for p in range(0, len(C.ops), 4)]
    for i, op in enumerate(ops):
        if op[0] == CCX and high_top in op[1:3] and op[3] not in live:
            later = [j for j in range(i+1, len(ops)) if ops[j] == op]
            if later:
                j = later[0]
                if j+1 < len(ops) and ops[j+1][0] == CCX and ops[j+1][3] == high_top:
                    candidates.append(j)
    assert len(candidates) == 1, candidates
    j = candidates[0]
    result = {}
    for name, edited in [('deleted', ops[:j]+ops[j+1:]),
                         ('late', ops[:j]+[ops[j+1],ops[j]]+ops[j+2:])]:
        state = sim.run(array('i', (x for op in edited for x in op)), C.num_wires, init, len(inputs))
        dirty = 0
        for w, v in enumerate(state):
            if w not in want: dirty |= v
        result[name] = dirty.bit_count()
    return {'lanes': len(inputs), 'failures': int(any(v == 0 for v in result.values())), 'dirty': result}


def bezout_checks():
    # Missing arithmetic, a wrong reduction band, a dirty carry or wrong
    # output rotation must fail these independent integer-map checks.
    assert hasattr(Mod, 'dbl_add'), 'missing exact Mod.dbl_add'
    out = {}
    primes = (3, 5, 7, 13, 61, 127, 251, 509, 1021, 163, 349, 547, 733)
    for q in primes:
        out[f'cell_exhaustive[{q}]'] = _bz_exhaust(q)
        out[f'cell_inverse_exhaustive[{q}]'] = _bz_exhaust(q, inverse=True)
    for q in (31, 61, 127, 251, 509):
        parts = [_bz_exhaust(q, m) for m in range(Mod(q).low_bits, q.bit_length()+1)]
        out[f'cell_every_m[{q}]'] = {'lanes': sum(x['lanes'] for x in parts), 'failures': sum(x['failures'] for x in parts), 'parts': parts}
    C, _, _ = _bz_build(2**256-2**32-977)
    ccx = sum(op == 2 for op in C.ops[::4])
    out['cell_count'] = {'lanes': 1, 'failures': int(ccx != 2086), 'CCX': ccx}
    for q in (127, 509):
        out[f'cell_census[{q}]'] = _bz_census(q)
        out[f'chain_order_control[{q}]'] = _bz_order_control(q)
    print(json.dumps(out, indent=1))
    return int(any(row['failures'] for row in out.values()))


if __name__ == "__main__":
    if "--bezout" in sys.argv:
        sys.exit(bezout_checks())
    main()
