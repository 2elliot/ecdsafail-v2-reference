"""Exhaustive check of SquareFold.compute under the module's current flags.

    python3.12 -I -B scripts/square_check.py <tree> <nmax> [broken | brokenfold | brokenmodadd] [MODULE.FLAG=value ...]

For every odd q < 2^nmax that SquareFold applies to (f = 2^n - q, f*f <= q) and every n-bit Y:
J[:n] must end Y^2 mod q, Y must be restored, and every wire outside Y and J must end 0. Also
prefix_xor and fan_into on every input up to nmax bits. SquareFold.modadd on every canonical a, b and
control for those q below 2^8, and cadd_lookahead and gt_lookahead on every input up to 6 bits. With "broken", prefix_xor drops its last CX;
with "brokenfold", the fold's first add drops its carry in; with "brokenmodadd", cadd_lookahead's
second X layer is uncontrolled. Each must fail.
"""
import ast, importlib, sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
from qec import arith, sim, square
from qec.builder import Circuit

NMAX = int(sys.argv[2])
args = sys.argv[3:]
if args and args[0] == "broken":
    args = args[1:]
    orig = arith.prefix_xor

    def broken(C, a):
        g = arith._Gates()
        orig(g, a)
        for o in g.ops[:-1]:
            C.cx(*o)
    arith.prefix_xor = broken
    square.prefix_xor = broken
if args and args[0] == "brokenfold":
    args = args[1:]
    orig_add = square._add

    def broken_add(C, a, b, cin=None, temps=()):
        orig_add(C, a, b, None, temps)
    square._add = broken_add
if args and args[0] == "brokenmodadd":
    args = args[1:]
    orig_emit = arith._emit

    def cadd_broken(C, a, b, cout, ctrl, temps, copies):
        (cin,) = C.alloc(1)
        up = arith._up_sweep(a, b, cin, temps)
        orig_emit(C, up.ops)
        if cout is not None:
            C.ccx(ctrl, a[-1], cout)
        C.ccx(ctrl, cin, b[0])
        for i in range(1, len(a)):
            C.ccx(ctrl, a[i - 1], b[i])
        for x in b:
            C.cx(ctrl, x)
        orig_emit(C, reversed(up.ops))
        for x in b:
            C.x(x)
        C.free([cin])
    arith.cadd_lookahead = cadd_broken
    square.cadd_lookahead = cadd_broken
for kv in args:
    k, v = kv.split("=", 1); mod, flag = k.rsplit(".", 1)
    setattr(importlib.import_module("qec." + mod), flag, ast.literal_eval(v))


def lanes_of(n):
    return {i: sum(1 << y for y in range(1 << n) if y >> i & 1) for i in range(n)}


bad = cases = 0
for n in range(1, NMAX + 1):
    C = Circuit()
    a = C.alloc(n)
    arith.prefix_xor(C, a)
    m = lanes_of(n)
    st = sim.run(C.ops, C.num_wires, {w: m[i] for i, w in enumerate(a)}, 1 << n)
    for y in range(1 << n):
        got = sum((st[w] >> y & 1) << i for i, w in enumerate(a))
        want = sum((bin(y & ((2 << i) - 1)).count("1") & 1) << i for i in range(n))
        bad += got != want
    for k in range(4):
        C = Circuit()
        src, T, cp = C.alloc(1), C.alloc(n), C.alloc(k)
        arith.fan_into(C, src[0], T, cp)
        ins = src + T
        m = lanes_of(n + 1)
        st = sim.run(C.ops, C.num_wires, {w: m[i] for i, w in enumerate(ins)}, 1 << (n + 1))
        for v in range(1 << (n + 1)):
            s = v & 1
            got = [st[w] >> v & 1 for w in ins]
            want = [s] + [(v >> (i + 1) & 1) ^ s for i in range(n)]
            bad += got != want
        bad += any(st[w] for w in cp)
    cases += 1 << n

for q in range(3, 1 << NMAX, 2):
    if not square.SquareFold.applies(q):
        continue
    mod = arith.Mod(q)
    fold = square.SquareFold(mod)
    n = mod.n
    C = Circuit()
    Y, J = C.alloc(n), C.alloc(fold.junk)
    fold.compute(C, Y, J)
    m = lanes_of(n)
    init = {w: m[i] for i, w in enumerate(Y)}
    st = sim.run(C.ops, C.num_wires, init, 1 << n)
    for y in range(1 << n):
        got = sum((st[w] >> y & 1) << i for i, w in enumerate(J[:n]))
        bad += got != y * y % q
    bad += any(st[w] != init[w] for w in Y)
    live = set(Y) | set(J)
    bad += any(st[w] for w in range(C.num_wires) if w not in live)
    cases += 1 << n

    if n > 8:
        continue
    C = Circuit()
    (k,), a, b = C.alloc(1), C.alloc(n), C.alloc(n)
    fold.modadd(C, a, b, k)
    vals = [(c, x, y) for c in (0, 1) for x in range(q) for y in range(q)]
    cs, xs, ys = zip(*vals)
    init = sim.pack(cs, [k]) | sim.pack(xs, a) | sim.pack(ys, b)
    st = sim.run(C.ops, C.num_wires, init, len(vals))
    want = sim.pack(cs, [k]) | sim.pack(xs, a) | sim.pack([(y + c * x) % q for c, x, y in vals], b)
    bad += sum(st[w] != want.get(w, 0) for w in range(C.num_wires))
    cases += len(vals)

for w in range(1, 7):
    for nt, ncp in ((1, 0), (2, 1), (w, 3), (2 * w, 0)):
        for kind in ("cadd", "gt"):
            C = Circuit()
            (k,), a, b, (o,), T, cp = C.alloc(1), C.alloc(w), C.alloc(w), C.alloc(1), C.alloc(nt), C.alloc(ncp)
            if kind == "cadd":
                arith.cadd_lookahead(C, a, b, o, k, T, cp)
            else:
                arith.gt_lookahead(C, a, b, o, k, T)
            ins = [k] + a + b + [o]
            m = lanes_of(len(ins))
            st = sim.run(C.ops, C.num_wires, {x: m[i] for i, x in enumerate(ins)}, 1 << len(ins))
            for v in range(1 << len(ins)):
                c, x, y, ov = v & 1, v >> 1 & ((1 << w) - 1), v >> (w + 1) & ((1 << w) - 1), v >> (2 * w + 1) & 1
                got = [st[z] >> v & 1 for z in ins]
                if kind == "cadd":
                    s = y + c * x
                    want_y, want_o = s % (1 << w), ov ^ (s >> w)
                else:
                    want_y, want_o = y, ov ^ (c & (x > y))
                want = [c] + [x >> i & 1 for i in range(w)] + [want_y >> i & 1 for i in range(w)] + [want_o]
                bad += got != want
            bad += any(st[z] for z in T + cp)
            cases += 1 << len(ins)
print("cases", cases, "bad", bad)
print("OK" if bad == 0 else "FAIL")
sys.exit(0 if bad == 0 else 1)
