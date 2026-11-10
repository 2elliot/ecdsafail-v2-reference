"""Phase names for the ops of chosen definitions of an emitted (and post-passed) circuit.

    python3.12 -I -B phasemap.py <tree> <circuit.json> <out_phases.json> <def_id> [def_id ...] [-- MODULE.FLAG=value ...]

Builds the tree's secp256k1 circuit with wrappers on the arithmetic helpers that tag the op range each
call emits, then matches each op of the circuit's definition (post-passed: reordered, some pairs
cancelled) to the k-th raw op with the same content. Each op gets the innermost tag path around its raw
op, or "body@k" (k = top-level tags before it) outside every tag.
"""
import ast, importlib, json, sys
from collections import defaultdict, deque
from pathlib import Path
argv = sys.argv[1:]
flags = argv[argv.index("--") + 1:] if "--" in argv else []
argv = argv[:argv.index("--")] if "--" in argv else argv
tree, circ, outp, ids = Path(argv[0]).resolve(), argv[1], argv[2], [int(x) for x in argv[3:]]
sys.path.insert(0, str(tree))
for kv in flags:
    k, v = kv.split("=", 1)
    mod, flag = k.rsplit(".", 1)
    setattr(importlib.import_module("qec." + mod), flag, ast.literal_eval(v))
from qec import builder, point_add, arith, square, two_bit
from qec.emit import definitions
from qec.curve import SECP256K1, SECP256K1_G

TAGS = []  # (abs_start, abs_end, path)
STACK = []
COUNT = defaultdict(int)

def wrap(mod, name, method=False):
    orig = getattr(mod, name)
    def f(*a, **k):
        C = a[1] if method else (a[0] if a else None)
        if not isinstance(C, builder.Circuit):
            return orig(*a, **k)
        parent = "/".join(STACK)
        COUNT[(parent, name)] += 1
        STACK.append(f"{name}#{COUNT[(parent, name)]}")
        s = len(C.prims) // 4
        try:
            return orig(*a, **k)
        finally:
            TAGS.append((s, len(C.prims) // 4, "/".join(STACK)))
            STACK.pop()
    setattr(mod, name, f)

for name in ("shared_difference_u6", "shared_difference_u1", "shared_difference_ripple", "_fan", "_unfan", "_emit",
             "_par_rotate", "_fanout", "_unfanout", "cswap", "adjacent_xor", "lookahead_sweep", "add_const"):
    if hasattr(point_add, name):
        wrap(point_add, name)
for name in ("add_lookahead", "add", "sub", "gt_into", "_carries", "_uncarry", "add_const"):
    wrap(square, name)
for name in ("lookahead_chain", "_emit"):
    wrap(arith, name)
for name in ("_rows", "_multiple"):
    wrap(square.SquareFold, name, method=True)
for name in ("eq_const_into", "load_const", "_restore_layout"):
    wrap(point_add, name)
for name in ("cmodadd_const", "creflect", "modadd", "cneg"):
    wrap(arith.Mod, name, method=True)
import functools
for name in ("_eq_flags", "_kernel"):
    orig_m = getattr(point_add.ControlledPointAdd, name)
    def mk(orig_m, name):
        def f(self, *a, **k):
            C = self.C
            parent = "/".join(STACK)
            COUNT[(parent, name)] += 1
            STACK.append(f"{name}#{COUNT[(parent, name)]}")
            s0 = len(C.prims) // 4
            try:
                return orig_m(self, *a, **k)
            finally:
                TAGS.append((s0, len(C.prims) // 4, "/".join(STACK)))
                STACK.pop()
        return f
    setattr(point_add.ControlledPointAdd, name, mk(orig_m, name))

orig_record = builder.Circuit.record
def record(self, fn, regs):
    start = len(self.prims) // 4
    saved = list(STACK); STACK.clear()  # tags inside a block are relative to it
    COUNT.clear()
    b = orig_record(self, fn, regs)
    STACK[:] = saved
    mine = [t for t in TAGS if t[0] >= start]
    TAGS[:] = [t for t in TAGS if t[0] < start]
    b.tags = [(s - start, e - start, p) for s, e, p in mine]
    b.label = getattr(fn, "__qualname__", "?")
    return b
builder.Circuit.record = record

gen = point_add.ControlledPointAdd(SECP256K1, SECP256K1_G)
C = gen.build(True)
order = definitions(C)
doc = json.load(open(circ)); defs = doc.get("circuit", doc)["definitions"]
NAMES = {builder.X: "X", builder.CX: "CX", builder.CCX: "CCX", builder.SWAP: "SWAP"}
def key_raw(op, a, b, c):
    o = NAMES[op]
    if o == "X": return ("X", a)
    if o == "CX": return ("CX", a, b)
    if o == "SWAP": return ("SWAP", min(a, b), max(a, b))
    return ("CCX", min(a, b), max(a, b), c)
def key_json(r):
    o = r["op"]
    if o == "X": return ("X", r["target"])
    if o == "CX": return ("CX", r["control"], r["target"])
    if o == "SWAP": return ("SWAP", min(r["wire1"], r["wire2"]), max(r["wire1"], r["wire2"]))
    return ("CCX", min(r["control1"], r["control2"]), max(r["control1"], r["control2"]), r["target"])
out = {}
TOPBLK = builder.Block(C.prims, [], [], C.calls)
TOPBLK.tags = list(TAGS); TOPBLK.label = "top"
for d in ids:
    blk = order[d] if d >= 0 else TOPBLK
    if d >= 0:
        assert not blk.calls and defs[d]["arity"] == len(blk.wires), (d, "has calls or pruned wires")
        local = {w: k for k, w in enumerate(sorted(blk.wires))}
    else:
        local = {w: w for w in blk.wires}
    n = len(blk.ops) // 4
    phase = ["body@0"] * n
    depth = [-1] * n
    tops = sorted((s, e) for s, e, p in blk.tags if "/" not in p)
    for s, e, p in blk.tags:
        dep = p.count("/")
        for i in range(s, e):
            if dep > depth[i]:
                depth[i], phase[i] = dep, p
    import bisect
    starts = [s for s, _ in tops]
    for i in range(n):
        if depth[i] < 0:
            phase[i] = f"body@{bisect.bisect_right(starts, i)}"
    q = defaultdict(deque)
    for i in range(n):
        op, a, b, c = blk.ops[4 * i:4 * i + 4]
        q[key_raw(op, local[a], local[b], local[c])].append(i)
    names = []
    for r in (defs[d]["operations"] if d >= 0 else doc.get("circuit", doc)["operations"]):
        if "call" in r:
            names.append("call"); continue
        dq = q[key_json(r)]
        names.append(phase[dq.popleft()] if dq else "unmatched")
    out[d] = names
    print(d, blk.label, "ops", len(names), "unmatched", names.count("unmatched"))
json.dump(out, open(outp, "w"))
