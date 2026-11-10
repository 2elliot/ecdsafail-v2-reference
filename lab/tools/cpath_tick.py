"""Critical path of a window of Bezout ticks, attributed to the code that emitted each op.

    python3.12 -I -B cpath.py <tree> [MODULE.FLAG=value ...] [--window=a:b] [--inv] [--depth]

Ops are tagged while each tick definition is recorded (the innermost frame in qec/two_bit.py or
qec/arith.py below TwoBitMod.tick, as function:line). The depth pass reorders ops; tags follow ops by
content in order. Each op on the critical path is charged its finish time minus its binding input's.
"""
import ast
import collections
import importlib
import sys
from pathlib import Path

tree = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(tree))
window, inv, depth, detail = (150, 160), False, False, False
for kv in sys.argv[2:]:
    if kv.startswith("--window="):
        a, b = kv.split("=", 1)[1].split(":")
        window = (int(a), int(b))
        continue
    if kv == "--inv":
        inv = True
        continue
    if kv == "--depth":
        depth = True
        continue
    if kv == "--segments":
        continue
    if kv == "--detail":
        detail = True
        continue
    k, v = kv.split("=", 1)
    mod, flag = k.rsplit(".", 1)
    m = importlib.import_module("qec." + mod)
    assert hasattr(m, flag), k
    setattr(m, flag, ast.literal_eval(v))

from qec import builder, depth_pass, two_bit  # noqa: E402
from qec.builder import CCX, CX, SWAP, X  # noqa: E402
from qec.curve import SECP256K1, SECP256K1_G  # noqa: E402
from qec.point_add import ControlledPointAdd  # noqa: E402

TAGS = {}          # id(block) -> list of tags in op order (raw)
cur = {"on": False, "tags": None}
SRC = {str(Path(two_bit.__file__).resolve()), str(Path(importlib.import_module("qec.arith").__file__).resolve())}


def tag_now():
    f = sys._getframe(2)
    chain = []
    while f is not None:
        co = f.f_code
        if co.co_name in ("tick", "_tick") and co.co_filename.endswith("two_bit.py") and not (co.co_name == "tick" and f.f_back and f.f_back.f_code.co_name == "_tick"):
            break
        if co.co_filename in SRC or co.co_filename.endswith("two_bit.py"):
            chain.append(f"{co.co_qualname}:{f.f_lineno}")
        f = f.f_back
    if f is None:
        return "?"
    # outermost function below tick, plus tick's own line
    top = chain[-1] if chain else "tick:0"
    return f"{top.split(':')[0]}@{f.f_lineno}#" + ">".join(reversed(chain[-3:]))


orig = {}
for name in ("x", "cx", "ccx", "swap"):
    orig[name] = getattr(builder.Circuit, name)

    def make(name):
        o = orig[name]

        def g(self, *a):
            n0 = len(self.prims)
            o(self, *a)
            if cur["on"] and len(self.prims) > n0:
                cur["tags"].append(tag_now())
        return g
    setattr(builder.Circuit, name, make(name))

orig_tick = two_bit.TwoBitMod.tick


def tick(self, C, *a, **k):
    start = len(C.prims)
    cur["on"], cur["tags"] = True, []
    try:
        return orig_tick(self, C, *a, **k)
    finally:
        cur["on"] = False
        cur["start"] = start
        cur["last"] = list(cur["tags"])


two_bit.TwoBitMod.tick = tick
orig_record = builder.Circuit.record


def record(self, fn, regs):
    blk = orig_record(self, fn, regs)
    if cur.get("last") is not None and len(cur["last"]) == len(blk.ops) // 4:
        TAGS[id(blk)] = cur["last"]
    cur["last"] = None
    return blk


builder.Circuit.record = record

gen = ControlledPointAdd(SECP256K1, SECP256K1_G)
C = gen.build(True)
M = gen.M
body = M.block
ticks_set = {id(b) for b in M._bezout_ticks.values()}
assert all(t in TAGS for t in ticks_set), "untagged tick definition"
idx = [i for i, c in enumerate(body.calls) if id(c.block) in ticks_set]
region = body.calls[idx[window[0]]:idx[min(window[1], len(idx)) - 1] + 1]

order, seen = [], set()


def visit(calls):
    for c in calls:
        if id(c.block) not in seen:
            seen.add(id(c.block))
            visit(c.block.calls)
            order.append(c.block)


visit(region)
ids = {b: i for i, b in enumerate(order)}
if depth:
    bodies = depth_pass.optimize(order, ids)
else:
    bodies = {b: (b.ops, b.calls) for b in order}


def retag(block, newops):
    """Tags for the reordered body by op content (controls unordered)."""
    old, tags = block.ops, TAGS[id(block)]
    q = collections.defaultdict(collections.deque)
    for i in range(len(old) // 4):
        k, a, b, c = old[4 * i:4 * i + 4]
        key = (k, min(a, b), max(a, b), c) if k == CCX else (k, a, b, c)
        q[key].append(tags[i])
    out = []
    for i in range(len(newops) // 4):
        k, a, b, c = newops[4 * i:4 * i + 4]
        key = (k, min(a, b), max(a, b), c) if k == CCX else (k, a, b, c)
        out.append(q[key].popleft())
    return out


ops, tags = [], []


def expand(block, table, inverse, tick_no):
    bops, calls = bodies[block]
    btags = retag(block, bops) if id(block) in TAGS else None
    m = dict(zip(sorted(block.wires), table))
    seq, pos = [], 0
    for c in calls:
        seq.append(("p", pos, c.pos))
        seq.append(("c", c))
        pos = c.pos
    seq.append(("p", pos, len(bops) // 4))
    for item in (reversed(seq) if inverse else seq):
        if item[0] == "p":
            _, a, b = item
            for i in (range(b - 1, a - 1, -1) if inverse else range(a, b)):
                k, x, y, z = bops[4 * i:4 * i + 4]
                ops.append((k, m[x], m[y], m[z]))
                tags.append(btags[i] if btags else "other")
        else:
            c = item[1]
            expand(c.block, [m[w] for w in c.wires], inverse != c.inverse, tick_no)


pos = region[0].pos
for n, c in enumerate(region):
    for i in range(pos, c.pos):
        ops.append(tuple(body.ops[4 * i:4 * i + 4]))
        tags.append("between")
    expand(c.block, c.wires, c.inverse, n)
    pos = c.pos
if inv:
    ops.reverse()
    tags.reverse()

r, last = {}, {}
fin, bind = [0] * len(ops), [-1] * len(ops)
for i, (k, a, b, c) in enumerate(ops):
    if k == X:
        src = a
        t = r.get(a, 0) + 1
        outs = ((a, t),)
    elif k in (CX, SWAP):
        src = a if r.get(a, 0) >= r.get(b, 0) else b
        t = max(r.get(a, 0), r.get(b, 0)) + (1 if k == CX else 3)
        outs = ((a, t), (b, t))
    else:
        ra, rb, rt = r.get(a, 0), r.get(b, 0), r.get(c, 0)
        s2 = max(rb, rt + 1) + 1
        s4 = max(ra, s2 + 1) + 1
        src = a if ra >= s2 + 1 else (b if rb >= rt + 1 else c)
        t = s4 + 7
        outs = ((c, s4 + 6), (a, s4 + 7), (b, s4 + 7))
    bind[i] = last.get(src, -1)
    fin[i] = t
    for w, tt in outs:
        r[w] = tt
        last[w] = i
end = max(range(len(ops)), key=lambda i: fin[i])
total = fin[end]
charge = collections.Counter()
count = collections.Counter()
i = end
while i >= 0:
    j = bind[i]
    d = fin[i] - (fin[j] if j >= 0 else 0)
    tg = tags[i].split("#")
    key = tg[-1] if detail else (tg[0].split("@")[0] if not tg[0].startswith("tick@") else tg[0])
    charge[key] += d
    count[key] += 1
    i = j
nt = window[1] - window[0]
print(f"window {window} {'inverse' if inv else 'forward'} depth_pass={depth} ticks {total} per_tick {total / nt:.0f} ops {len(ops)}")
for key, v in charge.most_common(60):
    print(f"  {v / nt:8.1f}  {100 * v / total:5.1f}%  ops {count[key] / nt:7.1f}  {key}")

if "--segments" in sys.argv:
    path = []
    i = end
    while i >= 0:
        path.append(i)
        i = bind[i]
    path.reverse()
    segs = []
    for i in path:
        tg = tags[i].split("#")
        key = tg[-1] if detail else (tg[0].split("@")[0] if not tg[0].startswith("tick@") else tg[0])
        st = fin[bind[i]] if bind[i] >= 0 else 0
        if segs and segs[-1][0] == key:
            segs[-1][2] = fin[i]
            segs[-1][3] += 1
        else:
            segs.append([key, st, fin[i], 1, ops[i], ops[bind[i]] if bind[i] >= 0 else None, tags[i].split("#")[-1]])
    lo, hi = total * 4 // 10, total * 6 // 10
    for key, a, b, n, op, bop, tg in segs:
        if lo <= a <= hi:
            print(f"    {a:7d} {b:7d} {b - a:6d} ops {n:4d}  {key:28s} first {op} <- {bop}  [{tg}]")
