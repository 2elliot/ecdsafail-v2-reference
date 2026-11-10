"""Screen a tree's Bezout pass: build the circuit (hierarchical, under flags), take M's body from its
first Bezout tick call to its last, flatten it with the depth pass applied to every definition it calls,
and run the board's tick model on it with all wires free at tick 0, forward and reversed (M runs once
each way). Also prints the whole circuit's wires and CCX from the build.

    python3.12 -I -B bzpass.py <tree> [MODULE.FLAG=value ...] [--window=a:b]

--window=a:b flattens only tick calls a..b-1 (counted in emission order) for a quicker read.
"""
import ast
import importlib
import sys
import time
from pathlib import Path

tree = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(tree))
window = None
for kv in sys.argv[2:]:
    if kv.startswith("--window="):
        a, b = kv.split("=", 1)[1].split(":")
        window = (int(a), int(b))
        continue
    k, v = kv.split("=", 1)
    if k == "two_bit.compare_blocks.max_blocks":
        import qec.two_bit as tb
        tb.compare_blocks.__defaults__ = (int(v),)
        continue
    mod, flag = k.rsplit(".", 1)
    m = importlib.import_module("qec." + mod)
    assert hasattr(m, flag), k
    setattr(m, flag, ast.literal_eval(v))

from qec import depth_pass  # noqa: E402
from qec.builder import CCX, CX, SWAP, X  # noqa: E402
from qec.curve import SECP256K1, SECP256K1_G  # noqa: E402
from qec.point_add import ControlledPointAdd  # noqa: E402

t0 = time.time()
gen = ControlledPointAdd(SECP256K1, SECP256K1_G)
C = gen.build(True)
memo = {}


def nccx(block):
    if id(block) not in memo:
        memo[id(block)] = block.ops[0::4].count(CCX) + sum(nccx(c.block) for c in block.calls)
    return memo[id(block)]


total_ccx = C.prims[0::4].count(CCX) + sum(nccx(c.block) for c in C.calls)
M = gen.M
body = M.block
ticks_set = {id(b) for b in M._bezout_ticks.values()}
idx = [i for i, c in enumerate(body.calls) if id(c.block) in ticks_set]
lo_i, hi_i = idx[0], idx[-1]
if window:
    lo_i, hi_i = idx[window[0]], idx[min(window[1], len(idx)) - 1]
region = body.calls[lo_i:hi_i + 1]

# definitions reachable from the region, callees first
order, seen = [], set()


def visit(calls):
    for c in calls:
        if id(c.block) not in seen:
            seen.add(id(c.block))
            visit(c.block.calls)
            order.append(c.block)


visit(region)
ids = {b: i for i, b in enumerate(order)}
bodies = depth_pass.optimize(order, ids)

out = []


def expand(block, table, inverse):
    ops, calls = bodies[block]
    wires = sorted(block.wires)
    # local (builder) wires -> caller wires
    m = dict(zip(wires, table))
    seq = []
    pos = 0
    for c in calls:
        seq.append(("p", pos, c.pos))
        seq.append(("c", c))
        pos = c.pos
    seq.append(("p", pos, len(ops) // 4))
    for item in (reversed(seq) if inverse else seq):
        if item[0] == "p":
            _, a, b = item
            rng = range(b - 1, a - 1, -1) if inverse else range(a, b)
            for i in rng:
                k, x, y, z = ops[4 * i:4 * i + 4]
                out.append((k, m[x], m[y], m[z]))
        else:
            c = item[1]
            expand(c.block, [m[w] for w in c.wires], inverse != c.inverse)


pos = region[0].pos
ntick = 0
for c in region:
    for i in range(pos, c.pos):
        out.append(tuple(body.ops[4 * i:4 * i + 4]))
    expand(c.block, c.wires, c.inverse)
    ntick += id(c.block) in ticks_set
    pos = c.pos


def run(ops):
    r = {}
    for k, a, b, c in ops:
        if k == X:
            r[a] = r.get(a, 0) + 1
        elif k == CX or k == SWAP:
            t = max(r.get(a, 0), r.get(b, 0)) + (1 if k == CX else 3)
            r[a] = r[b] = t
        else:
            ra, rb, rt = r.get(a, 0), r.get(b, 0), r.get(c, 0)
            s2 = max(rb, rt + 1) + 1
            s4 = max(ra, s2 + 1) + 1
            r[c] = s4 + 6
            r[a] = r[b] = s4 + 7
    return max(r.values(), default=0)


fwd = run(out)
inv = run(reversed(out))
ccx_region = sum(1 for o in out if o[0] == CCX)
print(f"wires {C.num_wires} CCX {total_ccx} ticks_called {ntick} region_ccx {ccx_region} "
      f"fwd {fwd} inv {inv} sum {fwd + inv} per_tick {(fwd + inv) / 2 / ntick:.1f} "
      f"[{time.time() - t0:.0f}s] {' '.join(sys.argv[2:])}", flush=True)
