"""Exhaustive check of the tree forms of the top-level constant adds (arith.CONST_TREE).

    python3.12 -I -B scripts/const_tree_check.py <tree> <nmax> [broken]

gt_const_tree: for n = 1..nmax, every value < 2^n, every b < 2^n, t and ctrl in {0, 1}, uncontrolled
too: t must end t ^ ctrl*[value > b], b and ctrl restored. Mod.cmodadd_const and Mod.creflect: for
every odd q with 3 <= q < 2^nmax, every constant < q, every canonical b and ctrl in {0, 1}: b must
end (b + ctrl*v) mod q, or (a - b) mod q when ctrl, and ctrl restored. Every other wire must end 0.
With "broken", _carry_tree skips its last join; the check must fail.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
from qec import arith, sim
from qec.builder import Circuit

NMAX = int(sys.argv[2])
arith.CONST_TREE = True
if len(sys.argv) > 3 and sys.argv[3] == "broken":
    orig = arith._carry_tree

    def broken(g, G, P, q):
        rec = arith._Gates()
        orig(rec, G, P, q)
        last = max((i for i, o in enumerate(rec.ops) if len(o) == 3 and o[2] == G[-1]), default=None)
        for i, o in enumerate(rec.ops):
            if i != last:
                (g.cx if len(o) == 2 else g.ccx)(*o)
    arith._carry_tree = broken


def run(build, regs_vals):
    """Run build on fresh registers for every lane; return (final values per register, dirty)."""
    C = Circuit()
    regs = build(C)
    lanes = len(regs_vals)
    init = {}
    for r, reg in enumerate(regs):
        for lane, vals in enumerate(regs_vals):
            for i, w in enumerate(reg):
                if vals[r] >> i & 1:
                    init[w] = init.get(w, 0) | 1 << lane
    st = sim.run(C.ops, C.num_wires, init, lanes)
    out = [sim.unpack(st, reg, lanes) for reg in regs]
    live = {w for reg in regs for w in reg}
    dirty = any(st[w] for w in range(C.num_wires) if w not in live)
    return [list(zip(*out))[lane] for lane in range(lanes)], dirty


bad = 0
cases = 0
for n in range(1, NMAX + 1):
    for value in range(1 << n):
        for ctrl in (True, False):
            def build(C, value=value, ctrl=ctrl, n=n):
                c, b, t = C.alloc(1), C.alloc(n), C.alloc(1)
                arith.gt_const_tree(C, value, b, t[0], c[0] if ctrl else None)
                return [c, b, t]
            vals = [(c, b, t) for c in ((0, 1) if ctrl else (0,)) for b in range(1 << n) for t in (0, 1)]
            outs, dirty = run(build, vals)
            for (c, b, t), got in zip(vals, outs):
                want = (c, b, t ^ ((1 if ctrl is False else c) & (value > b)))
                bad += got != want
            bad += dirty
            cases += len(vals)
print("gt_const_tree cases", cases, "bad", bad)

for q in range(3, 1 << NMAX, 2):
    n = q.bit_length()
    mod = arith.Mod(q)
    vals = [(c, b) for c in (0, 1) for b in range(q)]
    for v in range(q):
        def build(C, v=v):
            c, b = C.alloc(1), C.alloc(n)
            mod.cmodadd_const(C, c[0], v, b)
            return [c, b]
        outs, dirty = run(build, vals)
        bad += dirty + sum(got != (c, (b + c * v) % q) for (c, b), got in zip(vals, outs))

        def build_r(C, v=v):
            c, b = C.alloc(1), C.alloc(n)
            mod.creflect(C, c[0], v, b)
            return [c, b]
        outs, dirty = run(build_r, vals)
        bad += dirty + sum(got != (c, (v - b) % q if c else b) for (c, b), got in zip(vals, outs))
        cases += 2 * len(vals)
print("cases", cases, "bad", bad)
print("OK" if bad == 0 else "FAIL")
sys.exit(0 if bad == 0 else 1)
