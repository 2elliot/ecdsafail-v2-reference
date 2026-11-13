"""The fold's probe stretch, old against new, on every input.

    python3 -I -B scripts/probe_check.py [wmax]

RankOne.fold wrote e ^= top*[prod + lo >= 2^W] as: carry_into onto a fresh c1; ccx(top, c1, e);
carry_into again to clear c1. It now calls carry_into(prod, lo, e, ctl=top). The toy primes have no
high bits in the fold, so they never run the controlled form; this check does. For W = 1 to wmax
(default 6) it builds both stretches on the same wires and runs every (prod, lo, top, e): each must
leave e ^ top*carry, prod, lo and top as found and every other wire 0. Controls: the new call
without its control, and each of its CCX dropped or turned into a CX, must differ from that.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.builder import CCX, CX, Circuit
from qec.rank_one import carry_into


def old(C, prod, lo, top, e):
    (c1,) = C.alloc(1)
    carry_into(C, prod, lo, c1)
    C.ccx(top, c1, e)
    carry_into(C, prod, lo, c1)
    C.free([c1])


def new(C, prod, lo, top, e):
    carry_into(C, prod, lo, e, ctl=top)


def uncontrolled(C, prod, lo, top, e):
    carry_into(C, prod, lo, e)


def build(fn, W):
    C = Circuit()
    prod, lo, (top, e) = C.alloc(W), C.alloc(W), C.alloc(2)
    fn(C, prod, lo, top, e)
    return C, prod, lo, top, e


def failures(ops, wires, prod, lo, top, e):
    W = len(lo)
    lanes = [(p, l, t, x) for p in range(1 << W) for l in range(1 << W) for t in (0, 1) for x in (0, 1)]
    init = sim.pack([v[0] for v in lanes], prod)
    init.update(sim.pack([v[1] for v in lanes], lo))
    init.update(sim.pack([v[2] for v in lanes], [top]))
    init.update(sim.pack([v[3] for v in lanes], [e]))
    state = sim.run(ops, wires, init, len(lanes))
    dirty = 0
    for w in range(wires):
        if w not in prod and w not in lo and w not in (top, e):
            dirty |= state[w]
    got = list(zip(*(sim.unpack(state, r, len(lanes)) for r in (prod, lo, [top], [e]))))
    return sum(1 for j, (p, l, t, x) in enumerate(lanes)
               if got[j] != (p, l, t, x ^ (t & (p + l >= 1 << W))) or dirty >> j & 1)


def main():
    wmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    bad = mutants = missed = 0
    for W in range(1, wmax + 1):
        counts = []
        for fn in (old, new):
            C, *regs = build(fn, W)
            bad += failures(C.prims, C.num_wires, *regs)
            counts.append(C.prims[0::4].count(CCX))
        C, *regs = build(uncontrolled, W)
        mutants += 1
        missed += not failures(C.prims, C.num_wires, *regs)
        C, *regs = build(new, W)
        ops = C.prims
        for p in range(0, len(ops), 4):
            if ops[p] != CCX:
                continue
            _, x, y, t = ops[p:p + 4]
            for repl in ([], [CX, x, t, t], [CX, y, t, t]):
                mutants += 1
                missed += not failures(ops[:p] + type(ops)("i", repl) + ops[p + 4:], C.num_wires, *regs)
        print(f"W={W}: {4 << 2 * W} lanes, CCX {counts[0]} -> {counts[1]}", flush=True)
    print(f"{bad} failures; {mutants} mutants, {missed} unnoticed")
    sys.exit(1 if bad or missed else 0)


if __name__ == "__main__":
    main()
