"""Exhaustive check of RankOne._and_nonzero: t ^= ctl*[x != 0].

    python3 -I -B scripts/nonzero_check.py [nmax]

For x of 1 to nmax wires (default 8), runs every (x, ctl, t): t must end t ^ ctl*[x != 0], x and
ctl as found and every other wire 0. Controls: each CCX dropped, or turned into a CX from either
control, must fail.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.builder import CCX, CX, Circuit
from qec.rank_one import RankOne


def failures(ops, wires, x, ctl, t):
    lanes = [(v, c, u) for v in range(1 << len(x)) for c in (0, 1) for u in (0, 1)]
    init = sim.pack([l[0] for l in lanes], x)
    init.update(sim.pack([l[1] for l in lanes], [ctl]))
    init.update(sim.pack([l[2] for l in lanes], [t]))
    state = sim.run(ops, wires, init, len(lanes))
    dirty = 0
    for w in range(wires):
        if w not in x and w not in (ctl, t):
            dirty |= state[w]
    gx, gc, gt = (sim.unpack(state, r, len(lanes)) for r in (x, [ctl], [t]))
    return sum(1 for j, (v, c, u) in enumerate(lanes)
               if (gx[j], gc[j], gt[j]) != (v, c, u ^ (c & (v != 0))) or dirty >> j & 1)


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    bad = mutants = missed = 0
    for n in range(1, nmax + 1):
        C = Circuit()
        x, (ctl, t) = C.alloc(n), C.alloc(2)
        RankOne._and_nonzero(C, x, ctl, t)
        ops = C.prims
        bad += failures(ops, C.num_wires, x, ctl, t)
        for p in range(0, len(ops), 4):
            if ops[p] != CCX:
                continue
            _, a, b, c = ops[p:p + 4]
            for repl in ([], [CX, a, c, c], [CX, b, c, c]):
                mut = ops[:p] + type(ops)("i", repl) + ops[p + 4:]
                mutants += 1
                missed += not failures(mut, C.num_wires, x, ctl, t)
        print(f"n={n}: {4 << n} lanes, {ops[0::4].count(CCX)} CCX", flush=True)
    print(f"{bad} failures; {mutants} mutants, {missed} unnoticed")
    sys.exit(1 if bad or missed else 0)


if __name__ == "__main__":
    main()
