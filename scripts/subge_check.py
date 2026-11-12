"""Exhaustive check of rank_one.sub_if_ge against the compare and controlled subtract it replaces.

    python3 -I -B scripts/subge_check.py [nmax]

For len(b) = 1 to nmax (default 7), len(a) = len(b) and len(b) - 1, with and without ctl, builds
sub_if_ge(C, a, b, flag, ctl) and the old pair (greater; fix the flag; sub under the flag) and runs
every (a, b, ctl) at once with flag clean: both must leave flag = ctl*[b >= a], b - flag*a mod
2^len(b), a and ctl as found, and no other wire set. Controls: each CCX of sub_if_ge dropped, or
turned into a CX from either control, must fail.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.builder import CCX, CX, Circuit
from qec.rank_one import greater, sub, sub_if_ge


def old(C, a, b, flag, ctl):
    greater(C, a, b, flag, ctl=ctl)
    if ctl is None:
        C.x(flag)
    else:
        C.cx(ctl, flag)
    sub(C, a, b, ctl=flag)


def build(fn, la, lb, ctl):
    C = Circuit()
    a, b, (flag,) = C.alloc(la), C.alloc(lb), C.alloc(1)
    c = C.alloc(1)[0] if ctl else None
    fn(C, a, b, flag, c)
    return C, a, b, flag, c


def failures(ops, wires, a, b, flag, c):
    lanes = [(x, y, k) for x in range(1 << len(a)) for y in range(1 << len(b)) for k in ((0, 1) if c is not None else (1,))]
    init = sim.pack([l[0] for l in lanes], a)
    init.update(sim.pack([l[1] for l in lanes], b))
    if c is not None:
        init.update(sim.pack([l[2] for l in lanes], [c]))
    state = sim.run(ops, wires, init, len(lanes))
    live = [*a, *b, flag] + ([c] if c is not None else [])
    dirty = 0
    for w in range(wires):
        if w not in live:
            dirty |= state[w]
    ga, gb, gf = (sim.unpack(state, r, len(lanes)) for r in (a, b, [flag]))
    gc = sim.unpack(state, [c], len(lanes)) if c is not None else [1] * len(lanes)
    bad = 0
    for j, (x, y, k) in enumerate(lanes):
        f = k * (y >= x)
        bad += (ga[j], gb[j], gf[j], gc[j]) != (x, (y - f * x) % (1 << len(b)), f, k) or dirty >> j & 1
    return bad, len(lanes)


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    bad = mutants = missed = saved = 0
    for lb in range(1, nmax + 1):
        for la in sorted({lb, lb - 1} - {0}):
            for ctl in (False, True):
                counts = []
                for fn in (sub_if_ge, old):
                    C, a, b, flag, c = build(fn, la, lb, ctl)
                    fails, lanes = failures(C.prims, C.num_wires, a, b, flag, c)
                    bad += fails
                    counts.append(C.prims[0::4].count(CCX))
                    if fn is sub_if_ge and lb <= 5:
                        ops = C.prims
                        for p in range(0, len(ops), 4):
                            if ops[p] != CCX:
                                continue
                            _, x, y, t = ops[p:p + 4]
                            for repl in ([], [CX, x, t, t], [CX, y, t, t]):
                                mut = ops[:p] + type(ops)("i", repl) + ops[p + 4:]
                                mutants += 1
                                missed += not failures(mut, C.num_wires, a, b, flag, c)[0]
                saved += counts[1] - counts[0]
                print(f"len(b)={lb} len(a)={la} ctl={ctl}: {lanes} lanes, CCX {counts[1]} -> {counts[0]}", flush=True)
    print(f"{bad} failures; {mutants} mutants, {missed} unnoticed")
    sys.exit(1 if bad or missed else 0)


if __name__ == "__main__":
    main()
