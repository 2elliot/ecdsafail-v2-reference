"""Exhaustive check of rank_one.add's carry in, on the residue's N -> U - N under ctl.

    python3 -I -B scripts/addcin_check.py [nmax]

For n = 1 to nmax (default 6) builds the old stretch (complement N under ctl; inc under ctl; add U
under ctl) and the new one (complement; add U under ctl with cin = ctl) and runs every (U, N, ctl):
both must leave N = (U - N) mod 2^n under ctl and N otherwise, U and ctl as found, every other wire
0. It also runs add(a, b, cin=c) with no control against b + a + c. Controls: the new stretch
without its carry in, and each of its CCX dropped or turned into a CX, must fail.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.builder import CCX, CX, Circuit
from qec.rank_one import add, inc


def old(C, U, N, ctl):
    for x in N:
        C.cx(ctl, x)
    inc(C, N, ctl)
    add(C, U, N, ctl=ctl)


def new(C, U, N, ctl):
    for x in N:
        C.cx(ctl, x)
    add(C, U, N, ctl=ctl, cin=ctl)


def no_cin(C, U, N, ctl):
    for x in N:
        C.cx(ctl, x)
    add(C, U, N, ctl=ctl)


def plain(C, U, N, ctl):
    add(C, U, N, cin=ctl)


def failures(fn, n, want, mutate=None):
    C = Circuit()
    U, N, (ctl,) = C.alloc(n), C.alloc(n), C.alloc(1)
    fn(C, U, N, ctl)
    ops = C.prims if mutate is None else mutate(C.prims)
    lanes = [(u, v, c) for u in range(1 << n) for v in range(1 << n) for c in (0, 1)]
    init = sim.pack([l[0] for l in lanes], U)
    init.update(sim.pack([l[1] for l in lanes], N))
    init.update(sim.pack([l[2] for l in lanes], [ctl]))
    state = sim.run(ops, C.num_wires, init, len(lanes))
    dirty = 0
    for w in range(C.num_wires):
        if w not in U and w not in N and w != ctl:
            dirty |= state[w]
    got = list(zip(*(sim.unpack(state, r, len(lanes)) for r in (U, N, [ctl]))))
    bad = sum(1 for j, (u, v, c) in enumerate(lanes) if got[j] != (u, want(u, v, c, n), c) or dirty >> j & 1)
    return bad, C.prims


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    reflect = lambda u, v, c, n: (u - v) % (1 << n) if c else v
    bad = mutants = missed = 0
    for n in range(1, nmax + 1):
        b_old, ops_old = failures(old, n, reflect)
        b_new, ops = failures(new, n, reflect)
        b_plain, _ = failures(plain, n, lambda u, v, c, n: (u + v + c) % (1 << n))
        bad += b_old + b_new + b_plain
        if n > 1:
            mutants += 1
            missed += not failures(no_cin, n, reflect)[0]
        for p in range(0, len(ops), 4):
            if ops[p] != CCX:
                continue
            _, x, y, t = ops[p:p + 4]
            for repl in ([], [CX, x, t, t], [CX, y, t, t]):
                mutants += 1
                missed += not failures(new, n, reflect, lambda o: o[:p] + type(o)("i", repl) + o[p + 4:])[0]
        print(f"n={n}: {2 << 2 * n} lanes, CCX {ops_old[0::4].count(CCX)} -> {ops[0::4].count(CCX)}", flush=True)
    print(f"{bad} failures; {mutants} mutants, {missed} unnoticed")
    sys.exit(1 if bad or missed else 0)


if __name__ == "__main__":
    main()
