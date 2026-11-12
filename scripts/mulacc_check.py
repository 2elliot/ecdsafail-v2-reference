"""Exhaustive check of rank_one.mul_acc(..., sign=-1, clears=True).

    python3 -I -B scripts/mulacc_check.py [kmax]

For len(acc) = k = 1 to kmax (default 5), with a and b of k wires, runs every (a, b) with acc
holding a*b mod 2^k, the caller's promise: acc must end 0 with a and b as found and no other wire
set, as the plain mul_acc(sign=-1) leaves it. Controls: each CCX dropped, or turned into a CX from
either control, must fail.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.builder import CCX, CX, Circuit
from qec.rank_one import mul_acc


def failures(ops, wires, a, b, acc):
    k = len(acc)
    lanes = [(x, y) for x in range(1 << k) for y in range(1 << k)]
    init = sim.pack([x for x, _ in lanes], a)
    init.update(sim.pack([y for _, y in lanes], b))
    init.update(sim.pack([x * y % (1 << k) for x, y in lanes], acc))
    state = sim.run(ops, wires, init, len(lanes))
    dirty = 0
    for w in range(wires):
        if w not in a and w not in b:
            dirty |= state[w]
    ga, gb = sim.unpack(state, a, len(lanes)), sim.unpack(state, b, len(lanes))
    return sum(1 for j, l in enumerate(lanes) if (ga[j], gb[j]) != l or dirty >> j & 1)


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    bad = mutants = missed = 0
    for k in range(1, kmax + 1):
        counts = []
        for clears in (True, False):
            C = Circuit()
            a, b, acc = C.alloc(k), C.alloc(k), C.alloc(k)
            mul_acc(C, a, b, acc, sign=-1, clears=clears)
            ops = C.prims
            bad += failures(ops, C.num_wires, a, b, acc)
            counts.append(ops[0::4].count(CCX))
            if clears:
                for p in range(0, len(ops), 4):
                    if ops[p] != CCX:
                        continue
                    _, x, y, t = ops[p:p + 4]
                    for repl in ([], [CX, x, t, t], [CX, y, t, t]):
                        mut = ops[:p] + type(ops)("i", repl) + ops[p + 4:]
                        mutants += 1
                        missed += not failures(mut, C.num_wires, a, b, acc)
        print(f"k={k}: {1 << 2 * k} lanes, CCX {counts[1]} -> {counts[0]}", flush=True)
    print(f"{bad} failures; {mutants} mutants, {missed} unnoticed")
    sys.exit(1 if bad or missed else 0)


if __name__ == "__main__":
    main()
