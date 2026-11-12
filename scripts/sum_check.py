"""Exhaustive check of square._sum, the out-of-place S = A + B of a Karatsuba node.

    python3 -I -B scripts/sum_check.py [hmax]

For every width h = 1 to hmax (default 8) and len(B) = h or h - 1, the two shapes a node gives it,
builds _sum(C, A, B, S) on a clean S of h + 1 wires and runs every (A, B) at once: S must hold
A + B, A and B must be as they were, and the build must use h CCX at most and no other wire.
Controls: each CCX dropped in turn, and each turned into a CX from either control, must fail.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.builder import CCX, CX, Circuit
from qec.square import _sum


def failures(ops, wires, A, B, S):
    h, m = len(A), len(B)
    pairs = [(a, b) for a in range(1 << h) for b in range(1 << m)]
    init = sim.pack([a for a, _ in pairs], A)
    init.update(sim.pack([b for _, b in pairs], B))
    state = sim.run(ops, wires, init, len(pairs))
    got = zip(sim.unpack(state, A, len(pairs)), sim.unpack(state, B, len(pairs)), sim.unpack(state, S, len(pairs)))
    return sum(1 for (a, b), out in zip(pairs, got) if out != (a, b, a + b))


def main():
    hmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    bad = mutants = missed = 0
    for h in range(1, hmax + 1):
        for m in {h, max(h - 1, 0)}:
            C = Circuit()
            A, B, S = C.alloc(h), C.alloc(m), C.alloc(h + 1)
            _sum(C, A, B, S)
            ops = C.prims
            ccx = [p for p in range(0, len(ops), 4) if ops[p] == CCX]
            assert len(ccx) <= h and C.num_wires == 2 * h + m + 1
            fails = failures(ops, C.num_wires, A, B, S)
            bad += fails
            for p in ccx:
                _, x, y, t = ops[p:p + 4]
                for repl in ([], [CX, x, t, t], [CX, y, t, t]):
                    mut = ops[:p] + type(ops)("i", repl) + ops[p + 4:]
                    mutants += 1
                    missed += not failures(mut, C.num_wires, A, B, S)
            print(f"h={h} len(B)={m}: {len(ccx)} CCX, {1 << (h + m)} lanes, {fails} failures", flush=True)
    print(f"{bad} failures; {mutants} mutants, {missed} unnoticed")
    sys.exit(1 if bad or missed else 0)


if __name__ == "__main__":
    main()
