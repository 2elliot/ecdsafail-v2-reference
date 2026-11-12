"""The multiplier M(x, y) = (x, x*y mod q) at primes that take the rank-one start.

    python3 -I -B scripts/m_check.py [--fast]

Builds DialogMultiplier.apply on fresh wires and runs it on every x in [1, q) and every y in [0, q)
at q = 61 (J = 3) and 251 (J = 4), and on 40,000 sampled (x, y) with the edges at q = 4093 and
8191: Y must hold x*y mod q, X must be as found and every other wire 0. It prints each prime's J.
A control flips a bit after each sub_if_ge and must fail.
--fast skips 8191.
"""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import rank_one, sim
from qec.arith import Mod
from qec.builder import Circuit
from qec.point_add import DialogMultiplier


def run(q, lanes):
    mod = Mod(q)
    n = mod.n
    seen = []
    plain = rank_one.RankOne.convert

    def convert(self, *args):
        seen.append(self.J)
        return plain(self, *args)

    rank_one.RankOne.convert = convert
    try:
        C = Circuit()
        X, Y = C.alloc(n), C.alloc(n)
        X0, Y0 = list(X), list(Y)
        DialogMultiplier(C, mod).apply(X, Y)
    finally:
        rank_one.RankOne.convert = plain
    init = sim.pack([x for x, _ in lanes], X0)
    init.update(sim.pack([y for _, y in lanes], Y0))
    state = sim.run(C.ops, C.num_wires, init, len(lanes))
    gx, gy = sim.unpack(state, X, len(lanes)), sim.unpack(state, Y, len(lanes))
    dirty = 0
    live = set(X) | set(Y)
    for w in range(C.num_wires):
        if w not in live:
            dirty |= state[w]
    bad = sum(1 for j, (x, y) in enumerate(lanes) if (gx[j], gy[j]) != (x, x * y % q) or dirty >> j & 1)
    return bad, sorted(set(seen))


def main():
    rng = random.Random(4)
    bad = 0
    for q in (61, 251, 4093) + (() if "--fast" in sys.argv else (8191,)):
        if q < 1000:
            lanes = [(x, y) for x in range(1, q) for y in range(q)]
        else:
            edge = [1, 2, 3, q - 1, q - 2, (q - 1) // 2, (q + 1) // 2]
            lanes = [(x, y) for x in edge for y in edge + [0]]
            lanes += [(rng.randrange(1, q), rng.randrange(q)) for _ in range(40000 - len(lanes))]
        fails, js = run(q, lanes)
        bad += fails
        print(f"q={q}: rank-one start J={js}, {len(lanes)} lanes, {fails} failures", flush=True)
        if not js:
            bad += 1
    # Control: a flipped bit after each sub_if_ge must fail.
    plain = rank_one.sub_if_ge

    def late(C, a, b, flag, ctl=None):
        plain(C, a, b, flag, ctl)
        C.x(b[0])

    rank_one.sub_if_ge = late
    try:
        fails, _ = run(251, [(x, y) for x in range(1, 251) for y in range(0, 251, 5)])
    finally:
        rank_one.sub_if_ge = plain
    print(f"control, a wrong bit after each sub_if_ge at q=251: {fails} failures", flush=True)
    bad += not fails
    print("PASS" if not bad else "FAIL")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
