"""Exhaustive check of SquareFold.compute with an absorbed constant, at small primes.

    python3 -I -B scripts/fold_check.py [max_bits]

For every prime q = 2^n - f with f*f <= q, n = 4 to max_bits (default 11), and every c below q that
absorb takes, and once with no absorbed constant, builds compute(C, Y, J, kappa, cins) and its reverse and runs every n-bit Y at once:
J[:n] must hold (Y^2 - c) mod q after compute, Y must be as it was, and every wire of J must be 0
after the reverse. At secp256k1 it samples Y, with the edges, for c = 3a and two random c. It counts the builds whose second fold loads a nonzero constant (klo), which is
the path where _multiple copies E and adds klo as a constant. Controls, each of which must fail on
some build: that add one bit short, the constant's lowest set bit dropped, and a digit add that
takes its target as clean one bit below where the bound says.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim, square
from qec.arith import Mod
from qec.builder import Circuit
from qec.square import SquareFold


def is_prime(v):
    return v > 1 and all(v % p for p in range(2, int(v ** 0.5) + 1))


def run(q, c, ys=None):
    """(failures, loaded) over every Y, or over ys, or None when absorb does not take c."""
    mod = Mod(q)
    n = mod.n
    fold = SquareFold(mod)
    if c is None:
        # No absorbed constant: compute leaves (Y^2 - offset) mod q.
        cins, kappa, c = frozenset(), None, fold.offset
    else:
        plan = fold.absorb(c)
        if plan is None:
            return None
        cins, kappa = plan
    C = Circuit()
    Y, J = C.alloc(n), C.alloc(fold.junk)
    blk = C.record(lambda y, j: fold.compute(C, y, j, kappa, cins), [Y, J])
    ys = list(range(1 << n)) if ys is None else ys
    lanes = len(ys)
    C.apply(blk, [Y, J], inline=True)
    mid = len(C.prims)
    C.apply(blk, [Y, J], inverse=True, inline=True)
    state = sim.run(C.prims[:mid], C.num_wires, sim.pack(ys, Y), lanes)
    got, back = sim.unpack(state, J[:n], lanes), sim.unpack(state, Y, lanes)
    bad = sum(1 for i, y in enumerate(ys) if got[i] != (y * y - c) % q or back[i] != y)
    state = sim.run(C.prims, C.num_wires, sim.pack(ys, Y), lanes)
    dirty = 0
    for w in range(C.num_wires):
        if w not in Y:
            dirty |= state[w]
    bad += dirty.bit_count() + sum(1 for y, v in zip(ys, sim.unpack(state, Y, lanes)) if v != y)
    return bad, kappa is not None and bool(fold.split(kappa)[1])


def sweep(bits, mutate=None):
    builds = loaded = bad = 0
    for n in range(4, bits + 1):
        for f in range(1, 1 << (n // 2 + 1), 2):
            q = (1 << n) - f
            if not is_prime(q) or not SquareFold.applies(q):
                continue
            for c in [None, *(range(q) if n <= 8 else range(0, q, max(1, q // 97)))]:
                out = run(q, c)
                if out is None:
                    continue
                builds += 1
                loaded += out[1]
                bad += out[0]
    return builds, loaded, bad


def main():
    bits = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    builds, loaded, bad = sweep(bits)
    print(f"{builds} builds, {loaded} with a loaded constant, {bad} failures", flush=True)
    ok = not bad and loaded > 0
    import random
    from qec.curve import SECP256K1, SECP256K1_G
    q, rng = SECP256K1.q, random.Random(11)
    ys = [rng.randrange(1 << 256) for _ in range(500)] + [0, 1, 2, q - 1, q, q + 1, (1 << 256) - 1, 1 << 255]
    ys += [(1 << 256) - 1 - rng.randrange(1 << k) for k in range(1, 256, 5)] + [rng.randrange(1 << k) for k in range(1, 256, 5)]
    for c in (3 * SECP256K1_G[0] % q, rng.randrange(q), rng.randrange(q)):
        out = run(q, c, ys)
        print(f"secp256k1 c={c:#x}: " + ("not absorbed" if out is None else f"{len(ys)} lanes, loaded={out[1]}, {out[0]} failures"), flush=True)
        ok &= out is None or not out[0]
    plain = square._add_loaded

    def short(C, reg, value):
        plain(C, reg[:-1], value)

    def dropped(C, reg, value):
        plain(C, reg, value & (value - 1))

    clean = square._clean_from
    for name, attr, fn in (("constant add one bit short", "_add_loaded", short),
                           ("constant's lowest set bit dropped", "_add_loaded", dropped),
                           ("a digit's target taken as clean one bit early", "_clean_from", lambda v: max(1, clean(v) - 1))):
        saved = getattr(square, attr)
        setattr(square, attr, fn)
        try:
            _, _, missed = sweep(min(bits, 8))
        finally:
            setattr(square, attr, saved)
        print(f"control, {name}: {missed} failures", flush=True)
        ok &= missed > 0
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
