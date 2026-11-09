"""Exhaustive check of the Bezout tick's add-or-shift (port of 2elliot depth/v2's
scripts/add_or_shift_check.py to the bez-jack tree, where the add is a module function).

    python3.12 -I -B scripts/add_or_shift_check.py <tree> <nmax> [la|cdkm|tick ...]

la = two_bit.add_or_shift_la (PR8), cdkm = two_bit.add_or_shift_cdkm (port), tick = whatever
TwoBitMod._add_or_shift dispatches to under the tree's flags, wide = add_or_shift_la with every
composition of n - 1 into F blocks of 1 to 4 bits, wideb = the same for B. For n = 4..nmax builds
the add alone on J, R (n), word (n + 2) with fresh scratch, simulates all 2^(2n+3) assignments, and
requires on every lane: J and R unchanged, every scratch wire 0, word = word + (R >> 1) mod 2^(n+2)
when J = 0, and word[:n + 1] rotated left by one with word[n + 1] unchanged when J = 1. Expected
words are computed bit-sliced with no generator code.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
from qec import sim, two_bit  # noqa: E402
from qec.builder import Circuit  # noqa: E402

PRIMES = {4: 13, 5: 31, 6: 61, 7: 127, 8: 251, 9: 509, 10: 1021, 11: 2039}
impls = sys.argv[3:] or ["la", "cdkm"]


def lane_mask(i, N):
    L = 1 << N
    period = 1 << (i + 1)
    block = ((1 << (1 << i)) - 1) << (1 << i)
    reps = ((1 << L) - 1) // ((1 << period) - 1)
    return reps * block


def compositions(m, top=4):
    if m == 0:
        yield []
        return
    for r in range(1, min(top, m) + 1):
        for rest in compositions(m - r, top):
            yield [r] + rest


runs = []
for impl in impls:
    for n in range(4, int(sys.argv[2]) + 1):
        if impl in ("wide", "wideb"):
            runs += [(impl, n, sizes) for sizes in compositions(n - 1)]
        else:
            runs.append((impl, n, None))

ok = True
for impl, n, sizes in runs:
    C = Circuit()
    (J,) = C.alloc(1)
    R, word = C.alloc(n), C.alloc(n + 2)
    ins = [J] + R + word
    if impl == "la":
        two_bit.add_or_shift_la(C, J, R, word, n)
    elif impl == "wide":
        two_bit.add_or_shift_la(C, J, R, word, n, fsizes=sizes)
    elif impl == "wideb":
        two_bit.add_or_shift_la(C, J, R, word, n, bsizes=sizes)
    elif impl == "cdkm":
        two_bit.add_or_shift_cdkm(C, J, R, word, n)
    else:
        two_bit.TwoBitMod(PRIMES[n])._add_or_shift(C, J, R, word)
    N = len(ins)
    full = (1 << (1 << N)) - 1
    init = {w: lane_mask(i, N) for i, w in enumerate(ins)}
    st = sim.run(C.ops, C.num_wires, init, 1 << N)
    j, r, w = init[J], [init[x] for x in R], [init[x] for x in ins[n + 1:]]
    addend = r[1:] + [0, 0, 0]
    total, c = [], 0
    for x, y in zip(w, addend):
        total.append(x ^ y ^ c)
        c = (x & y) | (c & (x ^ y))
    rot = [w[n]] + w[:n] + [w[n + 1]]
    want = [(j & u) | (~j & full & v) for u, v in zip(rot, total)]
    outs = set(ins)
    scratch = [x for x in range(C.num_wires) if x not in outs]
    good = (st[J] == j and [st[x] for x in R] == r and [st[x] for x in word] == want
            and all(st[x] == 0 for x in scratch))
    ok &= good
    if impl not in ("wide", "wideb") or not good:
        print(f"{impl:4s} n={n} lanes=2^{N} scratch={len(scratch)} ccx={C.ops[0::4].count(2)} "
              f"map_and_clean={good}{' F ' + str(sizes) if sizes else ''}", flush=True)
for impl in ("wide", "wideb"):
    if impl in impls:
        print(impl, "layouts checked:", sum(1 for r in runs if r[0] == impl))
print("add_or_shift exhaustive:", "OK" if ok else "FAIL")
sys.exit(0 if ok else 1)
