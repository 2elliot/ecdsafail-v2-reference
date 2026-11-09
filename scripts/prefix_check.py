"""Checks of the tick's P, two_bit.PrefixR: after PrefixR(C, bits, sizes) its top holds AND(bits), and
write(t) leaves bits + t mod 2^h with t unchanged and every other wire 0.

    python3.12 -I -B scripts/prefix_check.py <tree> [hmax] [draws]

Exhaustive: every block-size list PrefixR accepts for h = 2..hmax (default 9), on all 2^(h+1) inputs.
At h = 221 with the tree's P_SIZES_221: for each position k of the lowest 0 bit (and all ones), t = 0
and 1, `draws` random fillings above it (default 24), which reach every block's carry and junk path.
"""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
from qec import sim, two_bit  # noqa: E402
from qec.builder import Circuit  # noqa: E402

hmax = int(sys.argv[2]) if len(sys.argv) > 2 else 9
draws = int(sys.argv[3]) if len(sys.argv) > 3 else 24


def build(h, sizes):
    C = Circuit()
    bits, (t,) = C.alloc(h), C.alloc(1)
    P = two_bit.PrefixR(C, bits, sizes=list(sizes))
    (top,) = C.alloc(1)
    C.cx(P.top, top)
    P.write(t)
    return C, bits, t, top


def run(h, sizes, cases):
    """cases: (bits value, t). Returns the number of failing lanes."""
    C, bits, t, top = build(h, sizes)
    init = sim.pack([v for v, _ in cases], bits)
    init.update(sim.pack([u for _, u in cases], [t]))
    st = sim.run(C.ops, C.num_wires, init, len(cases))
    got = sim.unpack(st, bits, len(cases))
    tops = sim.unpack(st, [top], len(cases))
    ts = sim.unpack(st, [t], len(cases))
    mask = (1 << h) - 1
    bad = 0
    others = [w for w in range(C.num_wires) if w not in set(bits) | {t, top}]
    for lane, (v, u) in enumerate(cases):
        ok = got[lane] == (v + u) & mask and ts[lane] == u and tops[lane] == int(v == mask)
        bad += not ok
    dirty = sum(1 for w in others if st[w])
    return bad + dirty


def size_lists(h):
    """Every sizes list PrefixR accepts at width h: block m <= carriers + 1, blocks fill h."""
    out = []

    def go(start, ntops, acc):
        if start == h:
            out.append(tuple(acc))
            return
        for m in range(1, min(start + ntops, h - start) + 1):
            go(start + m, ntops + 1, acc + [m])

    go(2, 1, [])
    return out


ok = True
for h in range(2, hmax + 1):
    lists = size_lists(h)
    cases = [(v, u) for u in (0, 1) for v in range(1 << h)]
    fails = sum(run(h, s, cases) > 0 for s in lists)
    ok &= fails == 0
    print(f"h={h} size lists {len(lists)} lanes {len(cases)} failing lists {fails}", flush=True)

h, sizes = 221, two_bit.P_SIZES_221
rng = random.Random(2026)
mask = (1 << h) - 1
cases = []
for k in range(h + 1):
    for u in (0, 1):
        for _ in range(draws if k < h else 1):
            hi = rng.getrandbits(h) & ~((1 << (k + 1)) - 1) & mask if k < h else 0
            cases.append(((1 << k) - 1 | hi, u))
bad = run(h, sizes, cases)
ok &= bad == 0
print(f"h=221 sizes {tuple(sizes)} lanes {len(cases)} failing {bad}", flush=True)
print("prefix check:", "OK" if ok else "FAIL")
sys.exit(0 if ok else 1)
