"""Checks of the tick's carry-select compare, two_bit.compare_blocks: t ^= z * [A > B], with A, B, z,
cin, the clean wires and the dirty wires (any values) all restored.

    python3.12 -I -B scripts/compare_check.py <tree> [wmax] [draws]

Exhaustive: for widths 2..wmax (default 6), every block bound list with up to 4 blocks, on all A, B,
z and every value of the dirty wires the blocks borrow. At width 255 with the tree's CUTS_255 for
MAX_BLOCKS blocks: A = B, and for each bit k where they first differ from the top, both orders, with
`draws` random lower bits (default 16), random dirty values and z = 0 and 1.
"""
import itertools
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
from qec import sim, two_bit  # noqa: E402
from qec.builder import Circuit  # noqa: E402

wmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
draws = int(sys.argv[3]) if len(sys.argv) > 3 else 16


def build(width, nclean, ndirty, bounds=None):
    C = Circuit()
    A, B = C.alloc(width), C.alloc(width)
    z, t, cin = C.alloc(3)
    clean, dirty = C.alloc(nclean), C.alloc(ndirty)
    saved = two_bit.weighted_bounds
    if bounds is not None:
        two_bit.weighted_bounds = lambda w, blocks: list(bounds)
    try:
        two_bit.compare_blocks(C, A, B, z, t, cin, clean, dirty, {}, max_blocks=two_bit.MAX_BLOCKS
                               if bounds is None else len(bounds) - 1)
    finally:
        two_bit.weighted_bounds = saved
    return C, A, B, z, t, cin, clean, dirty


def run(width, cases, nclean, ndirty, bounds=None):
    """cases: (a, b, z, dirty value). Returns failing lanes plus wires left dirty."""
    C, A, B, z, t, cin, clean, dirty = build(width, nclean, ndirty, bounds)
    n = len(cases)
    init = sim.pack([c[0] for c in cases], A)
    init.update(sim.pack([c[1] for c in cases], B))
    init.update(sim.pack([c[2] for c in cases], [z]))
    init.update(sim.pack([c[3] for c in cases], dirty))
    st = sim.run(C.ops, C.num_wires, init, n)
    keep = A + B + [z] + dirty
    want = {w: init.get(w, 0) for w in keep}
    tv = sim.unpack(st, [t], n)
    bad = sum(tv[i] != (c[2] & int(c[0] > c[1])) for i, c in enumerate(cases))
    bad += sum(st[w] != want[w] for w in keep)
    bad += sum(1 for w in [cin] + clean if st[w])
    bad += sum(1 for w in range(C.num_wires) if w not in set(keep) | {t, cin} | set(clean) and st[w])
    return bad


ok = True
for width in range(2, wmax + 1):
    lists = [(0, *cut, width) for k in range(0, min(3, width - 1) + 1)
             for cut in itertools.combinations(range(1, width), k)]
    fails = 0
    for bounds in lists:
        m = len(bounds) - 2
        # carry-ins all borrowed dirty, then all clean
        for nclean, ndirty in ((m, m), (2 * m, 0)):
            cases = [(a, b, zz, d) for a in range(1 << width) for b in range(1 << width) for zz in (0, 1)
                     for d in range(1 << ndirty)]
            fails += run(width, cases, nclean, ndirty, bounds) > 0
    ok &= fails == 0
    print(f"width={width} bound lists {len(lists)} failing {fails}", flush=True)

width, blocks = 255, two_bit.MAX_BLOCKS
rng = random.Random(2026)
cases = []
for zz in (0, 1):
    for _ in range(4):
        a = rng.getrandbits(width)
        cases.append((a, a, zz, rng.getrandbits(10)))
    for k in range(width):
        for _ in range(draws):
            hi = rng.getrandbits(width) >> (k + 1) << (k + 1)
            lo_a, lo_b = rng.getrandbits(k) if k else 0, rng.getrandbits(k) if k else 0
            a, b = hi | 1 << k | lo_a, hi | lo_b
            d = rng.getrandbits(10)
            cases += [(a, b, zz, d), (b, a, zz, d)]
bad = run(width, cases, 10, 10)
ok &= bad == 0
print(f"width=255 blocks {blocks} cuts {two_bit.CUTS_255[blocks]} lanes {len(cases)} failing {bad}", flush=True)
print("compare check:", "OK" if ok else "FAIL")
sys.exit(0 if ok else 1)
