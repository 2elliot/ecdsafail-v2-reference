"""Exhaustive digest of square._rows on every Y, for comparing two trees.

    python3 -I -B scripts/rows_check.py . [nmax]

For leaf widths n = 1..nmax, every output length L = 1..2n, copy in {True, False} and kept-carry
lists K of 0 or enough wires, builds _rows(C, Y, Q, copy, K) on fresh wires and simulates it on all
2^n values of Y at once (lanes), Q and K starting at 0. Prints one line per configuration with the
final Q values, whether Y came back, and the K values, so two trees' outputs can be diffed. Also
checks Q against the contract Y*(2^n - Y) (copy) or Y*(2^n - 1 - Y) mod 2^L directly.
"""
import hashlib, sys
from pathlib import Path

tree = Path(sys.argv[1]).resolve()
nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 10
sys.path.insert(0, str(tree))
from qec import builder, square  # noqa: E402

bad = 0
for n in range(1, nmax + 1):
    lanes = 1 << n
    ALL = (1 << lanes) - 1
    ybits = [sum(((y >> j) & 1) << y for y in range(lanes)) for j in range(n)]
    for L in range(1, 2 * n + 1):
        for copy in (True, False):
            for kn in (0, n * n):
                C = builder.Circuit()
                Y, Q = C.alloc(n), C.alloc(L)
                K = C.alloc(kn) if kn else None
                square._rows(C, Y, Q, copy, K)
                val = [0] * C.num_wires
                for j, w in enumerate(Y):
                    val[w] = ybits[j]
                ops = C.prims
                for i in range(0, len(ops), 4):
                    op, a, b, c = ops[i:i + 4]
                    if op == builder.X: val[a] ^= ALL
                    elif op == builder.CX: val[b] ^= val[a]
                    elif op == builder.CCX: val[c] ^= val[a] & val[b]
                    elif op == builder.SWAP: val[a], val[b] = val[b], val[a]
                    else: raise SystemExit(f"op {op}")
                yback = all(val[w] == ybits[j] for j, w in enumerate(Y))
                qv = [sum(((val[Q[k]] >> y) & 1) << k for k in range(L)) for y in range(lanes)]
                want = [(y * ((1 << n) - y) if copy else y * ((1 << n) - 1 - y)) % (1 << L) for y in range(lanes)]
                kvals = [val[w] for w in (K or [])]
                ok = yback and qv == want
                bad += not ok
                ccx = sum(1 for i in range(0, len(ops), 4) if ops[i] == builder.CCX)
                h = hashlib.sha256(repr((qv, yback, kvals)).encode()).hexdigest()[:16]
                print(f"n={n} L={L} copy={int(copy)} K={kn} ccx={ccx} ok={int(ok)} digest={h}")
print(f"configs failing the contract: {bad}")
