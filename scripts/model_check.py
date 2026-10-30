"""Exhaustive check of the classical dialog model: Bezout identity for all
(x, y) at small q, the proven iteration bound, and the Euclid width schedule."""
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec.model import inplace_mul, iteration_bound, width_schedule


def trace(q, x):
    u, v, sizes = q, x, []
    while v:
        sizes.append(max(u, v).bit_length())
        if v & 1 and u > v:
            u, v = v, u
        if v & 1:
            v -= u
        v >>= 1
    return sizes


rows = []
for q in [3, 5, 7, 11, 13, 31, 43, 61, 127, 251, 509, 1021, 4093]:
    B, W = iteration_bound(q), width_schedule(q, iteration_bound(q))
    worst = 0
    width_ok = True
    for x in range(q):
        t = trace(q, x)
        worst = max(worst, len(t))
        width_ok &= all(s <= W[i] for i, s in enumerate(t))
    ys = range(q) if q < 300 else random.Random(q).sample(range(q), 8)
    mul_ok = all(inplace_mul(q, x, y) == (x * y % q if x else y) for x in range(q) for y in ys)
    rows.append({"q": q, "worst_iterations": worst, "bound": B, "width_schedule_ok": width_ok, "multiplier_ok": mul_ok})
    print(json.dumps(rows[-1]))
q = 2**256 - 2**32 - 977
rng = random.Random(1)
ok = all(inplace_mul(q, x, y) == (x * y % q if x else y) for x, y in ((rng.randrange(q), rng.randrange(q)) for _ in range(2000)))
print(json.dumps({"q": "secp256k1", "bound": iteration_bound(q), "random_pairs": 2000, "multiplier_ok": ok}))
sys.exit(0 if ok and all(r["worst_iterations"] <= r["bound"] and r["width_schedule_ok"] and r["multiplier_ok"] for r in rows) else 1)
