"""Write harness lanes (c, R) with expected R + c*A from the reference curve law."""
import argparse
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec.curve import SECP256K1, SECP256K1_G, toy_secp_curves

ap = argparse.ArgumentParser()
ap.add_argument("--curve", default="secp256k1")
ap.add_argument("--addend-index", type=int, default=1)
ap.add_argument("--random", type=int, default=250)
ap.add_argument("--seed", type=int, default=99)
ap.add_argument("--out", type=Path, required=True)
args = ap.parse_args()
rng = random.Random(args.seed)
if args.curve == "secp256k1":
    curve, A = SECP256K1, SECP256K1_G
    twoA = curve.add(A, A)
    pts = [None, A, curve.neg(A), curve.neg(twoA), twoA] + [curve.mul(rng.randrange(1, 2**256), A) for _ in range(args.random)]
else:
    q = int(args.curve.split(":")[1])
    curve = next(c for c in toy_secp_curves(q.bit_length()) if c.q == q)
    pts = curve.points()
    A = pts[args.addend_index]
enc = lambda P: ("0x0", "0x0") if P is None else (hex(P[0]), hex(P[1]))
lanes = []
for c in (0, 1):
    for R in pts:
        out = curve.add(R, A) if c else R
        (x, y), (ex, ey) = enc(R), enc(out)
        lanes.append({"c": c, "x": x, "y": y, "ex": ex, "ey": ey})
pk = curve.mul(5, SECP256K1_G) if curve is SECP256K1 else SECP256K1_G
args.out.write_text(json.dumps({"public_key": [hex(pk[0]), hex(pk[1])], "lanes": lanes}))
print(len(lanes), "lanes")
