"""Emit the controlled point-addition circuit as v2-ir ArithmeticCircuit JSON."""
import argparse
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec.curve import SECP256K1, SECP256K1_G, toy_secp_curves
from qec.emit import definitions, write_circuit
from qec.point_add import ControlledPointAdd

ap = argparse.ArgumentParser()
ap.add_argument("--curve", default="secp256k1", help="secp256k1 or toy:<q>")
ap.add_argument("--addend-index", type=int, default=1, help="toy curves: index into the point list")
ap.add_argument("--out", type=Path, required=True)
ap.add_argument("--hier", action="store_true", help="emit recorded blocks as definitions and calls")
args = ap.parse_args()

if args.curve == "secp256k1":
    curve, A = SECP256K1, SECP256K1_G
else:
    q = int(args.curve.split(":")[1])
    curve = next(c for c in toy_secp_curves(q.bit_length()) if c.q == q)
    A = curve.points()[args.addend_index]
gen = ControlledPointAdd(curve, A)
C = gen.build(args.hier)
opener = gzip.open if args.out.suffix == ".gz" else open
with opener(args.out, "wt") as fh:
    nbytes = write_circuit(fh, C, gen.layout)
meta = {"curve": curve.name, "q": hex(curve.q), "addend": [hex(A[0]), hex(A[1])], "n": gen.n, "num_wires": C.num_wires}
if args.hier:
    defs = definitions(C)
    meta.update(definitions=len(defs), call_sites=len(C.calls) + sum(len(b.calls) for b in defs),
                primitive_ops=len(C) + sum(len(b.ops) // 4 for b in defs))
else:
    meta.update(ops=len(C), ccx=C.ops[0::4].count(2))
meta.update(json_bytes=nbytes, layout=gen.layout)
args.out.with_name(args.out.name.split(".json")[0] + ".meta.json").write_text(json.dumps(meta))
print(json.dumps({k: v for k, v in meta.items() if k != "layout"}))
