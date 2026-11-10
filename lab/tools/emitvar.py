"""Emit the port tree's secp256k1 circuit (hierarchical, with the tree's depth pass) under given flags.

    python3.12 -I -B emitvar.py <tree> <out.json> [MODULE.FLAG=value ...]

e.g. point_add.RIPPLE_MAX_W=0 two_bit.MARK_FROM_SW=False two_bit.CDKM_ADD=True
Values are Python literals. Same as scripts/emit.py --hier otherwise.
"""
import ast
import json
import sys
from pathlib import Path

tree, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
sys.path.insert(0, str(tree))
import importlib  # noqa: E402

for kv in sys.argv[3:]:
    k, v = kv.split("=", 1)
    if k == "two_bit.compare_blocks.max_blocks":
        import qec.two_bit as tb
        tb.compare_blocks.__defaults__ = (int(v),)
        continue
    mod, flag = k.rsplit(".", 1)
    m = importlib.import_module("qec." + mod)
    assert hasattr(m, flag), k
    setattr(m, flag, ast.literal_eval(v))

from qec.curve import SECP256K1, SECP256K1_G  # noqa: E402
from qec.emit import write_circuit  # noqa: E402
from qec.point_add import ControlledPointAdd  # noqa: E402

gen = ControlledPointAdd(SECP256K1, SECP256K1_G)
C = gen.build(True)
with open(out, "w") as fh:
    written = write_circuit(fh, C, gen.layout)
print(json.dumps({"flags": sys.argv[3:], "num_wires": C.num_wires, **written}))
