#!/bin/sh
# Pins the SHA-256 of each instance's `operations` JSON array, the bytes the
# benchmark's compute_operations_hash digests.
set -eu
cd "$(dirname "$0")/.."
PY=${PY:-~/ecdsa_venv/bin/python}
$PY -I - <<'PYEOF' > receipts/baseline/flat_ops.sha256
import hashlib, sys
sys.path.insert(0, ".")
from qec.curve import SECP256K1, SECP256K1_G, toy_secp_curves
from qec.emit import _op
from qec.point_add import ControlledPointAdd

def digest(curve, A):
    ops = ControlledPointAdd(curve, A).build().ops
    h = hashlib.sha256(b"[")
    it = iter(ops)
    h.update(",".join(_op(*t) for t in zip(it, it, it, it)).encode())
    h.update(b"]")
    return h.hexdigest(), len(ops) // 4

toy = next(c for c in toy_secp_curves(7) if c.q == 127)
for name, curve, A in (("toy127", toy, toy.points()[5]), ("secp256k1", SECP256K1, SECP256K1_G)):
    d, n = digest(curve, A)
    print(f"{d}  {name}  ops={n}")
PYEOF
cat receipts/baseline/flat_ops.sha256
