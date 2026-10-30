"""Negative control: the checker must reject single-gate mutations."""
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from check import check
from qec.builder import CCX, CX
from qec.curve import toy_secp_curves
from qec.point_add import ControlledPointAdd

rng = random.Random(7)
curve = toy_secp_curves(7)[3]
pts = curve.points()
A = pts[5]
gen = ControlledPointAdd(curve, A)
C = gen.build()
cases = [(c, R) for c in (0, 1) for R in pts]
assert not check(gen, C.ops, cases), "unmutated circuit must pass"
caught = 0
trials = 200
ccx_idx = [i for i in range(len(C)) if C.ops[4 * i] == CCX]
for _ in range(trials):
    ops = C.ops[:]
    i = rng.choice(ccx_idx)
    ops[4 * i] = CX
    ops[4 * i + 2] = ops[4 * i + 3]
    if check(gen, ops, cases):
        caught += 1
print(json.dumps({"q": curve.q, "addend": list(A), "mutation": "one CCX demoted to CX", "trials": trials, "caught": caught}))
