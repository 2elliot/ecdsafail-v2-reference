"""Gate check of the dialog multiplier with Euclid steps borrowing ring wires (point_add.RING_BORROW).

    python3.12 -I -B scripts/ring_borrow_check.py <tree> [broken]

For each q and setting, every x in [1, q) against a few y, forward and inverse, through
model_check.gate_pairs' simulation. U6_RING = 0 makes most steps want more ring than the tape
has, so the borrow path runs at these small q; a setting counts only if its recorded step shapes differ from
RING_BORROW = 0's. With "broken", each borrowed wire is flipped to 1 before the step; the check
must fail.
"""
import sys
from pathlib import Path
ROOT = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import io
import linecache
from qec import flatten, point_add, sim
from qec.builder import Circuit
from qec.arith import Mod
from qec.emit import write_circuit
from check import expand, structural


def gate_mul(q, ys, inverse=False):
    """model_check.gate_pairs over every x in [1, q): bad output or dirty-scratch lane bits."""
    pairs = [(x, y) for x in range(1, q) for y in ys]
    C = Circuit(hierarchical=True)
    mod = Mod(q)
    X, Y = C.alloc(mod.n), C.alloc(mod.n)
    inputs = (list(X), list(Y))
    point_add.DialogMultiplier(C, mod).apply(X, Y, inverse=inverse)
    buf = io.StringIO()
    write_circuit(buf, C, {"control": [], "x": X, "y": Y})
    ops = expand(flatten.load(io.StringIO(buf.getvalue())))
    structural(ops, C.num_wires)
    xs, values = zip(*pairs)
    state = sim.run(ops, C.num_wires, sim.pack(xs, inputs[0]) | sim.pack(values, inputs[1]), len(pairs))
    expected = [y * (pow(x, -1, q) if inverse else x) % q for x, y in pairs]
    want = sim.pack(xs, X) | sim.pack(expected, Y)
    bad = 0
    for w, actual in enumerate(state):
        bad |= actual ^ want.get(w, 0)
    return bad.bit_count()


if len(sys.argv) > 2 and sys.argv[2] == "broken":
    orig_alloc = Circuit.alloc

    def alloc(self, k=1):
        out = orig_alloc(self, k)
        frame = sys._getframe(1)
        if "borrow = (self.C.alloc" in linecache.getline(frame.f_code.co_filename, frame.f_lineno):
            for w in out:
                self.x(w)
        return out
    Circuit.alloc = alloc


def shapes(q):
    C = Circuit(hierarchical=True)
    mod = Mod(q)
    X, Y = C.alloc(mod.n), C.alloc(mod.n)
    m = point_add.DialogMultiplier(C, mod)
    m.apply(X, Y)
    return set(m._euclid_steps), C.num_wires


bad = 0
for q in (1021, 4093, 16381):
    for u6, borrow in ((193, 18), (0, 8), (0, 40)):
        point_add.U6_RING = u6
        point_add.RING_BORROW = 0
        base, w0 = shapes(q)
        point_add.RING_BORROW = borrow
        got, w1 = shapes(q)
        ys = [0, 1, q - 1, q // 3]
        b = gate_mul(q, ys) + gate_mul(q, ys, inverse=True)
        print(f"q={q} U6_RING={u6} RING_BORROW={borrow} borrowed={got != base} wires {w0}->{w1} bad={b}", flush=True)
        bad += b
point_add.U6_RING, point_add.RING_BORROW = 193, 0
print("OK" if bad == 0 else "FAIL")
sys.exit(0 if bad == 0 else 1)
