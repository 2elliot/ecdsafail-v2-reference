"""Exhaustive gate-level check of the dialog packer, qec.point_add.pack_trits.

For every group length k up to 7 and every valid dialog of k (b0, sw) pairs
(sw implies b0, so 3^k dialogs), it simulates pack_trits and requires the
packed register to hold model.pack_value(dialog) with every other wire zero,
then simulates unpack_trits after it and requires the original pairs back with
every other wire zero. Prints one JSON line per k and exits 1 on any failure.
"""
import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec import sim
from qec.builder import Circuit
from qec.model import TRIT_GROUP, pack_value, trit_width
from qec.point_add import pack_trits, unpack_trits

PAIRS = [(0, 0), (1, 0), (1, 1)]


def check(k):
    C = Circuit()
    trits = [C.alloc(2) for _ in range(k)]
    inputs = [w for g in trits for w in g]
    packed = pack_trits(C, [list(g) for g in trits])
    n_pack = len(C)
    back = unpack_trits(C, list(packed), k)
    outputs = [w for g in back for w in g]
    dialogs = list(itertools.product(PAIRS, repeat=k))
    lanes = len(dialogs)
    init = [sum((b0 | sw << 1) << 2 * i for i, (b0, sw) in enumerate(d)) for d in dialogs]
    raw = trit_width(k) >= 2 * k
    want = init if raw else [pack_value(d) for d in dialogs]
    failures = 0

    s = sim.run(C.ops[:4 * n_pack], C.num_wires, sim.pack(init, inputs), lanes)
    got = sim.unpack(s, packed, lanes)
    failures += sum(g != w for g, w in zip(got, want))
    dirty = [w for w in range(C.num_wires) if w not in set(packed) and s[w]]

    s = sim.run(C.ops, C.num_wires, sim.pack(init, inputs), lanes)
    failures += sum(g != w for g, w in zip(sim.unpack(s, outputs, lanes), init))
    dirty += [w for w in range(C.num_wires) if w not in set(outputs) and s[w]]

    ops = C.ops[:4 * n_pack]
    return {
        "k": k, "dialogs": lanes, "raw_wires": 2 * k, "packed_wires": len(packed),
        "pack_ccx": ops[0::4].count(2), "pack_ops": n_pack, "peak_wires": C.num_wires,
        "failures": failures, "dirty_wires": len(dirty), "group_used": k == TRIT_GROUP,
    }


def main():
    ok = True
    for k in range(1, 8):
        row = check(k)
        ok &= row["failures"] == 0 and row["dirty_wires"] == 0
        print(json.dumps(row))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
