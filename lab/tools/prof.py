"""PROFILING ONLY. Emit a generator tree's secp256k1 circuit with one part made free (its gates
dropped while wires are still allocated as before). The emitted circuit is WRONG; its tick count
is only an upper bound on what that part can save. Never submit or simulate its output as a
candidate.

    python3.12 -I -B prof.py <tree> <part[,part...]|none> <out.json>
"""
import json
import sys
from pathlib import Path

tree, parts, out = Path(sys.argv[1]).resolve(), sys.argv[2], Path(sys.argv[3])
sys.path.insert(0, str(tree))

from qec import builder, two_bit, point_add, square, arith  # noqa: E402
from qec.curve import SECP256K1, SECP256K1_G  # noqa: E402
from qec.emit import write_circuit  # noqa: E402

MUTE = builder.MUTE


def muted(fn):
    def g(*a, **k):
        MUTE[0] += 1
        try:
            return fn(*a, **k)
        finally:
            MUTE[0] -= 1
    return g


def free_method(cls, name):
    setattr(cls, name, muted(getattr(cls, name)))


def free_func(mod, name):
    setattr(mod, name, muted(getattr(mod, name)))


DM, TB, CPA = point_add.DialogMultiplier, two_bit.TwoBitMod, point_add.ControlledPointAdd
PARTS = {
    # Euclid
    "euclid_step": lambda: free_method(DM, "_euclid_step"),          # every general step (not tick0)
    "euclid_diff": lambda: free_method(DM, "_shared_difference"),    # compare + difference in each step
    "euclid_rotate": lambda: free_func(point_add, "_par_rotate"),     # the conditional rotation
    "euclid_tick0": lambda: free_method(DM, "_euclid_tick0"),
    "euclid_pass": lambda: free_method(DM, "_euclid_pass"),          # steps + sweep + block codes (shared code blocks too)
    # Bezout
    "bezout_tick": lambda: free_method(TB, "tick"),
    "bz_add": lambda: free_func(two_bit, "add_or_shift_la"),
    "bz_compare": lambda: free_func(two_bit, "compare_blocks"),
    "bz_prefixP": lambda: (free_method(two_bit.PrefixR, "__init__"), free_method(two_bit.PrefixR, "write")),
    "bz_prefixM": lambda: (free_method(two_bit.Prefix, "__init__"), free_method(two_bit.Prefix, "write"),
                           free_method(two_bit.Prefix, "release"), free_method(two_bit.Prefix, "restore")),
    "bz_chains": lambda: (free_method(TB, "_chain_up"), free_method(TB, "_chain_down")),
    "bz_final_swap": lambda: two_bit.PROF_FREE.add("final_swap"),
    "bz_b1_bracket": lambda: two_bit.PROF_FREE.add("b1_bracket"),
    "bz_decode": lambda: (free_method(DM, "_decode"),),
    # outside the multiplier
    "square": lambda: free_method(CPA, "_square_mulsub"),
    "consts": lambda: (free_method(arith.Mod, "cmodadd_const"), free_method(arith.Mod, "creflect")),
    "eq_flags": lambda: free_method(CPA, "_eq_flags"),
}
if parts != "none":
    for p in parts.split(","):
        PARTS[p]()

gen = CPA(SECP256K1, SECP256K1_G)
C = gen.build(True)
with open(out, "w") as fh:
    written = write_circuit(fh, C, gen.layout)
print(json.dumps({"parts": parts, "num_wires": C.num_wires, **{k: v for k, v in written.items()}}))
