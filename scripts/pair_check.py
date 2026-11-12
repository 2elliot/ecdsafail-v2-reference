"""The two-bit junk cell of the Euclid chain (_pair_up, _pair_down) against the single-bit cells.

    python -I -B scripts/pair_check.py [--below 260]

Cell, every input: a_0, b_0, a_1, b_1, c and (g0, g1) in {00, 10, 11}, the junk wire 0.
  - after _pair_up the junk wire holds c_2 = MAJ(a_1, b_1, MAJ(a_0, b_0, c)), and c and g0 are as
    found; _pair_up never touches g1;
  - after _pair_down A_i = a_i ^ g1*(a_i ^ b_i ^ 1), B_i = b_i ^ g0*(a_i ^ c_i), c, g0 and g1 are
    as found: the values the in-place and single junk descents of _shared_difference leave;
  - the cell is 6 CCX, and reversing both halves restores the input with the junk wire 0.
  Mutants: the cell with any one gate dropped, and with any one CCX turned into a CX from either
  control, must fail on some input; every one does.
Multiplier, every odd prime below --below with n >= 4 and every budget of junk wires from 0 to
  the plan's full count: M and its inverse on every x in 1..q-1 and six y, as gates, must give
  y*x and y/x with x kept and every other wire 0. The same build with JUNK_PAIRS off must agree,
  and the run counts the builds that hold a pair, a single beside pairs, and a pair under an
  in-place bit. A control with _pair_down's last CCX dropped must give bad lanes.
"""
import argparse
import io
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from check import expand, structural
from qec import flatten, sim
from qec.arith import Mod
from qec.builder import CCX, CX, Circuit
from qec.emit import write_circuit
from qec.point_add import DialogMultiplier

NAMES = ["A0", "B0", "A1", "B1", "Cw", "J", "g0", "g1"]


def cell_ops():
    """(up ops, down ops) on wires 0..7 in NAMES' order."""
    C = Circuit()
    w = C.alloc(8)
    DialogMultiplier._pair_up(C, *w[:6])
    cut = len(C.ops)
    DialogMultiplier._pair_down(C, *w)
    return list(C.ops[:cut]), list(C.ops[cut:])


def run(ops, state):
    state = list(state)
    for j in range(0, len(ops), 4):
        op, a, b, c = ops[j:j + 4]
        if op == CCX:
            state[c] ^= state[a] & state[b]
        elif op == CX:
            state[b] ^= state[a]
        else:
            assert op == 0
            state[a] ^= 1
    return state


def cell_bad(up, down):
    """Inputs on which the cell breaks its contract."""
    bad = 0
    for a0, b0, a1, b1, c in itertools.product((0, 1), repeat=5):
        c1 = (a0 & b0) | (a0 & c) | (b0 & c)
        c2 = (a1 & b1) | (a1 & c1) | (b1 & c1)
        for g0, g1 in ((0, 0), (1, 0), (1, 1)):
            start = [a0, b0, a1, b1, c, 0, g0, 0]
            mid = run(up, start)
            ok = mid[5] == c2 and mid[4] == c and mid[6] == g0 and mid[7] == 0
            mid[7] = g1
            end = run(down, mid)
            want = [a0 ^ (g1 & (a0 ^ b0 ^ 1)), b0 ^ (g0 & (a0 ^ c)), a1 ^ (g1 & (a1 ^ b1 ^ 1)),
                    b1 ^ (g0 & (a1 ^ c1)), c, None, g0, g1]
            ok = ok and all(x is None or x == y for x, y in zip(want, end))
            # The reversed halves are the inverse pass's.
            back = run(reversed_ops(down), end)
            back[7] = 0
            back = run(reversed_ops(up), back)
            ok = ok and back == start
            bad += not ok
    return bad


def reversed_ops(ops):
    out = []
    for j in range(len(ops) - 4, -1, -4):
        out += ops[j:j + 4]
    return out


def cell_rows():
    up, down = cell_ops()
    ccx = (up + down)[0::4].count(CCX)
    bad = cell_bad(up, down)
    mutants = caught = 0
    both = up + down
    for j in range(0, len(both), 4):
        variants = [both[:j] + both[j + 4:]]
        if both[j] == CCX:
            _, a, b, c = both[j:j + 4]
            variants += [both[:j] + [CX, a, c, c] + both[j + 4:], both[:j] + [CX, b, c, c] + both[j + 4:]]
        for v in variants:
            cut = len(up) - (4 if j < len(up) and len(v) < len(both) else 0)
            mutants += 1
            caught += cell_bad(v[:cut], v[cut:]) > 0
    print(f"cell: {ccx} CCX, {bad} bad inputs of 96; mutants caught {caught} of {mutants}", flush=True)
    return bad == 0 and ccx == 6 and caught == mutants


def bad_lanes(q, ys, inverse):
    pairs = [(x, y) for x in range(1, q) for y in ys]
    C = Circuit(hierarchical=True)
    mod = Mod(q)
    X, Y = C.alloc(mod.n), C.alloc(mod.n)
    inputs = (list(X), list(Y))
    DialogMultiplier(C, mod).apply(X, Y, inverse=inverse)
    buf = io.StringIO()
    write_circuit(buf, C, {"control": [], "x": X, "y": Y})
    ops = expand(flatten.load(io.StringIO(buf.getvalue())))
    structural(ops, C.num_wires)
    xs, values = zip(*pairs)
    init = sim.pack(xs, inputs[0]) | sim.pack(values, inputs[1])
    state = sim.run(ops, C.num_wires, init, len(pairs))
    expected = [y * (pow(x, -1, q) if inverse else x) % q for x, y in pairs]
    want = sim.pack(xs, X) | sim.pack(expected, Y)
    bad = 0
    for w, actual in enumerate(state):
        bad |= actual ^ want.get(w, 0)
    return bad.bit_count()


def shapes(q):
    """(has a pair, has a single beside a pair, has a pair under an in-place bit) for the plan."""
    m = DialogMultiplier(Circuit(True), Mod(q))
    pair = mixed = under = False
    for key in m._tick_plan():
        if key is None or key[3] or not key[4]:
            continue
        w, top_clear, narrow, kappa, high = key[0], key[1], key[2], key[5], key[6]
        bits = DialogMultiplier._carry_junk(w, top_clear, narrow, kappa if high else None)
        chain = kappa if high else w - (1 if w >= 2 and (top_clear or narrow) else 0)
        singles, pairs = DialogMultiplier._junk_split(key[4], bits)
        pair |= pairs > 0
        mixed |= pairs > 0 and singles > 0
        under |= pairs > 0 and singles + 2 * pairs < chain - 1
    return pair, mixed, under


def multiplier_rows(below):
    primes = [p for p in range(11, below, 2) if all(p % d for d in range(3, int(p ** 0.5) + 1, 2))]
    saved = DialogMultiplier.CARRY_JUNK, DialogMultiplier.JUNK_PAIRS
    bad = builds = off_bad = 0
    seen = [0, 0, 0]
    try:
        for q in primes:
            ys = sorted({0, 1, 2, q - 1, q // 2, 3 % q})
            DialogMultiplier.CARRY_JUNK, DialogMultiplier.JUNK_PAIRS = 10 ** 9, True
            probe = DialogMultiplier(Circuit(True), Mod(q))
            full = sum(key[4] for key in probe._tick_plan() if key)
            for budget in range(full + 1):
                DialogMultiplier.CARRY_JUNK, DialogMultiplier.JUNK_PAIRS = budget, True
                b = bad_lanes(q, ys, False) + bad_lanes(q, ys, True)
                for k, flag in enumerate(shapes(q)):
                    seen[k] += flag
                bad += b
                builds += 1
                if b:
                    print("FAIL", q, budget, b, flush=True)
            DialogMultiplier.CARRY_JUNK, DialogMultiplier.JUNK_PAIRS = 10 ** 9, False
            off_bad += bad_lanes(q, ys, False) + bad_lanes(q, ys, True)
        # Control: the pair's last CCX dropped.
        DialogMultiplier.CARRY_JUNK, DialogMultiplier.JUNK_PAIRS = 10 ** 9, True
        real = DialogMultiplier._pair_down

        def broken(C, *wires):
            at = len(C.prims)
            real(C, *wires)
            ops = C.prims[at:]
            last = max(j for j in range(0, len(ops), 4) if ops[j] == CCX)
            del C.prims[at + last:at + last + 4]

        DialogMultiplier._pair_down = staticmethod(broken)
        control = sum(bad_lanes(q, [1, 2, q - 1], False) for q in primes[:6])
        DialogMultiplier._pair_down = staticmethod(real)
    finally:
        DialogMultiplier.CARRY_JUNK, DialogMultiplier.JUNK_PAIRS = saved
    print(f"multiplier: {len(primes)} primes below {below}, {builds} builds, {bad} bad lanes; builds with a pair "
          f"{seen[0]}, with a single beside pairs {seen[1]}, with a pair under an in-place bit {seen[2]}; "
          f"pairs off: {off_bad} bad lanes; control: {control} bad lanes", flush=True)
    return bad == 0 and off_bad == 0 and control > 0 and all(seen)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--below", type=int, default=260)
    a = ap.parse_args()
    ok = cell_rows()
    ok = multiplier_rows(a.below) and ok
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
