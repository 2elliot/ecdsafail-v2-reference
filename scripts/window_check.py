"""The high word's top borrow written straight onto its target (TOP_BORROW_DIRECT) and the lazy
window as one cell (LAZY_CELL), as gates.

    python -I -B scripts/window_check.py

A whole high word writes its top borrow E onto g1 with the ladder's last CCX. A lazy window's ladder
writes it onto p, the pending borrow, and then takes g1 ^= (p ^ E) s, which is E s because p s = 0
when the tick reaches its high word: a pending p means the word has a bit above the window. With
LAZY_CELL the window is not a ladder at all but a fixed list of gates with 6 CCX, the same map on
the window's domain:
  U's two bits or V's are 0; z is 1 beside bits of U's, 0 beside bits of V's;
  zg = g0 z; g1 = 0; cw, g0 and s are free while p = 0;
  a set p comes with s = 0, z = 0, U's bits 0 and V's bits both 1 (cw and g0 free).

1. The cell against the ladder (both switches off) on every point of that domain, 68 of them:
   every wire must agree. The cell must have 6 CCX and touch no scratch wire, and every one of its
   gates must be needed: dropping any one gives a wrong point.
2. The domain is what the gates see. At the primes below marked with a count, M is run as gates on
   every x with a marker in front of each lazy window, and every tuple the window finds on the
   forward pass must be in the domain. The run counts the tuples with p set.
3. For each (prime, LAZY_END) below M and its inverse must give y*x and y/x, with x kept and every
   other wire 0: with the cell, with the direct ladder (LAZY_CELL off), and with both off. The lanes
   are every x in 1..q-1 at the primes up to 17 bits, and lag_check's sampled inputs at the larger
   ones, each against y in 1, 2, q - 1 and one more. A pending borrow is rare: no x sets p at the
   primes just under a power of two, and hundreds do at the ones just above, so six of those are
   here; the run counts, from model_check's arithmetic model of p, the inputs that set p and those
   that reach a later high word with it set.
   Controls, which must give bad lanes on every build that has the gate:
     cell k     the cell without its k-th CCX, k = 1 .. 6: all six where some x sets p, five elsewhere
                (one of them only acts where a borrow stays pending)
     no write   LAZY_CELL off: the lazy ladder's last CCX (onto p) is dropped
     no s       LAZY_CELL off: g1 takes p ^ E whole, a CX in place of the CCX with s (builds where p is set)
     no g1      a whole word's last CCX (onto g1) is dropped
4. secp256k1 with the plan as built, on 500 of lag_check's inputs: 0 bad lanes, and the cell without
   its first CCX must fail.
"""
import argparse
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lag_check import inputs
from pair_check import bad_lanes
from qec.arith import Mod
from qec.builder import CCX, CX, SWAP, X, Circuit
from qec.model import tick_pass, tick_tape
from qec.point_add import DialogMultiplier

# model_check.py runs its checks at import, so its functions are loaded from the text above them.
_source = (Path(__file__).resolve().parent / "model_check.py").read_text()
_model = {"__name__": "model_check_functions", "__file__": str(Path(__file__).resolve().parent / "model_check.py")}
exec(compile(_source[:_source.index("\nrows = []")], "model_check.py", "exec"), _model)
_lazy_pending = _model["_lazy_pending"]

# (prime, LAZY_END values, probe the window's tuples). At a prime just above a power of two a borrow
# does leave the window and stay pending.
CASES = [(1021, (None,), False), (4093, (None, 4), False), (4099, (None, 4), True), (8209, (None, 4), True),
         (16411, (None,), True), (32771, (None, 5), True), (65537, (None,), True), (65539, (None, 5), True),
         (1048573, (4,), False), (4294967291, (6,), False), (2**61 - 1, (None,), False)]
SECP = 2**256 - 2**32 - 977
NAMES = ("U0", "U1", "V0", "V1", "cw", "g0", "z", "zg", "s", "p")


def domain():
    """The window's inputs as (u, v, z, cw, g0, s, p), u and v the two bits of U and of V."""
    points = []
    for u, v in itertools.product(range(4), repeat=2):
        if u and v:
            continue
        for z in ((1,) if u else (0,) if v else (0, 1)):
            for cw, g0, s, p in itertools.product((0, 1), repeat=4):
                if p and (s or z or u or v != 3):
                    continue
                points.append((u, v, z, cw, g0, s, p))
    return points


def window(cell, direct):
    """The lazy window alone: (ops, wires by name)."""
    D = DialogMultiplier
    saved = D.LAZY_CELL, D.TOP_BORROW_DIRECT
    D.LAZY_CELL, D.TOP_BORROW_DIRECT = cell, direct
    try:
        C = Circuit()
        m = D(C, Mod(SECP))
        w = dict(zip(("U0", "U1", "V0", "V1", "cw", "g0", "g1", "z", "zg", "E0", "E1", "E2", "s", "p"), C.alloc(14)))
        m._high_word([w["U0"], w["U1"]], [w["V0"], w["V1"]], w["cw"], w["g0"], w["g1"], (w["z"], w["zg"]),
                     [w["E0"], w["E1"], w["E2"]], (w["s"], w["p"]))
    finally:
        D.LAZY_CELL, D.TOP_BORROW_DIRECT = saved
    return list(C.prims), w


def run(ops, state):
    for j in range(0, len(ops), 4):
        op, a, b, c = ops[j:j + 4]
        if op == CCX:
            state[c] ^= state[a] & state[b]
        elif op == CX:
            state[b] ^= state[a]
        elif op == X:
            state[a] ^= 1
        else:
            state[a], state[b] = state[b], state[a]
    return state


def cell_rows():
    """(points, wrong points of the cell, of the direct ladder, CCX of the cell, scratch wires it
    touches, its gates, the gates whose removal no point catches)."""
    ladder, w = window(False, False)
    direct, _ = window(False, True)
    cell, _ = window(True, True)
    starts = []
    for u, v, z, cw, g0, s, p in domain():
        state = [0] * 14
        for name, bit in (("U0", u & 1), ("U1", u >> 1), ("V0", 1 ^ v & 1), ("V1", 1 ^ v >> 1), ("cw", cw), ("g0", g0),
                          ("z", z), ("zg", g0 & z), ("s", s), ("p", p)):
            state[w[name]] = bit
        starts.append(state)
    want = [run(ladder, list(state)) for state in starts]

    def wrong(ops):
        return sum(run(ops, list(state)) != out for state, out in zip(starts, want))

    scratch = {w["E0"], w["E1"], w["E2"]}
    touched = sum(1 for j in range(0, len(cell), 4) if scratch & set(cell[j + 1:j + 4]))
    spare = sum(1 for j in range(0, len(cell), 4) if wrong(cell[:j] + cell[j + 4:]) == 0)
    return len(starts), wrong(cell), wrong(direct), cell[0::4].count(CCX), touched, len(cell) // 4, spare


def seen(q):
    """(tuples the lazy windows find on the forward pass over every x, those outside the domain, those
    with p set). A marker, five CCX pairs that cancel, names the window's wires in the flat circuit."""
    D = DialogMultiplier
    real = D._high_word

    def marked(self, U, V, cw, g0, g1, held, scratch, lazy=None):
        if lazy:
            ws = [U[0], U[1], V[0], V[1], cw, g0, held[0], held[1], lazy[0], lazy[1]]
            for i in range(0, 10, 2):
                for _ in range(2):
                    self.C.prims.extend((CCX, ws[i], ws[i + 1], scratch[0]))
        real(self, U, V, cw, g0, g1, held, scratch, lazy)

    D._high_word = marked
    try:
        C = Circuit()
        mod = Mod(q)
        Xr, Yr = C.alloc(mod.n), C.alloc(mod.n)
        x_in, y_in = list(Xr), list(Yr)
        m = D(C, mod)
        lazy_ticks = sum(1 for key in m._tick_plan() if key and key[9])
        m.apply(Xr, Yr)
    finally:
        D._high_word = real
    ops, lanes = C.prims, q - 1
    full = (1 << lanes) - 1
    val = [0] * C.num_wires
    for i in range(mod.n):
        val[x_in[i]] = sum(1 << (x - 1) for x in range(1, q) if x >> i & 1)
    val[y_in[0]] = full
    good = {(u & 1, u >> 1, 1 ^ v & 1, 1 ^ v >> 1, cw, g0, z, g0 & z, s, p) for u, v, z, cw, g0, s, p in domain()}
    tuples = outside = pending = marks = 0
    j = 0
    while j < len(ops):
        op, a, b, c = ops[j:j + 4]
        if (op == CCX and marks < lazy_ticks and j + 40 <= len(ops)
                and all(ops[j + 8 * k] == CCX and ops[j + 8 * k + 3] == c
                        and ops[j + 8 * k:j + 8 * k + 4] == ops[j + 8 * k + 4:j + 8 * k + 8] for k in range(5))):
            # Forward markers only: the inverse pass replays them reversed, on the same states.
            marks += 1
            vs = [val[w] for k in range(5) for w in ops[j + 8 * k + 1:j + 8 * k + 3]]
            for lane in range(lanes):
                found = tuple(v >> lane & 1 for v in vs)
                tuples += 1
                outside += found not in good
                pending += found[9]
            j += 40
            continue
        if op == CCX:
            val[c] ^= val[a] & val[b]
        elif op == CX:
            val[b] ^= val[a]
        elif op == X:
            val[a] ^= full
        elif op == SWAP:
            val[a], val[b] = val[b], val[a]
        j += 4
    assert marks == lazy_ticks and all(val[Yr[i]] == val[Xr[i]] for i in range(mod.n)), q
    return tuples, outside, pending


def pending(q, xs):
    """(lazy ticks, whole high words, inputs that set p, inputs that reach a lazy tick's high word with
    p set)."""
    m = DialogMultiplier(Circuit(True), Mod(q))
    plan = m._tick_plan()
    lazy = sum(1 for key in plan if key and key[9])
    whole = sum(1 for key in plan if key and key[6] and not key[9])
    sets = hit = 0
    if lazy:
        for x in xs:
            states, cells = tick_pass(q, x, m.iters)[0], tick_tape(q, x, m.iters)[2]
            p = _lazy_pending(m, plan, states, cells)
            sets += any(after for _, _, after in p.values())
            # A subtract's high word finds p as the raw bit left it; a strip's finds it as the tick leaves it.
            hit += any((mid if cells[t][0] else after) for t, (_, mid, after) in p.items())
    return lazy, whole, sets, hit


def main():
    argparse.ArgumentParser().parse_args()
    D = DialogMultiplier
    saved = D.LAZY_END, D.TOP_BORROW_DIRECT, D.LAZY_CELL
    real = D._high_word

    def edited(edit):
        def high_word(self, U, V, cw, g0, g1, held, scratch, lazy=None):
            at = len(self.C.prims)
            real(self, U, V, cw, g0, g1, held, scratch, lazy)
            ops = self.C.prims[at:]
            del self.C.prims[at:]
            count = 0
            for j in range(0, len(ops), 4):
                op = tuple(ops[j:j + 4])
                count += op[0] == CCX
                self.C.prims.extend(edit(op, g1, lazy, count))
        return high_word

    def cell_without(k):
        return lambda op, g1, lazy, count: () if lazy and op[0] == CCX and count == k else op

    def no_write(op, g1, lazy, count):
        return () if lazy and op[0] == CCX and op[3] == lazy[1] else op

    def no_s(op, g1, lazy, count):
        return (CX, lazy[1], g1, g1) if lazy and op[0] == CCX and op[3] == g1 else op

    def no_g1(op, g1, lazy, count):
        return () if not lazy and op[0] == CCX and op[3] == g1 else op

    def control(edit, q, ys, xs):
        D._high_word = edited(edit)
        try:
            return bad_lanes(q, ys, False, xs)
        finally:
            D._high_word = real

    points, cell_wrong, direct_wrong, cell_ccx, touched, gates, spare = cell_rows()
    print(f"the cell against the ladder: {points} points, {cell_wrong} wrong; the direct ladder: {direct_wrong} wrong; "
          f"the cell has {cell_ccx} CCX in {gates} gates, {touched} on scratch, {spare} that no point needs", flush=True)
    ok = (points, cell_wrong, direct_wrong, cell_ccx, touched, spare) == (68, 0, 0, 6, 0, 0)

    bad = off_bad = builds = lazy_builds = whole_builds = set_builds = sets = hits = 0
    tuples = outside = set_tuples = 0
    fails = {"cell": 0, "no write": 0, "no s": 0, "no g1": 0}
    try:
        for q, ends, probe in CASES:
            xs = inputs(q)
            ys = sorted({1, 2, q - 1, q // 3})
            for D.LAZY_END in ends:
                D.TOP_BORROW_DIRECT, D.LAZY_CELL = True, True
                lazy, whole, some, hit = pending(q, xs)
                assert lazy
                builds += 1
                lazy_builds += 1
                whole_builds += whole > 0
                set_builds += some > 0
                sets += some
                hits += hit
                b = bad_lanes(q, ys, False, xs) + bad_lanes(q, ys, True, xs)
                bad += b
                row = []
                if probe:
                    found = seen(q)
                    tuples, outside, set_tuples = tuples + found[0], outside + found[1], set_tuples + found[2]
                    row.append(f"{found[0]} tuples, {found[1]} outside the domain, {found[2]} with p set")
                cell = [control(cell_without(k), q, ys, xs) for k in range(1, 7)]
                # One of the cell's CCX only acts where a borrow stays pending.
                needed = sum(n > 0 for n in cell)
                fails["cell"] += needed == (6 if some else 5)
                row.append(f"cell {needed} of 6")
                if whole:
                    n = control(no_g1, q, ys, xs)
                    fails["no g1"] += n > 0
                    row.append(f"no g1 {n}")
                D.LAZY_CELL = False
                off_bad += bad_lanes(q, ys, False, xs) + bad_lanes(q, ys, True, xs)
                for name, edit, there in (("no write", no_write, True), ("no s", no_s, some)):
                    if there:
                        n = control(edit, q, ys, xs)
                        fails[name] += n > 0
                        row.append(f"{name} {n}")
                D.TOP_BORROW_DIRECT = False
                off_bad += bad_lanes(q, ys, False, xs) + bad_lanes(q, ys, True, xs)
                print(f"q = {q}, LAZY_END {D.LAZY_END}: {lazy} lazy ticks, {whole} whole high words, {len(xs)} x, "
                      f"{some} set p and {hit} reach a high word with it set, {b} bad lanes; {'; '.join(row)}",
                      flush=True)
        D.LAZY_END, D.TOP_BORROW_DIRECT, D.LAZY_CELL = saved
        xs = inputs(SECP)
        xs = xs[::max(1, len(xs) // 500)]
        b = bad_lanes(SECP, [SECP // 3], False, xs) + bad_lanes(SECP, [SECP // 3], True, xs)
        secp = control(cell_without(1), SECP, [SECP // 3], xs)
        print(f"secp256k1, the built plan: {len(xs)} x, {b} bad lanes; the cell without its first CCX: {secp} bad lanes",
              flush=True)
        bad += b
    finally:
        D.LAZY_END, D.TOP_BORROW_DIRECT, D.LAZY_CELL = saved
        D._high_word = real
    print(f"{builds} builds ({whole_builds} with whole high words, {set_builds} where p is set), {sets} inputs that set "
          f"p, {hits} that reach a high word with it set, {bad} bad lanes; switches off: {off_bad} bad lanes; "
          f"{tuples} tuples seen, {outside} outside the domain, {set_tuples} with p set; controls fail: cell "
          f"{fails['cell']} of {lazy_builds}, no write {fails['no write']} of {lazy_builds}, no s {fails['no s']} of "
          f"{set_builds}, no g1 {fails['no g1']} of {whole_builds}")
    ok = (ok and whole_builds > 0 and set_builds >= 6 and hits > 100 and bad == 0 and off_bad == 0
          and tuples > 0 and outside == 0 and set_tuples > 100
          and fails["cell"] == lazy_builds and fails["no write"] == lazy_builds and fails["no s"] == set_builds
          and fails["no g1"] == whole_builds and secp > 0)
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
