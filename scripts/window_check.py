"""The high word's top borrow written straight onto its target (TOP_BORROW_DIRECT), as gates.

    python -I -B scripts/window_check.py

A whole high word writes its top borrow E onto g1 with the ladder's last CCX. A lazy window writes
it onto p, the pending borrow, and then takes g1 ^= (p ^ E) s, which is E s because p s = 0 when
the tick reaches its high word: a pending p means the word has a bit above the window. If p s were
ever 1, g1 would be flipped at that tick, so the check is M itself.

For each (prime, LAZY_END) below M and its inverse must give y*x and y/x, with x kept and every other
wire 0, with TOP_BORROW_DIRECT on and off. The lanes are every x in 1..q-1 at the primes up to 16
bits, and lag_check's sampled inputs at the larger ones, each against y in 1, 2, q - 1 and one more.
The run counts, from model_check's arithmetic model of the pending borrow, the inputs that reach a
lazy tick's high word with p set, so the builds do exercise it: p is set for no x at all at the
primes just under a power of two, and for hundreds at the ones just above, so six of those are here. Controls, which must give bad lanes
on every build that has the gate:
  no write   the lazy ladder's last CCX (onto p) is dropped
  no s       g1 takes p ^ E whole, a CX in place of the CCX with s (on the builds where some x sets p)
  no g1      a whole word's last CCX (onto g1) is dropped
Last, secp256k1 with the plan as built, on 500 of lag_check's inputs, with the first control, which must
fail, and the third, which is only reported: those inputs end before the whole high words, the last 24 ticks.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lag_check import inputs
from pair_check import bad_lanes
from qec.arith import Mod
from qec.builder import CCX, CX, Circuit
from qec.model import tick_pass, tick_tape
from qec.point_add import DialogMultiplier


# model_check.py runs its checks at import, so its functions are loaded from the text above them.
_source = (Path(__file__).resolve().parent / "model_check.py").read_text()
_model = {"__name__": "model_check_functions", "__file__": str(Path(__file__).resolve().parent / "model_check.py")}
exec(compile(_source[:_source.index("\nrows = []")], "model_check.py", "exec"), _model)
_lazy_pending = _model["_lazy_pending"]


# At a prime just above a power of two a borrow does leave the window and stay pending: 30, 40, 87,
# 253, 516 and 510 of the x below set p, and 2, 4, 8, 47, 107 and 102 reach a later high word with it.
CASES = [(1021, (None,)), (4093, (None, 4)), (4099, (None, 4)), (8209, (None, 4)), (16411, (None,)),
         (32771, (None, 5)), (65537, (None,)), (65539, (None, 5)), (1048573, (4,)), (4294967291, (6,)),
         (2**61 - 1, (None,))]


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
            # A subtract's high word finds p as the raw bit left it; a strip's finds it as the tick leaves it.
            sets += any(after for _, _, after in p.values())
            hit += any((mid if cells[t][0] else after) for t, (_, mid, after) in p.items())
    return lazy, whole, sets, hit


def main():
    argparse.ArgumentParser().parse_args()
    D = DialogMultiplier
    saved = D.LAZY_END, D.TOP_BORROW_DIRECT
    real = D._high_word

    def edited(edit):
        def high_word(self, U, V, cw, g0, g1, held, scratch, lazy=None):
            at = len(self.C.prims)
            real(self, U, V, cw, g0, g1, held, scratch, lazy)
            ops = self.C.prims[at:]
            del self.C.prims[at:]
            for j in range(0, len(ops), 4):
                self.C.prims.extend(edit(tuple(ops[j:j + 4]), g1, lazy))
        return high_word

    def no_write(op, g1, lazy):
        return () if lazy and op[0] == CCX and op[3] == lazy[1] else op

    def no_s(op, g1, lazy):
        return (CX, lazy[1], g1, g1) if lazy and op[0] == CCX and op[3] == g1 else op

    def no_g1(op, g1, lazy):
        return () if not lazy and op[0] == CCX and op[3] == g1 else op

    bad = off_bad = builds = lazy_builds = whole_builds = set_builds = sets = hits = 0
    fails = {"no write": 0, "no s": 0, "no g1": 0}
    try:
        for q, ends in CASES:
            xs = inputs(q)
            ys = sorted({1, 2, q - 1, q // 3})
            for D.LAZY_END in ends:
                D.TOP_BORROW_DIRECT = True
                lazy, whole, some, hit = pending(q, xs)
                if not (lazy or whole):
                    continue
                builds += 1
                lazy_builds += lazy > 0
                whole_builds += whole > 0
                set_builds += some > 0
                sets += some
                hits += hit
                b = bad_lanes(q, ys, False, xs) + bad_lanes(q, ys, True, xs)
                bad += b
                row = []
                for name, edit, there in (("no write", no_write, lazy), ("no s", no_s, some), ("no g1", no_g1, whole)):
                    if there:
                        D._high_word = edited(edit)
                        n = bad_lanes(q, ys, False, xs)
                        D._high_word = real
                        fails[name] += n > 0
                        row.append(f"{name} {n}")
                D.TOP_BORROW_DIRECT = False
                off_bad += bad_lanes(q, ys, False, xs) + bad_lanes(q, ys, True, xs)
                print(f"q = {q}, LAZY_END {D.LAZY_END}: {lazy} lazy ticks, {whole} whole high words, {len(xs)} x, "
                      f"{some} set p and {hit} reach a high word with it set, {b} bad lanes; controls: {', '.join(row)}",
                      flush=True)
        D.LAZY_END, D.TOP_BORROW_DIRECT = saved
        q = 2**256 - 2**32 - 977
        xs = inputs(q)
        xs = xs[::max(1, len(xs) // 500)]
        b = bad_lanes(q, [q // 3], False, xs) + bad_lanes(q, [q // 3], True, xs)
        secp = []
        for edit in (no_write, no_g1):
            D._high_word = edited(edit)
            secp.append(bad_lanes(q, [q // 3], False, xs))
            D._high_word = real
        print(f"secp256k1, the built plan: {len(xs)} x, {b} bad lanes; no write: {secp[0]} bad lanes, "
              f"no g1: {secp[1]}", flush=True)
        bad += b
    finally:
        D.LAZY_END, D.TOP_BORROW_DIRECT = saved
        D._high_word = real
    print(f"{builds} builds ({lazy_builds} with lazy ticks, {whole_builds} with whole high words, {set_builds} "
          f"where p is set), {sets} inputs that set p, {hits} that reach a high word with it set, {bad} bad "
          f"lanes; TOP_BORROW_DIRECT off: {off_bad} bad lanes; controls fail: no write {fails['no write']} of "
          f"{lazy_builds}, no s {fails['no s']} of {set_builds}, no g1 {fails['no g1']} of {whole_builds}")
    ok = (lazy_builds > 0 and whole_builds > 0 and set_builds >= 6 and hits > 100 and bad == 0 and off_bad == 0
          and fails["no write"] == lazy_builds and fails["no s"] == set_builds and fails["no g1"] == whole_builds
          and secp[0] > 0)
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
