"""The cleared strips of lag blocks (_lag_clear), as gates.

    python -I -B scripts/lag_check.py

For each (prime, forced block length c) below (LAG_BLOCK = (c,), LAG_MIN_SAVING None, as
schedule_check.py forces them), whose plans have lag ticks that keep their threshold: M and its
inverse must give y*x and y/x, with x kept and every other wire 0, both with LAG_CLEAR on and off.
The lanes are every x in 1..q-1 at the primes up to 16 bits, and at the larger ones powers of two
and their neighbours, values near q, q/k and 3,000 seeded random x; each against y in 1, 2, q - 1
and one more. The run counts the cleared strips and the inputs whose strip was set
(two AND NOT flag) at such a tick. Controls: with the clearing CCX dropped, and with the flag's
sources left out of it (two in place of two AND NOT flag), the same builds must give bad lanes.
Then the wires the lazy run borrows: at q = 1021, 4093 and 8191, every x, with LAZY_END None, 4 and
5, and on the sampled inputs at 20 and 32 bits, where the run keeps fewer counter wires than the
counter has, the same must hold, and keeping one counter wire fewer than _lazy_kept says must give bad lanes.
Last, secp256k1 with the plan as built, on 500 of those inputs, with the dropped-CCX control.
"""
import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pair_check import bad_lanes
from qec.arith import Mod
from qec.builder import CCX, Circuit
from qec.model import tick_pass, tick_tape
from qec.point_add import DialogMultiplier


CASES = [(4093, 1), (8191, 1), (65521, 1), (65521, 2), (4194301, 3), (268435399, 3), (1099511627689, 5),
         (2**61 - 1, 7)]


def inputs(q):
    if q < 1 << 17:
        return list(range(1, q))
    rng = random.Random(q)
    xs = {1, 2, 3, q - 1, q - 2, (q - 1) // 2, (q + 1) // 2}
    for j in range(1, q.bit_length()):
        xs |= {1 << j, (1 << j) - 1, (1 << j) + 1, q - (1 << j), q >> j, (q >> j) | 1}
    for k in range(2, 40):
        xs |= {q // k, q // k + 1, 2 * (q // k) + 1}
    xs |= {rng.randrange(1, q) for _ in range(3000)}
    return sorted(x for x in xs if 0 < x < q)


def cleared(q, xs):
    """(ticks whose strip _lag_clear clears, inputs x with a set strip at one of them)."""
    m = DialogMultiplier(Circuit(True), Mod(q))
    plan = m._tick_plan()
    ticks = [t for t, key in enumerate(plan) if key and key[10] and key[7] == "same"]
    hit = 0
    for x in xs:
        states, cells = tick_pass(q, x, m.iters)[0], tick_tape(q, x, m.iters)[2]
        hit += any(not cells[t][0] and cells[t][1] and not states[t][0] >> plan[t][5] for t in ticks)
    return len(ticks), hit


def main():
    argparse.ArgumentParser().parse_args()
    D = DialogMultiplier
    saved = D.LAG_BLOCK, D.LAG_MIN_SAVING, D.LAG_CLEAR
    real = D._lag_clear
    bad = off_bad = builds = strips = hits = dropped = blind = 0

    def no_ccx(self, after, flag, borrowed):
        at = len(self.C.prims)
        real(self, after, flag, borrowed)
        ops = self.C.prims[at:]
        keep = [j for j in range(0, len(ops), 4) if ops[j] != CCX]
        del self.C.prims[at:]
        for j in keep:
            self.C.prims.extend(ops[j:j + 4])

    def no_flag(self, after, flag, borrowed):
        # NOT g0 alone in place of NOT (g0 ^ flag left): two, not two AND NOT flag.
        for e, wire, (g0, b) in borrowed:
            self.C.x(g0)
            self.C.ccx(b, g0, wire)
            self.C.x(g0)
            self.C.free([wire])

    try:
        for q, c in CASES:
            if True:
                xs = inputs(q)
                ys = sorted({1, 2, q - 1, q // 3})
                D.LAG_BLOCK, D.LAG_MIN_SAVING, D.LAG_CLEAR = (c,), None, True
                ticks, hit = cleared(q, xs)
                assert ticks, (q, c)
                builds += 1
                strips += ticks
                hits += hit
                b = bad_lanes(q, ys, False, xs) + bad_lanes(q, ys, True, xs)
                bad += b
                print(f"q = {q}, blocks of {c}: {ticks} cleared strips, {len(xs)} x, {hit} with one set, {b} bad lanes", flush=True)
                D._lag_clear = no_ccx
                dropped += bad_lanes(q, ys, False, xs) > 0
                D._lag_clear = no_flag
                blind += bad_lanes(q, ys, False, xs) > 0
                D._lag_clear = real
                D.LAG_CLEAR = False
                off_bad += bad_lanes(q, ys, False, xs) + bad_lanes(q, ys, True, xs)
        # The lazy run's borrowed wires (V's extra wire, the counter above _lazy_kept): toy primes,
        # every x, with the run ending at 4 and 5 bits (the counter stays) and at 3 (it is cleared).
        D.LAG_BLOCK, D.LAG_MIN_SAVING, D.LAG_CLEAR = saved
        end, kept_real = D.LAZY_END, D._lazy_kept
        runs = long_runs = lazy_bad = lazy_control = short = 0
        for q, ends in ((1021, (None, 4, 5)), (4093, (None, 4, 5)), (8191, (None, 4, 5)), (1048573, (4,)),
                        (4294967291, (6,))):
            ys = sorted({1, 2, q - 1, q // 3})
            xs = inputs(q)
            for D.LAZY_END in ends:
                probe = D(Circuit(True), Mod(q))
                plan = probe._tick_plan()
                if not any(key and key[9] for key in plan):
                    continue
                runs += 1
                lazy_bad += bad_lanes(q, ys, False, xs) + bad_lanes(q, ys, True, xs)
                if probe._lazy_kept(plan):
                    # Control: one counter wire too few kept, so a wire that still counts is given back.
                    long_runs += 1
                    short += probe._lazy_kept(plan) < probe._lazy_width(plan)
                    D._lazy_kept = lambda self, plan: kept_real(self, plan) - 1
                    lazy_control += bad_lanes(q, ys, False, xs) > 0
                    D._lazy_kept = kept_real
        D.LAZY_END = end
        print(f"lazy runs: {runs} builds, {long_runs} that keep counter wires ({short} of them fewer than the "
              f"counter has), {lazy_bad} bad lanes; "
              f"control (one wire fewer kept) fails {lazy_control} of {long_runs}", flush=True)
        bad += lazy_bad
        dropped += 0 if runs and short and lazy_control == long_runs else -10 ** 6
        # secp256k1, the plan as built: 500 of the inputs above, one y each.
        D.LAG_BLOCK, D.LAG_MIN_SAVING, D.LAG_CLEAR = saved
        q = 2**256 - 2**32 - 977
        xs = inputs(q)
        xs = xs[::max(1, len(xs) // 500)]
        ticks, hit = cleared(q, xs)
        b = bad_lanes(q, [q // 3], False, xs) + bad_lanes(q, [q // 3], True, xs)
        D._lag_clear = no_ccx
        secp_control = bad_lanes(q, [q // 3], False, xs)
        D._lag_clear = real
        print(f"secp256k1, the built plan: {ticks} cleared strips, {len(xs)} x, {hit} with one set, {b} bad lanes; "
              f"no CCX: {secp_control} bad lanes", flush=True)
        bad += b
        dropped += 0 if secp_control else -10 ** 6
    finally:
        D.LAG_BLOCK, D.LAG_MIN_SAVING, D.LAG_CLEAR = saved
        D._lag_clear = real
    print(f"{builds} builds with a cleared strip ({strips} strips, {hits} inputs with one set), {bad} bad lanes; "
          f"LAG_CLEAR off: {off_bad} bad lanes; controls: no CCX fails {dropped} of {builds} builds, "
          f"no flag fails {blind} of {builds}")
    ok = builds == len(CASES) and hits > 0 and bad == 0 and off_bad == 0 and dropped == builds and blind > 0
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
