"""The high-word schedule against the joint bound, with no use of qec.model.

A flagged tick's kappa must be min(w, K_t), K_t the least K with 2^(3K + 2t) > 2 q^2 (q - 1), and a
high tick's must be K_t itself. On the semantic Euclid recurrence, with the held offset tracked
beside it, no state before a flagged tick has both registers at 2^kappa or more, and a subtract
there leaves U below 2^kappa.

The lazy ticks must be the high ticks with more than WINDOW bits above kappa, in one run that starts
at WINDOW + 1 bits after a drop at the same width and ends at WINDOW + 1 bits before a narrower tick
at the same threshold; raw must say that the tick before kept its threshold; the counter must be
wide enough. On the same recurrence, with the pending borrow p followed through the lazy rules, the
window is all ones under p once the tick's crossing bits have met it, p is never set beside U's
word, physical V fits w + 1 bits, and the last tick's top bit and extra wire are what its exit
assumes.

The lag blocks must be runs of c lazy ticks, none the lazy run's first or last, with fronts F, F - 1, ..
and the last one marked; F must keep the holes clear of the cells a lazy tick reads, F - 2 e >= cut_e + 2
(cut_e + 1 at the block's last tick, whose holes are closed before the next tick reads),
and above a U whose flag came on inside the block, F - (c - 1) - e >= kappa_(e-1) for e >= 1, leave
a cell at or above the front after the block, and stay inside V's wires at every tick, F - e <= w_e + 1 + c. On the recurrence, with k the strips counted so far
(two AND NOT flag): F - e - k >= cut + 2 at every tick but the last, cut + 1 there; under the flag U is below 2^(F - e - k) when
k > 0 and V is below 2^kappa.

    python3 -I -B scripts/schedule_check.py     # toy primes, every x; secp256k1, sampled x
"""
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from qec.arith import Mod
from qec.builder import Circuit
from qec.point_add import DialogMultiplier

SECP = 2**256 - 2**32 - 977
TOY = [3, 5, 7, 11, 13, 31, 43, 61, 127, 151, 251, 509, 1021, 4093, 8117, 65521]
WINDOW = 2


def joint_threshold(q, t):
    k = 0
    while 1 << (3 * k + 2 * t) <= 2 * q * q * (q - 1):
        k += 1
    return k


def plan_of(q, lag=None):
    """The plan, or with lag = c the plan that has a lag block of c ticks wherever one fits."""
    saved = DialogMultiplier.LAG_BLOCK, DialogMultiplier.LAG_MIN_SAVING
    if lag:
        DialogMultiplier.LAG_BLOCK, DialogMultiplier.LAG_MIN_SAVING = (lag,), None
    try:
        return DialogMultiplier(Circuit(), Mod(q))._tick_plan()
    finally:
        DialogMultiplier.LAG_BLOCK, DialogMultiplier.LAG_MIN_SAVING = saved


def thresholds_ok(q, plan):
    for t, key in enumerate(plan):
        if key is None or key[5] is None:
            if key is not None and key[6]:
                return False
            continue
        w, kappa, high = key[0], key[5], key[6]
        K = joint_threshold(q, t)
        if kappa != min(w, K) or (high and not 1 <= kappa == K < w):
            return False
    return True


def states_ok(q, plan, xs):
    for x in xs:
        u, v, e = q, x, 0
        for key in plan:
            lim = None if key is None or key[5] is None else 1 << key[5]
            V = v - e
            if V < 0 or (lim and u >= lim and (V ^ e) >= lim):
                return False
            if v & 1:
                e = int(u > v)
                u, v = min(u, v), abs(u - v) // 2
                if lim and u >= lim:
                    return False
            else:
                v //= 4 if v % 4 == 0 else 2
        if (u, v) != (1, 0):
            return False
    return True


def lazy_fields_ok(q, plan):
    lazy = [t for t, key in enumerate(plan) if key and key[9]]
    long_high = [t for t, key in enumerate(plan) if key and key[6] and key[0] - key[5] > WINDOW]
    if not lazy:
        return True
    first, final = lazy[0], lazy[-1]
    if lazy != long_high or lazy != list(range(first, final + 1)) or final + 1 >= len(plan) or first < 1:
        return False
    K = [None if key is None else min(key[0], joint_threshold(q, t)) for t, key in enumerate(plan)]
    w = [None if key is None else key[0] for key in plan]
    if w[first - 1] != w[first] or K[first - 1] != K[first] + 1 or w[first] - K[first] != WINDOW + 1:
        return False
    if w[final + 1] != w[final] - 1 or K[final + 1] != K[final] or w[final] - K[final] != WINDOW + 1:
        return False
    width = 0
    for t in lazy:
        raw, m_in, m_out, leave = plan[t][9]
        if raw != (t > first and K[t - 1] == K[t]) or (t > first and not raw and K[t - 1] != K[t] + 1):
            return False
        if m_in != width or leave != (t == final) or plan[t][5] != K[t]:
            return False
        if w[t] - K[t] - WINDOW + 1 > 1 << (m_out - 1) or m_out < 1:
            return False
        width = m_out
    return True


def lazy_states_ok(q, plan, xs):
    lazy = {t for t, key in enumerate(plan) if key and key[9]}
    if not lazy:
        return True
    first, final = min(lazy), max(lazy)
    mask = (1 << WINDOW) - 1
    for x in xs:
        u, v, e = q, x, 0
        p, phys = 0, None
        for t, key in enumerate(plan):
            V = v - e
            odd, two = v & 1, int(v % 4 == 0)
            if odd:
                un, vn, en = min(u, v), abs(u - v) // 2, int(u > v)
            else:
                un, vn, en = u, v // (4 if two else 2), e
            Vn = vn - en
            if t in lazy:
                w, k, (raw, _, _, leave) = key[0], key[5], key[9]
                cut = k + WINDOW
                if t == first:
                    phys, p = V, 0
                    if V >> (w - 1):
                        return False
                if raw:
                    r = phys >> (cut - 1) & 1
                    phys ^= p << (cut - 1)
                    p &= 1 - r
                if phys - (p << cut) != V or phys >> (w + 1):
                    return False
                z, hi = u >> k != 0, phys >> k
                if not odd and two:
                    hi >>= 1
                    r = hi >> (WINDOW - 1) & 1
                    hi ^= p << (WINDOW - 1)
                    p &= 1 - r
                if (z and (hi or p)) or (p and hi & mask != mask):
                    return False
                if odd:
                    word = u >> k if z else hi
                    borrow = en if word == 0 else (word if z else word - (p << WINDOW)) - (Vn >> (k - 1))
                    if borrow not in (0, 1):
                        return False
                    if borrow and word & mask == 0 and (hi >> WINDOW or u >> cut):
                        p = 1
                phys = Vn + (p << (cut - 1))
                if leave:
                    r = phys >> (cut - 1) & 1
                    if phys >> (w - 1) != p & (1 - r) or (phys ^ (p << (cut - 1)) ^ (p & (1 - r)) << (w - 1)) != Vn:
                        return False
            u, v, e = un, vn, en
    return True


def lag_blocks(plan):
    """[(first tick, c, front)] or None when the far fields are not whole blocks."""
    blocks, t = [], 0
    while t < len(plan):
        far = plan[t] and plan[t][10]
        if not far:
            t += 1
            continue
        front, c, _ = far
        for e in range(c):
            key = plan[t + e] if t + e < len(plan) else None
            if not (key and key[10] == (front - e, c, e == c - 1)):
                return None
        blocks.append((t, c, front))
        t += c
    return blocks


def lag_fields_ok(q, plan):
    blocks = lag_blocks(plan)
    if blocks is None:
        return False
    lazy = [t for t, key in enumerate(plan) if key and key[9]]
    for t, c, front in blocks:
        ticks = range(t, t + c)
        if not all(lazy[0] < s < lazy[-1] and plan[s][9] and plan[s][5] == joint_threshold(q, s) for s in ticks):
            return False
        K = [joint_threshold(q, s) for s in ticks]
        if any(front - 2 * e < K[e] + WINDOW + 2 - (e == c - 1) for e in range(c)):
            return False
        if any(front - (c - 1) - e < K[e - 1] for e in range(1, c)):
            return False
        if plan[t + c] is None or plan[t + c][0] + c < front - c or front - 2 * c < 0:
            return False
        if any(front - e > plan[t + e][0] + 1 + c for e in range(c)):
            return False
    return True


def lag_states_ok(q, plan, xs):
    blocks = lag_blocks(plan)
    if not blocks:
        return blocks is not None
    start = {t: (c, front) for t, c, front in blocks}
    for x in xs:
        u, v, e = q, x, 0
        left = 0
        for t, key in enumerate(plan):
            if t in start:
                left, front = start[t]
                k = 0
            odd, two = v & 1, int(v % 4 == 0)
            if left:
                kappa = key[5]
                z = u >> kappa != 0
                if front - k < kappa + WINDOW + 2 - (left == 1) or (z and k and u >> (front - k)):
                    return False
                if z and (v - e) >> kappa:
                    return False
                k += (1 - odd) & two & (1 - z)
                front -= 1
                left -= 1
            if odd:
                u, v, e = min(u, v), abs(u - v) // 2, int(u > v)
            else:
                v //= 4 if two else 2
    return True


def secp_inputs():
    q, rng = SECP, random.Random(2026)
    xs = {x for k in range(256) for d in (-1, 0, 1) for x in ((1 << k) + d, q - (1 << k) + d, (q >> k) + d)}
    a, b = 1, 1
    while b < 4 * q:
        xs |= {q * a // b + d for d in (-1, 0, 1)}
        a, b = b, a + b
    xs |= {rng.randrange(1, q) for _ in range(2000)}
    return sorted(x for x in xs if 0 < x < q)


def schedule_ok(q, xs, lag=None):
    plan = plan_of(q, lag)
    return (thresholds_ok(q, plan) and states_ok(q, plan, xs)
            and lazy_fields_ok(q, plan) and lazy_states_ok(q, plan, xs)
            and lag_fields_ok(q, plan) and lag_states_ok(q, plan, xs))


def main():
    rows = [{"q": q, "high_ticks": sum(1 for key in plan_of(q) if key and key[6]),
             "lazy_ticks": sum(1 for key in plan_of(q) if key and key[9]),
             "schedule_ok": schedule_ok(q, range(1, q))} for q in TOY]
    for lag in (1, 2, 3, 7):
        for q in TOY + [15361, 16381]:
            ticks = sum(1 for key in plan_of(q, lag) if key and key[10])
            if ticks:
                rows.append({"q": q, "lag_block": lag, "lag_ticks": ticks, "schedule_ok": schedule_ok(q, range(1, q), lag)})
    xs = secp_inputs()
    for lag in (7, 15):
        rows.append({"q": "secp256k1", "lag_block": lag, "lag_ticks": sum(1 for key in plan_of(SECP, lag) if key and key[10]),
                     "inputs": len(xs), "schedule_ok": schedule_ok(SECP, xs, lag)})
    rows.append({"q": "secp256k1", "high_ticks": sum(1 for key in plan_of(SECP) if key and key[6]),
                 "lazy_ticks": sum(1 for key in plan_of(SECP) if key and key[9]),
                 "lag_ticks": sum(1 for key in plan_of(SECP) if key and key[10]), "inputs": len(xs), "schedule_ok": schedule_ok(SECP, xs)})
    for row in rows:
        print(json.dumps(row))
    return 0 if all(row["schedule_ok"] for row in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
