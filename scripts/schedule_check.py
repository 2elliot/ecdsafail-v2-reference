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


def plan_of(q):
    return DialogMultiplier(Circuit(), Mod(q))._tick_plan()


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


def secp_inputs():
    q, rng = SECP, random.Random(2026)
    xs = {x for k in range(256) for d in (-1, 0, 1) for x in ((1 << k) + d, q - (1 << k) + d, (q >> k) + d)}
    a, b = 1, 1
    while b < 4 * q:
        xs |= {q * a // b + d for d in (-1, 0, 1)}
        a, b = b, a + b
    xs |= {rng.randrange(1, q) for _ in range(2000)}
    return sorted(x for x in xs if 0 < x < q)


def schedule_ok(q, xs):
    plan = plan_of(q)
    return (thresholds_ok(q, plan) and states_ok(q, plan, xs)
            and lazy_fields_ok(q, plan) and lazy_states_ok(q, plan, xs))


def main():
    rows = [{"q": q, "high_ticks": sum(1 for key in plan_of(q) if key and key[6]),
             "lazy_ticks": sum(1 for key in plan_of(q) if key and key[9]),
             "schedule_ok": schedule_ok(q, range(1, q))} for q in TOY]
    xs = secp_inputs()
    rows.append({"q": "secp256k1", "high_ticks": sum(1 for key in plan_of(SECP) if key and key[6]),
                 "lazy_ticks": sum(1 for key in plan_of(SECP) if key and key[9]), "inputs": len(xs), "schedule_ok": schedule_ok(SECP, xs)})
    for row in rows:
        print(json.dumps(row))
    return 0 if all(row["schedule_ok"] for row in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
