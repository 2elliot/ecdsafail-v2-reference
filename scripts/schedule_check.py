"""The high-word schedule against the joint bound, with no use of qec.model.

A flagged tick's kappa must be min(w, K_t), K_t the least K with 2^(3K + 2t) > 2 q^2 (q - 1), and a
high tick's must be K_t itself. On the semantic Euclid recurrence, with the held offset tracked
beside it, no state before a flagged tick has both registers at 2^kappa or more, and a subtract
there leaves U below 2^kappa.

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
TOY = [3, 5, 7, 11, 13, 31, 43, 61, 127, 151, 251, 509, 1021, 4093, 8117]


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
    return thresholds_ok(q, plan) and states_ok(q, plan, xs)


def main():
    rows = [{"q": q, "high_ticks": sum(1 for key in plan_of(q) if key and key[6]),
             "schedule_ok": schedule_ok(q, range(1, q))} for q in TOY]
    xs = secp_inputs()
    rows.append({"q": "secp256k1", "high_ticks": sum(1 for key in plan_of(SECP) if key and key[6]),
                 "inputs": len(xs), "schedule_ok": schedule_ok(SECP, xs)})
    for row in rows:
        print(json.dumps(row))
    return 0 if all(row["schedule_ok"] for row in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
