"""Bit-sliced basis-state simulator: wire w is a Python int whose bit j is the
value of w in test lane j. Shares no arithmetic with the generator."""
from .builder import CCX, CX, X


def run(ops, num_wires, init, lanes):
    s = [0] * num_wires
    for w, v in init.items():
        s[w] = v
    lanes_all = (1 << lanes) - 1
    it = iter(ops)
    for op, a, b, c in zip(it, it, it, it):
        if op == CCX:
            s[c] ^= s[a] & s[b]
        elif op == CX:
            s[b] ^= s[a]
        elif op == X:
            s[a] ^= lanes_all
        else:
            s[a], s[b] = s[b], s[a]
    return s


def pack(values, wires):
    """Transpose per-lane integers into per-wire lane masks."""
    out = {w: 0 for w in wires}
    for lane, v in enumerate(values):
        for i, w in enumerate(wires):
            if v >> i & 1:
                out[w] |= 1 << lane
    return out


def unpack(state, wires, lanes):
    vals = [0] * lanes
    for i, w in enumerate(wires):
        m = state[w]
        lane = 0
        while m:
            if m & 1:
                vals[lane] |= 1 << i
            m >>= 1
            lane += 1
    return vals
