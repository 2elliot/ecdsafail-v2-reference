"""Classical integer model of the exact dialog in-place multiplier.

M(x, y) = (x, x*y mod q) for x != 0 and M(0, y) = (0, y). The reference Euclid dialog is the oracle for the iteration bound and Bezout
identity. The held-offset recurrence and sweep describe the alternative
physical representation used by the shared-chain circuit.
"""


def iteration_bound(q):
    """Proven cap on active Euclid iterations: u*v at least halves per step."""
    return (q * (q - 1)).bit_length()


def euclid(q, x, iters):
    u, v, dialog = q, x, []
    for _ in range(iters):
        b0 = v & 1
        sw = b0 & (u > v)
        dialog.append((b0, sw))
        if sw:
            u, v = v, u
        if b0:
            v -= u
        v >>= 1
    return u, v, dialog


def offset_pass(q: int, x: int, iters: int) -> tuple[list[tuple[int, int, int]], list[tuple[int, int]]]:
    """States before each step and at the end, and the unswept tape cells."""
    U, V, e = q, x, 0
    states = [(U, V, e)]
    cells = []
    for _ in range(iters):
        V ^= e
        odd = V & 1
        if odd:
            if U > V:
                U, V, e = V, U - V - 2, 1
            else:
                V, e = V - U, 0
        V //= 2
        cells.append((odd, e))
        states.append((U, V, e))
    return states, cells


def sweep_tape(cells: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return the reference (parity, swap) tape without changing the input."""
    out = list(cells)
    for i in range(len(out) - 1, 0, -1):
        odd, held = out[i]
        out[i] = (odd, held ^ ((1 - odd) * out[i - 1][1]))
    return out


def bezout(q, y, dialog):
    r, s = y, 0
    for b0, sw in reversed(dialog):
        s = 2 * s % q
        if b0:
            s = (s + r) % q
        if sw:
            r, s = s, r
    return r, s


def inplace_mul(q, x, y, iters=None):
    iters = iters or iteration_bound(q)
    u, v, dialog = euclid(q, x or 1, iters)
    assert u == 1 and v == 0, "euclid did not reach (1, 0) within the bound"
    r, s = bezout(q, y, dialog)
    assert r == 0
    return s


def width_schedule(q, iters):
    """Bits needed before iteration i while v != 0.

    For nonzero canonical x the terminal state (1, 0) also fits every
    emitted step. The multiplier remaps x = 0 to this nonzero domain.
    """
    n = q.bit_length()
    return [min(n, (q * (q - 1) >> i).bit_length()) for i in range(iters)]
