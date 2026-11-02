"""Classical integer model of the exact dialog in-place multiplier.

M(x, y) = (x, x*y mod q) for x != 0 and M(0, y) = (0, y). The circuit in
point_add.py mirrors these loops gate for gate, so this model is the oracle
for the iteration bound and the Bezout identity.
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
