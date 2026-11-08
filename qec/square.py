"""Y^2 mod q on a fold of 2^n = f, for q = 2^n - f with f*f <= q.

Write ~Y = 2^n - 1 - Y. Each pair i < j gives y_i*~y_j + y_j*~y_i = y_i ^ y_j, and the
diagonal gives 0, so Y*~Y = sum over i < j of (y_i ^ y_j)*2^(i + j), a sum of rows with no
control on any of them. Q = Y + Y*~Y = Y*(2^n - Y), so Y^2 = 2^n*(Y - hi) - lo for
Q = lo + 2^n*hi, and Y - hi = D is at least 0. The rest folds 2^n to f:

    Z = (2^n - lo) + f*D = lo1 + 2^n*E, E <= f
    V = lo1 + f*E, V = Y^2 + f mod q, 0 <= V - f, V < 2^n + f*f

If V >= 2^n, dropping that bit subtracts q and f together and leaves V - 2^n <= f*f - 1 < q.
Otherwise V - f < q, so subtracting f is enough. So no comparison with q is needed.
"""
from .arith import _carries, _uncarry, add, add_const, gt_into
from .two_bit import Prefix


def signed_digits(v):
    """Non-adjacent form of v, as (position, sign)."""
    out, s = [], 0
    while v:
        if v & 1:
            d = 2 - (v & 3)
            out.append((s, d))
            v -= d
        v >>= 1
        s += 1
    return out


class SquareFold:
    def __init__(self, mod):
        self.mod, self.n = mod, mod.n
        self.f = f = (1 << mod.n) - mod.q
        self.b, self.m = f.bit_length(), (f * f).bit_length()
        # The low digit is +1 whatever f is, so its carry in can hold Z's + 1.
        self.digits = [(0, 1)] + signed_digits(f - 1)
        self.junk = 2 * mod.n + self.b + self.m + 3

    @staticmethod
    def applies(q):
        f = (1 << q.bit_length()) - q
        return f * f <= q

    def modadd(self, C, a, b, ctrl):
        """Mod.modadd with _reduce's carries above f's window held by a Prefix, so it runs
        beside the wires compute leaves live. The same flag, the same write, the same CCX."""
        f, w = self.f, self.b
        top, t = C.alloc(2)
        add(C, a, b, top, ctrl)
        cs = _carries(C, b[:w + 1], f)
        carry = Prefix(C, b[w:])
        # s's top bit and the carry into it are never both 1, so their XOR is the flag.
        C.ccx(cs[w], carry.top, top)
        C.ccx(top, cs[w], t)
        carry.write(t)
        C.ccx(top, cs[w], t)
        for i in reversed(range(w)):
            _uncarry(C, b, f, cs, i + 1)
            if cs[i] is not None:
                C.ccx(top, cs[i], b[i])
            if f >> i & 1:
                C.cx(top, b[i])
        gt_into(C, a, b, top, ctrl)
        C.free([top, t])

    def _rows(self, C, Y, Q):
        """Q = Y*(2^n - Y) on clean wires. After row i, Q < 2^(n + i + 1), so each carry
        lands on a wire no earlier row reached. Y's wires walk through their Gray code."""
        n = self.n
        for y, q in zip(Y, Q):
            C.cx(y, q)
        for i in range(n - 1):
            for w in Y[i + 1:]:
                C.cx(Y[i], w)
            add(C, Y[i + 1:], Q[2 * i + 1:n + i], cout=Q[n + i])
        for i in range(1, n):
            C.cx(Y[i - 1], Y[i])

    def _multiple(self, C, a, T, a_max, t_max=0, cin=None):
        """T += f*a mod 2^len(T), plus 1 if cin holds 1, for a <= a_max and T <= t_max on entry.

        Around each digit T holds t + d*a, plus 1 once cin is in, with d the digits so far.
        Those values at the corners of the box bound every input. A subtract adds onto ~T, so
        it bounds -v - 1 instead. Where every bound lies in [-2^i, 2^i) for a bit i at or above
        a's top, the digit's bits from i up all hold the sign. A sign that is always 0 never
        moves, so the digit adds below bit i only. Otherwise bit i takes the carry out and the
        bits above copy it. A digit onto a T that is 0, with no carry in, copies a."""
        inc, done, d = int(cin is not None), 0, 0
        for s, sign in self.digits:
            box = [(t, x) for t in (0, t_max) for x in (0, a_max)]
            before = [t + done + x * d for t, x in box]
            d += sign << s
            after = [t + inc + x * d for t, x in box]
            done = inc
            target = T[s:]
            if not target:
                continue
            if not any(before) and cin is None and sign > 0:
                for x, w in zip(a, target):
                    C.cx(x, w)
                continue
            ends = [v >> s for v in before + after]
            lo, hi = (min(ends), max(ends)) if sign > 0 else (-max(ends) - 1, -min(ends) - 1)
            top = next((i for i in range(len(a), len(target)) if -(1 << i) <= lo and hi < 1 << i),
                       len(target))
            carry = lo < 0 and top < len(target)
            above = target[top + 1:] if carry else []
            flips = [] if sign > 0 else target if carry else target[:top]
            for w in flips:
                C.x(w)
            for w in above:
                C.cx(target[top], w)
            pad = C.alloc(max(0, top - len(a)))
            add(C, a[:top] + pad, target[:top], cout=target[top] if carry else None, cin=cin)
            C.free(pad)
            for w in above:
                C.cx(target[top], w)
            for w in flips:
                C.x(w)
            cin = None

    def compute(self, C, Y, J):
        """J clean on entry. On exit J[:n] holds Y^2 mod q, canonical, for any n-bit Y.
        The other wires of J hold values fixed by Y, for the reverse of this block."""
        n, f, b, m = self.n, self.f, self.b, self.m
        lo, hi, E, T = J[:n], J[n:2 * n], J[2 * n:2 * n + b], J[2 * n + b:2 * n + b + m]
        t, c, k = J[2 * n + b + m:]
        self._rows(C, Y, lo + hi)
        # hi = Y - hi = ~(hi + ~Y), and lo = ~lo.
        for w in Y:
            C.x(w)
        add(C, Y, hi)
        for w in Y:
            C.x(w)
        for w in hi + lo:
            C.x(w)
        (one,) = C.alloc(1)
        C.x(one)
        # hi is any n-bit D, and lo + E holds ~lo, below 2^n.
        self._multiple(C, hi, lo + E, (1 << n) - 1, (1 << n) - 1, one)
        C.x(one)
        C.free([one])
        self._multiple(C, E, T, f)
        high = lo[m:]
        add(C, T, lo[:m], cout=t)
        if high:
            carry = Prefix(C, high)
            C.ccx(t, carry.top, c)
            carry.write(t)
        else:
            C.cx(t, c)
        # With c = 0, lo >= f. Take f off the low bits, and the borrow off the rest.
        C.x(c)
        add_const(C, lo[:b], (1 << b) - f, c, k)
        C.cx(c, k)
        high = lo[b:]
        for w in high:
            C.x(w)
        Prefix(C, high).write(k)
        for w in high:
            C.x(w)
        C.x(c)
