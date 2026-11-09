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
from .arith import _carries, _uncarry, add, add_const, add_lookahead, fan_into, gt_into, prefix_xor
from .two_bit import Prefix

# Add the square's rows as lookahead trees, on ROW_TEMPS clean wires beside the fold's own.
ROW_TREE = True
ROW_TEMPS = 120
# Fan each row's Gray-code step out through at most this many copies on clean temps, and undo the walk
# by a prefix network; 0 runs both as CX chains.
ROW_FAN_COPIES = 63
# Run the fold's adds as add_lookahead trees, on FOLD_TEMPS clean wires and the fold's clean junk.
FOLD_TREE = True
FOLD_TEMPS = 86


def _add(C, a, b, cin=None, temps=()):
    """b += a mod 2^len(b), plus 1 if cin holds 1: add_lookahead's tree on temps when FOLD_TREE."""
    if FOLD_TREE:
        add_lookahead(C, a, b, None, list(temps), tree=True, cin=cin)
    else:
        add(C, a, b, cin=cin)


def _fan_copies(targets, cap):
    """The number of copies, 2^j - 1 <= cap, that fans src into targets in the fewest CX layers:
    j to make them and ceil(targets / 2^j) to fan out (unmaking them overlaps the add)."""
    best, j = (targets, 0), 1
    while (1 << j) - 1 <= cap:
        best = min(best, (j + -(-targets // (1 << j)), (1 << j) - 1))
        j += 1
    return best[1]


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

    def _rows(self, C, Y, Q, temps):
        """Q = Y*(2^n - Y) on clean wires. After row i, Q < 2^(n + i + 1), so each carry
        lands on a wire no earlier row reached. Y's wires walk through their Gray code.
        The rows add in series, so each adds with add_lookahead, on clean temps, as a tree
        with ROW_TEMPS more clean wires when ROW_TREE."""
        n = self.n
        extra = C.alloc(ROW_TEMPS) if ROW_TREE else []
        temps = list(temps) + extra
        for y, q in zip(Y, Q):
            C.cx(y, q)
        for i in range(n - 1):
            if ROW_FAN_COPIES:
                # The tree takes temps from the front, so the copies sit at the back.
                k = _fan_copies(n - 1 - i, min(ROW_FAN_COPIES, len(temps)))
                fan_into(C, Y[i], Y[i + 1:], temps[len(temps) - k:])
            else:
                for w in Y[i + 1:]:
                    C.cx(Y[i], w)
            add_lookahead(C, Y[i + 1:], Q[2 * i + 1:n + i], Q[n + i], temps, tree=ROW_TREE)
        if ROW_FAN_COPIES:
            prefix_xor(C, Y)
        else:
            for i in range(1, n):
                C.cx(Y[i - 1], Y[i])
        C.free(extra)

    def _multiple(self, C, a, T, cin=None, temps=()):
        """T += f*a mod 2^len(T), plus 1 if cin holds 1. temps: clean wires for the adds' trees."""
        for s, sign in self.digits:
            target = T[s:]
            if not target:
                continue
            pad = C.alloc(max(0, len(target) - len(a)))
            b = a[:len(target)] + pad
            if sign > 0:
                _add(C, b, target, cin, temps)
            else:
                for w in target:
                    C.x(w)
                _add(C, b, target, None, temps)
                for w in target:
                    C.x(w)
            C.free(pad)
            cin = None

    def compute(self, C, Y, J):
        """J clean on entry. On exit J[:n] holds Y^2 mod q, canonical, for any n-bit Y.
        The other wires of J hold values fixed by Y, for the reverse of this block."""
        n, f, b, m = self.n, self.f, self.b, self.m
        lo, hi, E, T = J[:n], J[n:2 * n], J[2 * n:2 * n + b], J[2 * n + b:2 * n + b + m]
        t, c, k = J[2 * n + b + m:]
        # E, T, t, c and k are clean until the rows are done
        self._rows(C, Y, lo + hi, J[2 * n:])
        fold = C.alloc(FOLD_TEMPS) if FOLD_TREE else []
        # hi = Y - hi = ~(hi + ~Y), and lo = ~lo.
        for w in Y:
            C.x(w)
        _add(C, Y, hi, None, J[2 * n:] + fold)
        for w in Y:
            C.x(w)
        for w in hi + lo:
            C.x(w)
        (one,) = C.alloc(1)
        C.x(one)
        self._multiple(C, hi, lo + E, one, T + [t, c, k] + fold)
        C.x(one)
        C.free([one])
        self._multiple(C, E, T, None, [t, c, k] + fold)
        high = lo[m:]
        if FOLD_TREE:
            add_lookahead(C, T, lo[:m], t, [c, k] + fold, tree=True)
        else:
            add(C, T, lo[:m], cout=t)
        C.free(fold)
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
