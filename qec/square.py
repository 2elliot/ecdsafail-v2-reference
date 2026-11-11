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
from .arith import _carries, _carry_into, _uncarry, add, add_const, gt_into, load_const
from .two_bit import Prefix


def _keep_add(C, a, b, cout, K):
    """b += a, cout ^= carry out if given, each inner carry left on a junk wire popped from K."""
    c = None
    for k, (x, y) in enumerate(zip(a, b)):
        if k == len(a) - 1 and cout is None:
            C.cx(x, y)
            if c is not None:
                C.cx(c, y)
            break
        t = cout if k == len(a) - 1 else K.pop(0)
        if c is None:
            C.ccx(x, y, t)
            C.cx(x, y)
        else:
            # MAJ(x, y, c) = c ^ (x ^ c)(y ^ c); the sum x ^ y ^ c lands on y.
            C.cx(c, x)
            C.cx(c, y)
            C.ccx(x, y, t)
            C.cx(c, t)
            C.cx(x, y)
            C.cx(c, y)
            C.cx(c, x)
        c = t


def _rows(C, Y, Q, copy=True, K=None):
    """Q = Y*(2^w - Y) mod 2^len(Q) on len(Q) <= 2w clean wires, w = len(Y), or Y*(2^w - 1 - Y)
    without the copy. After row i, Q < 2^(w + i + 1), so each carry lands on a wire no earlier row
    reached; the uncopied sums are smaller. Y's wires walk through their Gray code. K, if given, is
    a list of junk wires, and a row whose inner carries fit there keeps them, at one CCX a bit."""
    n, L = len(Y), len(Q)
    if copy:
        for y, q in zip(Y, Q):
            C.cx(y, q)
    for i in range(n - 1):
        for w in Y[i + 1:]:
            C.cx(Y[i], w)
        width = min(n - 1 - i, L - 2 * i - 1)
        if width <= 0:
            continue
        a, b = Y[i + 1:i + 1 + width], Q[2 * i + 1:2 * i + 1 + width]
        cout = Q[n + i] if n + i < L else None
        if K is not None and len(K) >= width - 1:
            _keep_add(C, a, b, cout, K)
        else:
            add(C, a, b, cout=cout)
    for i in range(1, n):
        C.cx(Y[i - 1], Y[i])


def _change(C, a, b, sign=1, cin=None, signed=False):
    """b += sign*a (+ cin) mod 2^len(b); a zero-extended, or sign-extended when signed."""
    assert len(a) <= len(b)
    pad = C.alloc(len(b) - len(a))
    if signed:
        for p in pad:
            C.cx(a[-1], p)
    if sign < 0:
        for w in b:
            C.x(w)
    add(C, list(a) + pad, b, cin=cin)
    if sign < 0:
        for w in b:
            C.x(w)
    if signed:
        for p in pad:
            C.cx(a[-1], p)
    C.free(pad)


def _add_const_cin(C, reg, value, cin, cout=None, keep=()):
    """reg += value + cin mod 2^len(reg), cout ^= the carry out if given, cin a wire left as it was.
    The carry chain costs what add_const's controlled increment does: one CCX a bit each way. The
    carries into bits 1 to len(keep) land on keep's clean wires and stay there as junk."""
    cs = [cin]
    for i in range(len(reg) - 1):
        (t,) = keep[i:i + 1] or C.alloc(1)
        _carry_into(C, reg[i], value >> i & 1, cs[i], t)
        cs.append(t)
    if cout is not None:
        _carry_into(C, reg[-1], value >> (len(reg) - 1) & 1, cs[-1], cout)
    for i in reversed(range(len(reg))):
        if value >> i & 1:
            C.x(reg[i])
        C.cx(cs[i], reg[i])
        if i > len(keep):
            _carry_into(C, reg[i - 1], value >> (i - 1) & 1, cs[i - 1], cs[i])
            C.free([cs[i]])


class _Lend:
    """Lends clean wires to the allocator for a stretch of the build and takes the same wires
    back. Every allocation in the stretch is returned clean before it ends, so the lent wires
    are clean and free again, and the scratch the stretch needs comes out of them first."""

    def __init__(self, C, wires):
        self.C, self.wires = C, list(wires)

    def __enter__(self):
        self.C.free(self.wires)

    def __exit__(self, *exc):
        stack = self.C.free_stack
        for w in self.wires:
            stack.remove(w)


def _value_bits(w, kind, o=0):
    """Bits of the largest kind(Y) + o over w-bit Y: R(Y) <= 2^(2w - 2), T(Y) < 2^(2w - 2). A
    negative o wraps mod 2^(2w), so then all 2w."""
    if o < 0:
        return 2 * w
    top = 1 << (2 * w - 2) if kind == "R" else ((1 << w) - 1) ** 2 // 4
    return (top + o).bit_length()


def _widths(w, thr, kind, exact, L):
    """(V's width, W's width) at a node of width w over len(Q) = L. V = A*~B + B*~A (+ B for R)
    lies in [0, 2^w), so an exact R node's V fits w bits; T nodes and inexact R nodes carry a
    constant on V and keep w + 1. Q is only needed mod 2^L, so V is only needed mod 2^(L - h).
    The S child's Q must hold its B child's T exactly, which takes 2ws - 2 bits."""
    h = (w + 1) // 2
    ws = h + 1
    lv = min(w + (kind == "T" or not exact), L - h)
    return lv, min(2 * ws, max(lv, 2 * ws - 2 if ws > thr else 0))


def karatsuba_junk(w, thr, kind="R", exact=True, L=None):
    """Junk wires _karatsuba takes for width w with leaves of width at most thr, over L wires of Q."""
    if w <= thr:
        return 0
    L = 2 * w if L is None else L
    h = (w + 1) // 2
    inner = exact or kind == "T"
    lv, lw = _widths(w, thr, kind, inner, L)
    return (karatsuba_junk(h, thr, kind, inner, min(2 * h, L)) + karatsuba_junk(w - h, thr, "T", True, L - 2 * h)
            + h + 1 + lw + karatsuba_junk(h + 1, thr, "R" if kind == "T" else "T", inner, lw))


def karatsuba_offset(w, thr, kind="R", exact=True, cins=frozenset(), path=""):
    """The constant o that _karatsuba adds to its Q, or None where an inexact plan breaks a range.

    An exact node leaves o = 0. An inexact T node skips its 2^ws - 1 fixup, so its V carries
    delta = 2^ws - 1 more. An inexact R node takes its A and S children inexact and its B child
    exact, so delta = o(S) - o(A), and every o is o(A) + 2^(2h) o(B) + 2^h delta. A child's value
    plus its o must stay below 2^(2w), which holds while o < 3*2^(2w - 2), as R and T are at most
    2^(2w - 2), and o >= 0; V's true value lies in [0, 2^w), so V + delta fits w + 1 bits while |delta| <= 2^w,
    signed when delta < 0.

    An R node whose path (A for an A child, S for an S child, from the top's "") is in cins takes a
    carry in of 1 on Q[h:] += V, which adds 2^h to its o. T nodes take none and need exact children."""
    if w <= thr:
        return 0
    h = (w + 1) // 2
    m, ws = w - h, h + 1
    inner = exact or kind == "T"
    oa = karatsuba_offset(h, thr, kind, inner, cins, path + "A")
    os = karatsuba_offset(ws, thr, "R" if kind == "T" else "T", inner, cins, path + "S")
    if oa is None or os is None:
        return None
    if not (0 <= oa < 3 << (2 * h - 2) and 0 <= os < 3 << (2 * ws - 2)):
        return None
    if kind == "T":
        if oa or os or path in cins:
            return None
        return 0 if exact else ((1 << ws) - 1) << h
    delta = os - oa
    if abs(delta) > 1 << w:
        return None
    return oa + ((delta + (path in cins)) << h)


def _karatsuba(C, Y, Q, J, thr, kind="R", K=None, exact=True, cins=frozenset(), path=""):
    """Q = R(Y) = Y*(2^w - Y) (kind R) or T(Y) = Y*(2^w - 1 - Y) (kind T), plus
    karatsuba_offset(w, thr, kind, exact, cins, path), mod 2^len(Q), on len(Q) <= 2w clean wires.

    Y = A + 2^h*B with h = ceil(w/2), m = w - h, and S = A + B on ws = h + 1 bits. Then
        R(Y) = R(A) + 2^(2h)*T(B) + 2^h*V,   V = T(S) - R(A) - T(B) - [h != m]*2^m*S
        T(Y) = T(A) + 2^(2h)*T(B) + 2^h*V,   V = R(~S) - T(A) - T(B) - [h != m]*2^m*S
                                                  - (2^ws - 1)
    where ~S = 2^ws - 1 - S, so R(~S) = T(S) - S + 2^ws - 1. V is formed mod 2^len(V) (_widths).
    R(A) and T(B) go straight onto Q's halves; S and the S child's output take wires from J, a
    list of clean junk wires, and keep values fixed by Y. Leaves of width at most thr are the XOR
    rows. Bits of Q above the value's top (_value_bits) stay 0 and are not touched."""
    w, L = len(Y), len(Q)
    top = min(L, _value_bits(w, kind, karatsuba_offset(w, thr, kind, exact, cins, path)))
    if w <= thr:
        _rows(C, Y, Q[:top], kind == "R", K)
        return
    h = (w + 1) // 2
    m = w - h
    inner = exact or kind == "T"
    lv, lw = _widths(w, thr, kind, inner, L)
    # Q mod 2^top needs T(B) mod 2^(top - h), and Q[2h:] holds T(B) mod 2^(L - 2h).
    assert _value_bits(m, "T") <= L - 2 * h or top - h <= L - 2 * h
    A, B = Y[:h], Y[h:]
    _karatsuba(C, A, Q[:min(2 * h, L)], J, thr, kind, K, inner, cins, path + "A")
    _karatsuba(C, B, Q[2 * h:], J, thr, "T", K)
    ws = h + 1
    S = [J.pop(0) for _ in range(ws)]
    for x, s in zip(A, S):
        C.cx(x, s)
    pad = C.alloc(h - m)
    add(C, B + pad, S[:h], cout=S[h])
    C.free(pad)
    W = [J.pop(0) for _ in range(lw)]
    if kind == "T":
        for x in S:
            C.x(x)
    _karatsuba(C, S, W, J, thr, "R" if kind == "T" else "T", K, inner, cins, path + "S")
    if kind == "T":
        for x in S:
            C.x(x)
    V = W[:lv]
    _change(C, Q[:min(2 * h, lv)], V, -1)
    _change(C, Q[2 * h:L][:lv], V, -1)
    if h != m and lv > m:
        _change(C, S[:lv - m], V[m:], -1)
    out = V[:top - h]
    if kind == "R":
        delta = (karatsuba_offset(ws, thr, "T", inner, cins, path + "S")
                 - karatsuba_offset(h, thr, "R", inner, cins, path + "A"))
        one = C.alloc(1) if path in cins else [None]
        if one[0] is not None:
            C.x(one[0])
        _change(C, out, Q[h:top], cin=one[0], signed=delta < 0)
        if one[0] is not None:
            C.x(one[0])
            C.free(one)
        return
    if not exact:
        _change(C, out, Q[h:top])
        return
    # V holds the cross term plus 2^ws - 1. Take 2^ws off its top bits, which leaves at
    # least -1 on w + 1 bits, then add it sign-extended with a carry in of 1.
    if lv > ws:
        add_const(C, V[ws:], (1 << (lv - ws)) - 1)
    (one,) = C.alloc(1)
    C.x(one)
    _change(C, out, Q[h:top], cin=one, signed=True)
    C.x(one)
    C.free([one])


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
    # Kept row carries: each junk wire saves one CCX each way. 886 fills the 4,096-wire cap at
    # n = 256, where the square's own peak is 3,229 wires and M's peak lies outside it.
    KEEP = 886
    # The last fold ladder's low carries stay on junk wires through the modadd, one CCX saved each
    # way per wire. The modadd's 41 scratch wires and these 181 fit in the 222 that ladder itself
    # peaks at, n = 256, so they cost no row carries.
    LADDER_KEEP = 181
    def __init__(self, mod):
        self.mod, self.n = mod, mod.n
        self.f = f = (1 << mod.n) - mod.q
        self.b, self.m = f.bit_length(), (f * f).bit_length()
        # The low digit is +1 whatever f is, so its carry in can hold Z's + 1.
        self.digits = [(0, 1)] + signed_digits(f - 1)
        # Karatsuba leaves of at most n // 6 bits: 42 at n = 256, and 3 at the toy sizes so
        # the recursion runs there. KEEP junk wires hold row carries, one CCX saved per wire.
        self.thr = max(3, mod.n // 6)
        self.lkeep = min(self.LADDER_KEEP, max(0, mod.n - self.b - 1))
        self.fixed_junk = 2 * mod.n + self.b + self.m + 3 + self.lkeep
        self.junk = self.fixed_junk + self.karatsuba_junk(self.thr) + self.KEEP

    def karatsuba_junk(self, thr):
        exact = karatsuba_offset(self.n, thr, "R", False) != self._offset(thr)
        return karatsuba_junk(self.n, thr, "R", exact)

    def _offset(self, thr):
        o = karatsuba_offset(self.n, thr, "R", False)
        return o if o is not None and 0 <= o < 1 << self.n else 0

    @property
    def offset(self):
        """The o with compute's J[:n] = (Y^2 - o) mod q: the inexact Karatsuba's constant, when it
        lies below 2^n, so that hi = floor(Q / 2^n) stays at most Y; else 0 from an exact one."""
        o = karatsuba_offset(self.n, self.thr, "R", False)
        return o if o is not None and 0 <= o < 1 << self.n else 0

    def absorb(self, c):
        """(cins, kappa) with compute(C, Y, J, kappa, cins) leaving J[:n] = (Y^2 - c) mod q, or None.

        compute adds kappa to V before the canonical step, so J = (Y^2 - o + kappa) mod q for the
        Karatsuba's o, and kappa = (o - c) mod q. That holds for kappa <= q - f*f, as V + kappa
        then stays below 2^n + q. Below bit m, kappa = f*kE + klo: kE is added onto E, and klo
        loaded onto the clean T, and
        T stays below 2^m while f*f + (kappa mod 2^m) < 2^m and E + kE fits b bits; the bits from
        m up ride the carry ladder. The A child's carry in moves kappa mod 2^m by about 2^m / 2,
        for the half of the cs that miss that window with none."""
        n, f, q = self.n, self.f, self.mod.q
        o_nat = karatsuba_offset(n, self.thr, "R", False)
        exact = self.offset != o_nat
        for cins in (frozenset(), frozenset({"A"})):
            o = karatsuba_offset(n, self.thr, "R", exact, cins)
            if o is None or not 0 <= o < 1 << n:
                continue
            kappa = (o - c) % q
            if self.split(kappa) is not None:
                return cins, kappa
        return None

    def split(self, kappa):
        """(kE, klo, khi) for compute's kappa, or None outside the window absorb states."""
        f, m, n = self.f, self.m, self.n
        x = kappa % (1 << m)
        kE, klo, khi = x // f, x % f, kappa >> m
        if kappa > self.mod.q - f * f or f * f + x >= 1 << m or f + kE >= 1 << self.b:
            return None
        if khi and n <= m:
            return None
        return kE, klo, khi

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

    def compute(self, C, Y, J, kappa=None, cins=frozenset()):
        """J clean on entry. On exit J[:n] holds (Y^2 - offset) mod q, canonical, for any n-bit Y, or
        with kappa and cins from absorb(c), (Y^2 - c) mod q. The other wires of J hold values fixed by
        Y, for the reverse of this block."""
        n, f, b, m = self.n, self.f, self.b, self.m
        lo, hi, E, T = J[:n], J[n:2 * n], J[2 * n:2 * n + b], J[2 * n + b:2 * n + b + m]
        t, c, k = J[2 * n + b + m:2 * n + b + m + 3]
        ladder = list(J[2 * n + b + m + 3:self.fixed_junk])
        junk = list(J[self.fixed_junk:self.junk - self.KEEP])
        keep = list(J[self.junk - self.KEEP:]) if self.KEEP else None
        o = karatsuba_offset(self.n, self.thr, "R", False)
        kE, klo, khi = (0, 0, 0) if kappa is None else self.split(kappa)
        # The ladder's kept-carry wires are clean until the last ladder, so everything before it
        # takes its scratch from them first.
        with _Lend(C, ladder):
            self._fold(C, Y, J, junk, keep, o, cins, kE, klo, khi)
        high = lo[b:]
        for w in high:
            C.x(w)
        _add_const_cin(C, high, 0, k, keep=ladder)
        for w in high:
            C.x(w)
        C.x(c)

    def _fold(self, C, Y, J, junk, keep, o, cins, kE, klo, khi):
        n, f, b, m = self.n, self.f, self.b, self.m
        lo, hi, E, T = J[:n], J[n:2 * n], J[2 * n:2 * n + b], J[2 * n + b:2 * n + b + m]
        t, c, k = J[2 * n + b + m:2 * n + b + m + 3]
        _karatsuba(C, Y, lo + hi, junk, self.thr, "R", keep, self.offset != o, cins)
        assert not junk
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
        if kE:
            add_const(C, E, kE)
        load_const(C, T, klo)
        self._multiple(C, E, T, f + kE, klo)
        high = lo[m:]
        add(C, T, lo[:m], cout=t)
        if high:
            _add_const_cin(C, high, khi, t, c)
        else:
            C.cx(t, c)
        # With c = 0, lo >= f. Take f off the low bits, and the borrow off the rest.
        C.x(c)
        add_const(C, lo[:b], (1 << b) - f, c, k)
        C.cx(c, k)
