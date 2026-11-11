"""The rank-one start of a Bezout pass: its first J cells on unreduced integers.

The pass starts at (R, S) = (y, 0), so after its first j cells (R, S) = (U*y, V*y) for the integer
cofactor pair (U, V) of those cells, which is the Euclid pair before Euclid's last j ticks, so
U, V < 2^w with w = tick_widths(q, T)[T - j]. Unreduced, a tick is its add-or-shift and its swap:
no Prefix, low chains, compare or k machinery. One conversion after cell J reduces both registers:

  fold      X = hi*2^n + lo -> (X mod q, m = floor(X/q)), for q = 2^n - f with (2^w + 2)f <= 2^n:
            m = hi + e, e = [lo + (hi + 1)f >= 2^n], X mod q = lo + m*f mod 2^n
  cofactors (U, V) from the cells, on w-bit registers
  D         D = (U*s - V*r)/q, an integer in (-V, U), is (U*s - V*r)*q^-1 mod 2^(w + 1)
  mS        V*mR - U*mS = D, so U*mS - V*mR + D = 0 and mS -> U*mS mod 2^w in place ends at 0
  mR        the cells' tick mod U maps (mR, 0) to (0, V*mR mod U) = (0, D mod U) in place
  undo      D mod U, D and the cofactors

The middle cells use the tree's add-or-shift and held offset on wider words. The last cell is
semantic (a B1 doubles, adding nothing), so the reduced tick after it starts with no held lane.
docs/rank-one.md states the premises and the checks."""
from .arith import _maj, _uma, cswap
from .two_bit import Prefix, TwoBitMod


def plan(q, T, widths, J_max):
    """The largest J <= J_max the fold allows, or 0."""
    n = q.bit_length()
    f = (1 << n) - q
    for J in range(min(J_max, T - 2), 2, -1):
        w = widths[T - J]
        if (2 ** w + 2) * f <= 2 ** n and ((2 ** w) * f).bit_length() <= n:
            return J
    return 0


def _uma_ctl(C, ctl, c, b, a):
    C.ccx(c, b, a)
    C.cx(a, c)
    C.cx(a, b)
    C.cx(c, a)
    C.ccx(ctl, a, b)
    C.cx(c, a)


def add(C, a, b, ctl=None):
    """b += a (or ctl*a) mod 2^len(b), a zero-extended. 2 CCX a bit, 3 under a control."""
    n = len(b)
    a = list(a)
    ext = C.alloc(n - len(a)) if len(a) < n else []
    a = a + ext
    (c0,) = C.alloc(1)
    carry = [c0] + a[:-1]
    for i in range(n - 1):
        _maj(C, carry[i], b[i], a[i])
    i = n - 1
    C.cx(carry[i], a[i])
    if ctl is None:
        C.cx(a[i], b[i])
    else:
        C.ccx(ctl, a[i], b[i])
    C.cx(carry[i], a[i])
    for i in reversed(range(n - 1)):
        if ctl is None:
            _uma(C, carry[i], b[i], a[i])
        else:
            _uma_ctl(C, ctl, carry[i], b[i], a[i])
    C.free([c0] + ext)


def sub(C, a, b, ctl=None):
    for w in b:
        C.x(w)
    add(C, a, b, ctl)
    for w in b:
        C.x(w)


def carry_into(C, a, b, flag, ctl=None):
    """flag ^= [a + b >= 2^len(b)] (times ctl), a and b as found."""
    n = len(b)
    a = list(a)
    ext = C.alloc(n - len(a)) if len(a) < n else []
    a = a + ext
    (c0,) = C.alloc(1)
    carry = [c0] + a[:-1]
    for i in range(n):
        _maj(C, carry[i], b[i], a[i])
    if ctl is None:
        C.cx(a[n - 1], flag)
    else:
        C.ccx(ctl, a[n - 1], flag)
    for i in reversed(range(n)):
        C.ccx(carry[i], b[i], a[i])
        C.cx(a[i], carry[i])
        C.cx(a[i], b[i])
    C.free([c0] + ext)


def greater(C, a, b, flag, ctl=None):
    """flag ^= [a > b] for any widths, the carry of a + ~b."""
    k = max(len(a), len(b))
    pa, pb = C.alloc(k - len(a)), C.alloc(k - len(b))
    a, b = list(a) + pa, list(b) + pb
    for w in b:
        C.x(w)
    carry_into(C, a, b, flag, ctl)
    for w in b:
        C.x(w)
    C.free(pa + pb)


def inc(C, x, ctl):
    (one,) = C.alloc(1)
    C.cx(ctl, one)
    add(C, [one], x)
    C.cx(ctl, one)
    C.free([one])


def dec(C, x, ctl):
    (one,) = C.alloc(1)
    C.cx(ctl, one)
    sub(C, [one], x)
    C.cx(ctl, one)
    C.free([one])


def const_reg(C, k, value):
    r = C.alloc(k)
    for i in range(k):
        if value >> i & 1:
            C.x(r[i])
    return r


def free_const(C, r, value):
    for i in range(len(r)):
        if value >> i & 1:
            C.x(r[i])
    C.free(r)


def mul_acc(C, a, b, acc, sign=1):
    """acc += a*b (sign -1: acc -= a*b) mod 2^len(acc), row i on acc[i:]."""
    if sign < 0:
        for w in acc:
            C.x(w)
    for i, ai in enumerate(a[:len(acc)]):
        add(C, b[:len(acc) - i], acc[i:], ctl=ai)
    if sign < 0:
        for w in acc:
            C.x(w)


def mul_odd(C, u, x):
    """x = u*x mod 2^len(x) in place, u odd: x[i + 1:] += x_i*(u >> 1), i from the top down, so each
    step reads a bit that no later step changes."""
    k = len(x)
    for i in range(k - 2, -1, -1):
        part = u[1:][:k - i - 1]
        if part:
            add(C, part, x[i + 1:], ctl=x[i])


def mul_odd_const(C, x, c):
    """x = c*x mod 2^len(x) in place, c an odd constant."""
    k = len(x)
    half = c >> 1 & ((1 << max(k - 1, 0)) - 1)
    if not half:
        return
    reg = const_reg(C, k - 1, half)
    for i in range(k - 2, -1, -1):
        add(C, reg[:k - i - 1], x[i + 1:], ctl=x[i])
    free_const(C, reg, half)


def _csd(c):
    out, i = [], 0
    while c:
        if c & 1:
            d = 2 - (c & 3)
            out.append((i, d))
            c -= d
        c >>= 1
        i += 1
    return out


def _times(C, x, out, c, inverse):
    """out += c*x mod 2^len(out), out 0 before: a copy at one +1 digit of c, an add or subtract at
    each other."""
    W = len(out)
    terms = [t for t in _csd(c) if t[0] < W]
    pivot = next(t for t in terms if t[1] > 0)
    rest = [t for t in terms if t is not pivot]
    s0 = pivot[0]
    if not inverse:
        for i, w in enumerate(x[:W - s0]):
            C.cx(w, out[s0 + i])
        for s, d in rest:
            (add if d > 0 else sub)(C, x[:W - s], out[s:])
    else:
        for s, d in reversed(rest):
            (sub if d > 0 else add)(C, x[:W - s], out[s:])
        for i, w in enumerate(x[:W - s0]):
            C.cx(w, out[s0 + i])


def split(c, k):
    """(c_lo, t) when c = 2^t + c_lo with c_lo*x < 2^t for every k-bit x, else None. Then c*x is
    c_lo*x on its own wires and x itself at t."""
    t = c.bit_length() - 1
    lo_c = c - (1 << t)
    return (lo_c, t) if lo_c << k <= 1 << t else None


class RankOne:
    def __init__(self, q, J, w):
        self.q, self.J, self.w = q, J, w
        self.n = n = q.bit_length()
        self.f = (1 << n) - q
        self.W = ((2 ** w) * self.f).bit_length()
        # R and S through the unreduced cells: below 2^(w + 1)*q on a held lane, and two top wires
        # that the add-or-shift's word needs
        self.L = n + w + 3
        self.nar = TwoBitMod(q)

    # -- the unreduced cells -----------------------------------------------------------------------

    def tick(self, C, g, prev, p, R, S):
        """A middle cell on the held offset, as two_bit's tick without its reduction: the word takes
        R >> 1 off J (under a complement on a held lane) or doubles on J, then S = 2*word + g0*R0.
        Returns [this cell is B1] on the sw wire."""
        g0, b = g
        (sw,) = C.alloc(1)
        C.ccx(g0, b, sw)
        C.cx(sw, b)
        if p is not None:
            for x in S:
                C.cx(p, x)
        shim = _Width(len(R) - 2)
        TwoBitMod._add_or_shift(shim, C, b, R[:len(R) - 2], S, held_lane=p is not None)
        if p is not None:
            for x in S:
                C.cx(p, x)
            self.nar._drop_after_b1(C, prev, p)
        (l0,) = C.alloc(1)
        C.ccx(g0, R[0], l0)
        top = S[-1]
        S[:] = [l0] + S[:-1]
        C.free([top])
        C.cx(sw, b)
        cswap(C, sw, R, S)
        C.x(sw)
        C.cx(g0, sw)
        C.cx(b, sw)
        return sw

    def last_tick(self, C, g, prev, p, R, S):
        """The segment's last cell, semantic: 2S + R, 4S on J, and 2S on B1 with nothing added. On a
        held lane (g0 = 1) the word is H - R0 for H = S + R, and the complement makes 2H - R."""
        g0, b = g
        if p is not None:
            for x in S:
                C.cx(p, x)
        add(C, R[1:], S, ctl=g0)
        if p is not None:
            for x in S:
                C.cx(p, x)
            self.nar._drop_after_b1(C, prev, p)
        (j,) = C.alloc(1)
        C.x(g0)
        C.ccx(g0, b, j)
        C.x(g0)
        for i in range(len(S) - 1, 0, -1):
            C.cx(S[i - 1], S[i])
            C.ccx(j, S[i], S[i - 1])
            C.cx(S[i - 1], S[i])
        C.x(g0)
        C.ccx(g0, b, j)
        C.x(g0)
        C.free([j])
        (l0,) = C.alloc(1)
        C.ccx(g0, R[0], l0)
        top = S[-1]
        S[:] = [l0] + S[:-1]
        C.free([top])
        (sw,) = C.alloc(1)
        C.ccx(g0, b, sw)
        cswap(C, sw, R, S)
        C.ccx(g0, b, sw)
        C.free([sw])

    # -- the conversion ------------------------------------------------------------------------------

    def fold(self, C, X):
        """X (n + w wires) -> X mod q on X[:n] and floor(X/q) on X[n:]. The probe finds
        e = [lo + (hi + 1)f >= 2^n] against Prefix's AND of the high bits; the write adds m*f and its
        window carry clears e and moves the high bits."""
        n, f, W = self.n, self.f, self.W
        lo = X[:n]
        spare = C.alloc(1)
        hi = X[n:] + spare
        H = lo[W:]
        pre = Prefix(C, H) if H else None
        (one,) = C.alloc(1)
        C.x(one)
        inc(C, hi, one)
        C.x(one)
        C.free([one])
        cut = split(f, len(hi))
        if cut:
            lo_c, t = cut
            low = C.alloc((lo_c * ((1 << len(hi)) - 1)).bit_length())
            own = low + C.alloc(t - len(low))
            prod = own + hi
            assert len(prod) == W
        else:
            lo_c, low = f, C.alloc(W)
            own = prod = low

        def times(inverse=False):
            if low:
                _times(C, hi, low, lo_c, inverse)

        times()
        (e,) = C.alloc(1)
        (c1,) = C.alloc(1)
        carry_into(C, prod, lo[:W], c1)
        if H:
            C.ccx(pre.top, c1, e)
        else:
            C.cx(c1, e)
        carry_into(C, prod, lo[:W], c1)
        C.free([c1])
        # off e, m = hi - 1: the product follows
        times(inverse=True)
        C.x(e)
        dec(C, hi, e)
        C.x(e)
        times()
        a, b = list(prod), lo[:W]
        (c0,) = C.alloc(1)
        carry = [c0] + a[:-1]
        for i in range(W):
            _maj(C, carry[i], b[i], a[i])
        cw = a[W - 1]
        if H:
            C.ccx(pre.top, cw, e)
            pre.write(cw)
        else:
            C.cx(cw, e)
        for i in reversed(range(W)):
            _uma(C, carry[i], b[i], a[i])
        C.free([c0, e])
        times(inverse=True)
        C.free(own + spare)

    def cofactors(self, C, cells, U, V):
        """(U, V) of the cells from (1, 0), semantic: g0: V = 2V + U, swap on b; J: V = 4V;
        B1: V = 2V. The first cell is the terminal one, A0 or J, where V = g0."""
        C.x(U[0])
        C.cx(cells[0][0], V[0])
        for g0, b in cells[1:]:
            V[:] = V[-1:] + V[:-1]
            add(C, U, V, ctl=g0)
            (j,) = C.alloc(1)
            C.x(g0)
            C.ccx(g0, b, j)
            C.x(g0)
            for i in range(len(V) - 1, 0, -1):
                C.cx(V[i - 1], V[i])
                C.ccx(j, V[i], V[i - 1])
                C.cx(V[i - 1], V[i])
            C.x(g0)
            C.ccx(g0, b, j)
            C.x(g0)
            C.free([j])
            (sw,) = C.alloc(1)
            C.ccx(g0, b, sw)
            cswap(C, sw, U, V)
            C.ccx(g0, b, sw)
            C.free([sw])

    def diff(self, C, U, V, r, s, D):
        """D (w + 1 wires, 0 before) = (U*s - V*r)/q, two's complement."""
        k = self.w + 1
        mul_acc(C, U, s[:k], D)
        mul_acc(C, V, r[:k], D, sign=-1)
        mul_odd_const(C, D, pow(self.q, -1, 1 << k))

    def residue(self, C, U, D, Q):
        """D (w + 1 wires, two's complement) -> D mod U on D[:w] in place, its sign kept on D[w];
        Q (w wires, 0 before) takes the quotient bits of |D| by U."""
        w = self.w
        sigma, N = D[w], D[:w]
        for x in N:
            C.cx(sigma, x)
        inc(C, N, sigma)
        for i in range(w - 1, -1, -1):
            # Q[i] = [N >> i >= U]; when it is 1, U < 2^(w - i)
            greater(C, U, N[i:], Q[i])
            C.x(Q[i])
            sub(C, U[:w - i], N[i:], ctl=Q[i])
        (ctl,) = C.alloc(1)
        nz = self._nonzero(C, N)
        C.ccx(sigma, nz, ctl)
        self._nonzero(C, N, nz)
        for x in N:
            C.cx(ctl, x)
        inc(C, N, ctl)
        add(C, U, N, ctl=ctl)
        nz = self._nonzero(C, N)
        C.ccx(sigma, nz, ctl)
        self._nonzero(C, N, nz)
        C.free([ctl])

    def _nonzero(self, C, x, nz=None):
        """A fresh wire holding [x != 0], or with nz given, its clear."""
        for w in x:
            C.x(w)
        chain = []
        z = x[0]
        for w in x[1:]:
            (a,) = C.alloc(1)
            C.ccx(z, w, a)
            chain.append(a)
            z = a
        if nz is None:
            (out,) = C.alloc(1)
            C.cx(z, out)
            C.x(out)
        else:
            C.x(nz)
            C.cx(z, nz)
            C.free([nz])
            out = None
        prev = [x[0]] + chain[:-1]
        for i in reversed(range(len(chain))):
            C.ccx(prev[i], x[i + 1], chain[i])
        C.free(chain)
        for w in x:
            C.x(w)
        return out

    def _mod_double(self, C, U, s, ctl=None):
        """s = 2s mod U (or under ctl), s < U odd on w + 1 wires with top 0. Its flag [2s >= U] is
        then the new s0."""
        (beta,) = C.alloc(1)
        if ctl is None:
            s[:] = s[-1:] + s[:-1]
            greater(C, U, s, beta)
            C.x(beta)
            sub(C, U, s, ctl=beta)
            C.cx(s[0], beta)
        else:
            for i in range(len(s) - 1, 0, -1):
                C.cx(s[i - 1], s[i])
                C.ccx(ctl, s[i], s[i - 1])
                C.cx(s[i - 1], s[i])
            greater(C, U, s, beta, ctl=ctl)
            C.cx(ctl, beta)
            sub(C, U, s, ctl=beta)
            C.ccx(ctl, s[0], beta)
        C.free([beta])

    def _mod_add(self, C, U, r, s, ctl):
        """s = s + ctl*r mod U; the flag is ctl*[s < r] after."""
        (beta,) = C.alloc(1)
        add(C, r, s, ctl=ctl)
        greater(C, U, s, beta)
        C.x(beta)
        sub(C, U, s, ctl=beta)
        (t,) = C.alloc(1)
        greater(C, r, s, t)
        C.ccx(ctl, t, beta)
        greater(C, r, s, t)
        C.free([t, beta])

    def small_pass(self, C, cells, U, r, s):
        """(z, 0) -> (U*z mod U, V*z mod U) = (0, V*z mod U): the cells' semantic tick mod U."""
        w = self.w
        for i in range(w):
            C.ccx(cells[0][0], r[i], s[i])
        for g0, b in cells[1:]:
            self._mod_double(C, U, s)
            self._mod_add(C, U, r, s, g0)
            (j,) = C.alloc(1)
            C.x(g0)
            C.ccx(g0, b, j)
            C.x(g0)
            self._mod_double(C, U, s, ctl=j)
            C.x(g0)
            C.ccx(g0, b, j)
            C.x(g0)
            C.free([j])
            (sw,) = C.alloc(1)
            C.ccx(g0, b, sw)
            cswap(C, sw, r, s[:w])
            C.ccx(g0, b, sw)
            C.free([sw])

    def convert(self, C, term, cells, R, S):
        """R, S (L wires each) from (U*y, V*y) to (U*y mod q, V*y mod q) on their low n wires; the
        wires above are 0 after. term is the first cell's g0, cells the others' (g0, b)."""
        n, w = self.n, self.w
        C.free(R[n + w:] + S[n + w:])
        del R[n + w:], S[n + w:]
        self.fold(C, R)
        self.fold(C, S)
        r, mR, s, mS = R[:n], R[n:n + w], S[:n], S[n:n + w]
        flat = [term] + [x for c in cells for x in c]
        U, V, D = C.alloc(w), C.alloc(w), C.alloc(w + 1)

        def first(flat, r, s, U, V, D):
            cs = [flat[:1] + [None]] + [flat[k:k + 2] for k in range(1, len(flat), 2)]
            self.cofactors(C, cs, U, V)
            self.diff(C, U, V, r, s, D)

        regs1 = [flat, r[:w + 1], s[:w + 1], U, V, D]
        b1 = C.run(first, regs1)
        mul_odd(C, U, mS)
        mul_acc(C, V, mR, mS, sign=-1)
        add(C, D[:w], mS)
        C.free(mS)
        del S[n:]
        Q = C.alloc(w)
        regs2 = [U, D, Q]
        b2 = C.run(lambda U, D, Q: self.residue(C, U, D, Q), regs2)
        sp = C.alloc(w + 1)
        self.small_pass(C, [[term, None]] + cells, U, mR, sp)
        for i in range(w):
            C.cx(D[i], sp[i])
        C.free(sp)
        C.apply(b2, regs2, inverse=True)
        C.free(Q)
        C.apply(b1, regs1, inverse=True)
        C.free(U + V + D + mR)
        del R[n:]


class _Width:
    """two_bit's add-or-shift reads only its word width from self."""

    def __init__(self, n):
        self.n = n
