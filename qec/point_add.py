"""Exact, complete, controlled fixed-addend point addition:
|c>|R>|0> -> |c>|R + c*A>|0> for every R in E(F_q) encoded as little-endian
affine (x, y) with infinity = (0, 0)."""
from .arith import Mod, cswap, gt_into, load_const, sub
from .builder import Circuit
from .curve import Curve
from .model import TRIT_GROUP, dialog_groups, iteration_bound, trit_width, width_schedule


def _div3_window(C, a, b, c):
    """One restoring-division-by-3 step on three adjacent wires. With
    r' = a + 2b + 4c in [0, 6), it leaves r' - 3q in (a, b) and q = [r' >= 3]
    in c. Exact on those six inputs; 2 CCX, found by exhaustive search."""
    C.ccx(a, b, c)
    C.ccx(a, c, b)
    C.cx(c, a)


def _mul3_window(C, a, b, c):
    """Inverse of _div3_window: (r, q) with r < 3 -> r + 3q in (a, b, c)."""
    C.cx(c, a)
    C.ccx(a, c, b)
    C.ccx(a, b, c)


def pack_trits(C, trits):
    """Pack Euclid dialog pairs into base 3, in place.

    Each pair is (b0, sw) with sw = b0 AND [u > v], so it takes only the
    three values 00, 10, 11. CX(sw, b0) turns it into the trit t = b0 + sw,
    stored as (t & 1, t >> 1). Horner's rule P <- 3P + t then absorbs the
    trits one at a time: with t in the low two wires and P above them, the
    register holds (remainder t, quotient P) of a division of 3P + t by 3, so
    running the division windows backwards, low to high, leaves 3P + t. The
    top window of that division always sees r' < 2 and is skipped, and the
    top wire it would touch is zero and is freed.

    Returns the trit_width(k) wires that hold sum t_j * 3^(k-1-j). On the
    promise domain (every pair valid) every freed wire is zero. Groups whose
    packing frees nothing (k <= 2) are returned raw."""
    k = len(trits)
    if trit_width(k) >= 2 * k:
        return [w for g in trits for w in g]
    P = []
    for j, g in enumerate(trits, 1):
        C.cx(g[1], g[0])
        if j == 1:
            P = list(g)
            continue
        wj = trit_width(j)
        R = list(g) + P + C.alloc(wj - len(P) - 1)
        for i in range(wj - 1):
            _mul3_window(C, R[i], R[i + 1], R[i + 2])
        C.free([R[wj]])
        P = R[:wj]
    return P


def unpack_trits(C, P, k):
    """Exact inverse of pack_trits: k pairs (b0, sw) in dialog order."""
    if trit_width(k) >= 2 * k:
        return [P[2 * i:2 * i + 2] for i in range(k)]
    trits = []
    for j in range(k, 1, -1):
        wj = trit_width(j)
        R = list(P) + C.alloc(1)
        for i in reversed(range(wj - 1)):
            _div3_window(C, R[i], R[i + 1], R[i + 2])
        g, P = R[:2], R[2:]
        # The quotient is below 3^(j-1), so its wires above trit_width(j-1) are zero.
        C.free(P[trit_width(j - 1):])
        P = P[:trit_width(j - 1)]
        C.cx(g[1], g[0])
        trits.append(g)
    C.cx(P[1], P[0])
    trits.append(list(P))
    return trits[::-1]


class DialogMultiplier:
    """M(x, y) = (x, x*y mod q) for x != 0, M(0, y) = (0, y); see model.py.

    The dialog is stored packed: each group of TRIT_GROUP Euclid iterations
    keeps its pairs in trit_width(TRIT_GROUP) = 8 wires instead of 10."""

    def __init__(self, C, mod):
        self.C, self.mod = C, mod
        q, n = mod.q, mod.n
        self.iters = iteration_bound(q)
        self.widths = width_schedule(q, self.iters)
        self.groups = dialog_groups(self.iters)
        self.zbit = next(j for j in range(1, n) if q >> j & 1)
        self._euclid_steps = {}
        self._bezout_step = None
        self.block = None

    def _euclid_step(self, w):
        if w not in self._euclid_steps:
            C = self.C

            def step(U, V, g):
                g0, g1 = g
                C.cx(V[0], g0)
                (t,) = C.alloc(1)
                gt_into(C, U[:w], V[:w], t)
                C.ccx(g0, t, g1)
                gt_into(C, U[:w], V[:w], t)
                C.free([t])
                cswap(C, g1, U[:w], V[:w])
                sub(C, U[:w], V[:w], ctrl=g0)
                V[:] = V[1:] + V[:1]

            n = self.mod.n
            self._euclid_steps[w] = C.record_shape(step, [n, n, 2])
        return self._euclid_steps[w]

    def _bezout(self):
        if self._bezout_step is None:
            C, mod, n = self.C, self.mod, self.mod.n

            def step(g, R, S):
                g0, g1 = g
                mod.moddbl(C, S)
                mod.modadd(C, R, S, ctrl=g0)
                cswap(C, g1, R, S)

            self._bezout_step = C.record_shape(step, [2, n, n])
        return self._bezout_step

    def _body(self, Xr, Yr):
        C, q, n, T = self.C, self.mod.q, self.mod.n, self.iters
        U = C.alloc(n)
        load_const(C, U, q)
        # G starts empty and grows by one packed group per TRIT_GROUP iterations,
        # so the dialog never occupies more than its packed width.
        G = []

        def euclid(U, V, G):
            for s, k in self.groups:
                trits = []
                for i in range(s, s + k):
                    g = C.alloc(2)
                    C.apply(self._euclid_step(self.widths[i]), [U, V, g])
                    trits.append(g)
                G.extend(pack_trits(C, trits))

        E = C.run(euclid, [U, Xr, G])
        assert len(G) == sum(trit_width(k) if trit_width(k) < 2 * k else 2 * k for _, k in self.groups)
        S = C.alloc(n)
        bz = self._bezout()
        offsets, off = [], 0
        for _, k in self.groups:
            w = trit_width(k) if trit_width(k) < 2 * k else 2 * k
            offsets.append((off, w))
            off += w
        # Bezout reads the dialog last iteration first: unpack one group, run
        # its steps, and pack it again before moving to the previous group.
        for (s, k), (off, w) in zip(reversed(self.groups), reversed(offsets)):
            trits = unpack_trits(C, G[off:off + w], k)
            for g in reversed(trits):
                C.apply(bz, [g, Yr, S])
            G[off:off + w] = pack_trits(C, trits)
        cswap(C, U[self.zbit], Yr, S)
        C.apply(E, [U, Xr, G], inverse=True)
        load_const(C, U, q)
        C.free(U + G + Yr)
        Yr[:] = S

    def apply(self, Xr, Yr, inverse=False):
        if self.block is None:
            self.block = self.C.record(self._body, [Xr, Yr])
        self.C.apply(self.block, [Xr, Yr], inverse)


def eq_const_into(C, ctrl, bits, value, target):
    """target ^= ctrl AND [bits == value]."""
    flips = [w for i, w in enumerate(bits) if not value >> i & 1]
    for w in flips:
        C.x(w)
    chain = []
    acc = ctrl
    for w in bits[:-1]:
        (a,) = C.alloc(1)
        C.ccx(acc, w, a)
        chain.append((acc, w, a))
        acc = a
    C.ccx(acc, bits[-1], target)
    for acc, w, a in reversed(chain):
        C.ccx(acc, w, a)
    C.free([a for _, _, a in chain])
    for w in flips:
        C.x(w)


def encode(P, n):
    x, y = (0, 0) if P is None else P
    return x | y << n


class ControlledPointAdd:
    def __init__(self, curve: Curve, A):
        assert A is not None and curve.on_curve(A)
        self.curve, self.A = curve, A
        self.mod = Mod(curve.q)
        self.n = self.mod.n
        neg = curve.neg
        twoA = curve.add(A, A)
        # (exceptional input, its sum with A); the generic kernel divides by zero on these.
        self.exceptions = [(None, A), (A, twoA), (neg(A), None), (neg(twoA), neg(A))]
        enc = [encode(e, self.n) for e, _ in self.exceptions]
        assert len(set(enc)) == 4, "addend order must exceed 3"

    def _mulsub(self, gen, Yr, Xr):
        """X = X - gen * Y^2 mod q, Y untouched."""
        C, mod, n = self.C, self.mod, self.n
        T = C.alloc(n)
        for y, t in zip(Yr, T):
            C.cx(y, t)
        if self._modadd_blk is None:
            self._modadd_blk = C.record_shape(lambda k, a, b: mod.modadd(C, a, b, ctrl=k[0]), [1, n, n])
            self._dbl = C.record_shape(lambda s: mod.moddbl(C, s), [n])
        for i in range(n):
            (k,) = C.alloc(1)
            C.ccx(gen, Yr[i], k)
            C.apply(self._modadd_blk, [[k], T, Xr], inverse=True)
            C.ccx(gen, Yr[i], k)
            C.free([k])
            C.apply(self._dbl, [T])
        for _ in range(n):
            C.apply(self._dbl, [T], inverse=True)
        for y, t in zip(Yr, T):
            C.cx(y, t)
        C.free(T)

    def _kernel(self, gen, Xr, Yr):
        C, mod, q = self.C, self.mod, self.curve.q
        a, b = self.A
        mod.cmodadd_const(C, gen, q - a, Xr)
        mod.cmodadd_const(C, gen, q - b, Yr)
        self.M.apply(Xr, Yr, inverse=True)
        mod.cmodadd_const(C, gen, 3 * a, Xr)
        self._mulsub(gen, Yr, Xr)
        self.M.apply(Xr, Yr)
        mod.cneg(C, gen, Xr)
        mod.cmodadd_const(C, gen, q - b, Yr)
        mod.cmodadd_const(C, gen, a, Xr)

    def build(self, hierarchical=False):
        self.C = C = Circuit(hierarchical)
        n = self.n
        (c,) = C.alloc(1)
        Xr, Yr = C.alloc(n), C.alloc(n)
        declared = Xr + Yr
        self.M = DialogMultiplier(C, self.mod)
        self._modadd_blk = None

        flags = C.alloc(4)
        for f, (e, _) in zip(flags, self.exceptions):
            eq_const_into(C, c, Xr + Yr, encode(e, n), f)
        e, gen = C.alloc(2)
        for f in flags:
            C.cx(f, e)
        C.x(e)
        C.ccx(c, e, gen)
        for f, (e_in, e_out) in zip(flags, self.exceptions):
            load_const(C, Xr + Yr, encode(e_in, n) ^ encode(e_out, n), f)

        self._kernel(gen, Xr, Yr)

        C.ccx(c, e, gen)
        C.x(e)
        for f in flags:
            C.cx(f, e)
        C.free([e, gen])
        for f, (_, e_out) in zip(flags, self.exceptions):
            eq_const_into(C, c, Xr + Yr, encode(e_out, n), f)
        C.free(flags)

        _restore_layout(C, Xr + Yr, declared)
        self.layout = {"control": [c], "x": declared[:n], "y": declared[n:]}
        return C


def _restore_layout(C, cur, declared):
    """SWAP logical bit i from wire cur[i] onto declared[i]. A displaced wire
    holds either another live bit, which is tracked, or a zero ancilla."""
    cur = list(cur)
    where = {w: i for i, w in enumerate(cur)}
    for i, d in enumerate(declared):
        s = cur[i]
        if s == d:
            continue
        C.swap(s, d)
        j = where.pop(d, None)
        del where[s]
        if j is not None:
            cur[j] = s
            where[s] = j
        cur[i] = d
        where[d] = i
