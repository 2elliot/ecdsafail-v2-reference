"""Exact reversible integer and modular arithmetic. Registers are wire lists,
least significant bit first. Every routine returns its ancillas clean on every
input inside its stated domain."""


def _maj(C, c, b, a):
    C.cx(a, b)
    C.cx(a, c)
    C.ccx(c, b, a)


def _uma(C, c, b, a):
    C.ccx(c, b, a)
    C.cx(a, c)
    C.cx(c, b)


def add(C, a, b, cout=None):
    """b += a mod 2^w (CDKM). cout ^= carry out."""
    w = len(a)
    assert len(b) == w
    if w == 0:
        return
    (cin,) = C.alloc(1)
    carries = [cin] + a[:-1]
    for i in range(w):
        _maj(C, carries[i], b[i], a[i])
    if cout is not None:
        C.cx(a[-1], cout)
    for i in reversed(range(w)):
        _uma(C, carries[i], b[i], a[i])
    C.free([cin])


def sub(C, a, b, borrow=None):
    """b -= a mod 2^w. borrow ^= [a > b]."""
    for w in b:
        C.x(w)
    add(C, a, b, borrow)
    for w in b:
        C.x(w)


def gt_into(C, a, b, t):
    """t ^= [a > b], a and b restored. Carry out of (~b + a) is [a > b]."""
    w = len(a)
    if w == 0:
        return
    (cin,) = C.alloc(1)
    for x in b:
        C.x(x)
    carries = [cin] + a[:-1]
    for i in range(w):
        _maj(C, carries[i], b[i], a[i])
    C.cx(a[-1], t)
    for i in reversed(range(w)):
        C.ccx(carries[i], b[i], a[i])
        C.cx(a[i], carries[i])
        C.cx(a[i], b[i])
    for x in b:
        C.x(x)
    C.free([cin])


def and_reg(C, ctrl, a, t):
    for ai, ti in zip(a, t):
        C.ccx(ctrl, ai, ti)


def load_const(C, reg, value, ctrl=None):
    for i, w in enumerate(reg):
        if value >> i & 1:
            if ctrl is None:
                C.x(w)
            else:
                C.cx(ctrl, w)


def add_const(C, reg, value, ctrl=None, cout=None):
    """reg += ctrl*value mod 2^w."""
    value %= 1 << len(reg)
    if value == 0 and cout is None:
        return
    k = C.alloc(len(reg))
    load_const(C, k, value, ctrl)
    add(C, k, reg, cout)
    load_const(C, k, value, ctrl)
    C.free(k)


def cswap(C, ctrl, a, b):
    for ai, bi in zip(a, b):
        C.cx(bi, ai)
        C.ccx(ctrl, ai, bi)
        C.cx(bi, ai)


def csub(C, ctrl, a, b):
    """b -= ctrl*a mod 2^w."""
    t = C.alloc(len(a))
    and_reg(C, ctrl, a, t)
    sub(C, t, b)
    and_reg(C, ctrl, a, t)
    C.free(t)


class Mod:
    """Arithmetic modulo an odd q < 2^n on canonical values in [0, q)."""

    def __init__(self, q):
        self.q = q
        self.n = q.bit_length()

    def modadd(self, C, a, b):
        """b = (a + b) mod q."""
        q, n = self.q, self.n
        (top,) = C.alloc(1)
        add(C, a, b, top)
        ext = b + [top]
        add_const(C, ext, (1 << (n + 1)) - q)
        add_const(C, b, q, ctrl=top)
        gt_into(C, a, b, top)
        C.x(top)
        C.free([top])

    def moddbl(self, C, s):
        """s = 2s mod q, as a relabel plus one conditional subtraction."""
        q, n = self.q, self.n
        (lsb,) = C.alloc(1)
        ext = [lsb] + s[:-1] + [s[-1]]
        (wrap,) = C.alloc(1)
        add_const(C, ext, (1 << (n + 1)) - q, cout=wrap)
        C.x(wrap)
        add_const(C, ext, q, ctrl=wrap)
        C.x(wrap)
        C.cx(ext[0], wrap)
        C.free([wrap, ext[n]])
        s[:] = ext[:n]

    def cmodadd(self, C, ctrl, a, b):
        """b = (b + ctrl*a) mod q."""
        t = C.alloc(len(a))
        and_reg(C, ctrl, a, t)
        self.modadd(C, t, b)
        and_reg(C, ctrl, a, t)
        C.free(t)

    def cmodadd_const(self, C, ctrl, value, b):
        """b = (b + ctrl*value) mod q."""
        value %= self.q
        if value == 0:
            return
        k = C.alloc(self.n)
        load_const(C, k, value, ctrl)
        self.modadd(C, k, b)
        load_const(C, k, value, ctrl)
        C.free(k)

    def cneg(self, C, ctrl, b):
        """b = q - b when ctrl, valid for b in [1, q). A bijection on all of
        {0,1}^n either way, since it is ~b + (q + 1) mod 2^n."""
        for w in b:
            C.cx(ctrl, w)
        add_const(C, b, self.q + 1, ctrl=ctrl)
