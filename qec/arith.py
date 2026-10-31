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


def _cuma(C, ctrl, c, b, a):
    C.ccx(c, b, a)
    C.ccx(ctrl, c, b)
    C.cx(a, c)
    # Left ungated because it restores b when ctrl is 0, so _cuma must not become a gated _uma.
    C.cx(a, b)


def _cx(C, ctrl, src, dst):
    if ctrl is None:
        C.cx(src, dst)
    else:
        C.ccx(ctrl, src, dst)


def add(C, a, b, cout=None, ctrl=None):
    """b += ctrl*a mod 2^w (CDKM). cout ^= ctrl*carry out."""
    w = len(a)
    assert len(b) == w
    assert ctrl is None or ctrl not in (*a, *b, cout)
    if w == 0:
        return
    (cin,) = C.alloc(1)
    carries = [cin] + a[:-1]
    for i in range(w):
        _maj(C, carries[i], b[i], a[i])
    if cout is not None:
        _cx(C, ctrl, a[-1], cout)
    for i in reversed(range(w)):
        if ctrl is None:
            _uma(C, carries[i], b[i], a[i])
        else:
            _cuma(C, ctrl, carries[i], b[i], a[i])
    C.free([cin])


def sub(C, a, b, borrow=None, ctrl=None):
    """b -= ctrl*a mod 2^w. borrow ^= ctrl*[a > b]."""
    for w in b:
        C.x(w)
    add(C, a, b, borrow, ctrl)
    for w in b:
        C.x(w)


def gt_into(C, a, b, t, ctrl=None):
    """t ^= ctrl*[a > b], a and b restored. Carry out of (~b + a) is [a > b]."""
    assert ctrl is None or ctrl not in (*a, *b, t)
    w = len(a)
    if w == 0:
        return
    (cin,) = C.alloc(1)
    for x in b:
        C.x(x)
    carries = [cin] + a[:-1]
    for i in range(w):
        _maj(C, carries[i], b[i], a[i])
    _cx(C, ctrl, a[-1], t)
    for i in reversed(range(w)):
        C.ccx(carries[i], b[i], a[i])
        C.cx(a[i], carries[i])
        C.cx(a[i], b[i])
    for x in b:
        C.x(x)
    C.free([cin])


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


class Mod:
    """Arithmetic modulo an odd q < 2^n on canonical values in [0, q)."""

    def __init__(self, q):
        self.q = q
        self.n = q.bit_length()

    def modadd(self, C, a, b, ctrl=None):
        """b = (b + ctrl*a) mod q, for b < q, and a < q when ctrl is 1."""
        q, n = self.q, self.n
        (top,) = C.alloc(1)
        add(C, a, b, top, ctrl)
        ext = b + [top]
        add_const(C, ext, (1 << (n + 1)) - q)
        add_const(C, b, q, ctrl=top)
        # ctrl = 0 leaves top = 1 after the reduction (b < q) and gates off gt_into, so x(top) clears it.
        gt_into(C, a, b, top, ctrl)
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
