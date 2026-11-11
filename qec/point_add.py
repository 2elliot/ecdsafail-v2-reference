"""Exact, complete, controlled fixed-addend point addition:
|c>|R>|0> -> |c>|R + c*A>|0> for every R in E(F_q) encoded as little-endian
affine (x, y) with infinity = (0, 0)."""
from .arith import Mod, _maj, cswap, load_const
from .builder import Circuit, flatten
from .curve import Curve
from .square import SquareFold
from .zero_prefix import ZeroPrefix
from .two_bit import TwoBitMod
from .model import tick_bound, tick_joint, tick_widths
from . import block0_code, rank_one, tail_code, tape_block, two_bit
from .tape_block import FREE, TICKS


class DialogMultiplier:
    """M(x, y) = (x, x*y mod q) for canonical 1 <= x < q; see model.py."""

    # The tape stays decoded: no block, block-zero or tail code, and so no pebbled ticks. Every
    # code call below is gated on blocks or tail, and each encode was paired with a decode that
    # nothing between reads past, so the passes are the same maps on the same states.
    RAW_TAPE = True
    # Chain carries the forward pass leaves as junk for the inverse pass. Each is held through the
    # Bezout pass and saves one CCX in each Euclid pass, so this is the one knob that trades wires
    # for CCX; 1,157 puts the climb stack at 4,096 wires at n = 256.
    CARRY_JUNK = 1157
    # A tick whose registers reach this many bits above its joint threshold runs those bits as one
    # high word, 3 CCX a bit, on a flag that says which register holds it. 6 is 148 CCX cheaper at
    # n = 256; 2 keeps the high word in every toy prime of 8 bits or more.
    HIGH_MIN = 2
    # A high word longer than this is lazy: its lowest LAZY_WINDOW bits keep the borrow chain, and
    # the bits above only move from U to V. One pending borrow wire is enough from 2 up: see
    # _lazy_plan.
    LAZY_WINDOW = 2
    # Lazy ticks run in blocks with a far front: the cells from the front up are one register for
    # U's word and V's, which stays put through the block and lags behind V's shifts; one shift
    # after the block's last tick brings it down by the count. These are the block lengths the
    # plan chooses from; () turns it off. A run of blocks is planned where it saves at least
    # LAG_MIN_SAVING CCX; None plans the first length that fits at every tick, for toy sizes.
    LAG_BLOCK = (15, 7, 3)
    LAG_MIN_SAVING = 1
    # The Bezout pass's first cells run on unreduced integers, up to this many, then one conversion
    # (rank_one.py); 0 turns it off.
    RANK_J = 10

    def __init__(self, C, mod):
        self.C, self.mod = C, mod
        q, n = mod.q, mod.n
        self.iters = tick_bound(q)
        self.widths = tick_widths(q, self.iters)
        self._euclid_steps = {}
        self._realigns = {}
        self._euclid_tick0_block = None
        self._encode_block = None
        self._encode_block0 = None
        self._bezout_step = None
        self._bezout_ticks = {}
        self.narrow = TwoBitMod(q)
        self.block = None
        self._halves = None
        self._cutting = False
        self.carry_junk = self.CARRY_JUNK

    @property
    def pass_tail_width(self):
        return self.widths[self.iters - 1]

    @property
    def tape_external(self):
        return 2 * self.iters - (self.mod.n - self.pass_tail_width) - self.u_lent

    @property
    def blocks(self):
        return 0 if self.RAW_TAPE else self.iters // TICKS

    @property
    def u_lent(self):
        """U's wires that late tape cells take once they are zero, one per zero wire the pass
        ends with: V's tail and each coded block's free wires, which U takes back."""
        coded_free = sum(
            len(self._block_free(j, self.swept)) for j in range(self.blocks))
        return min(self.pass_tail_width + coded_free, self.mod.n - self.pass_tail_width)

    @property
    def lent_after_pass(self):
        coded_free = sum(
            len(self._block_free(j, self.swept)) for j in range(self.blocks))
        return self.pass_tail_width + coded_free - self.u_lent

    @property
    def swept(self):
        """Whether the tape at rest is the B1-swept one, which the code reads last cell first."""
        return self.narrow.merged

    def with_free_slots(self, wires):
        """A coded tape from its wires, with None where each block's free wires were."""
        G = list(wires)
        slots = sorted(slot for j in range(self.blocks) for slot in self._free_slots(j, self.swept))
        for slot in slots:
            G.insert(slot, None)
        return G

    def coded_tape(self, pool):
        """The tape after a pass."""
        return self.with_free_slots(pool[self.lent_after_pass:])

    def _block_free(self, j, swept):
        return block0_code.FREE if swept and j == 0 else FREE

    def _free_slots(self, j, swept):
        return [2 * ((j + 1) * TICKS - 1 - f // 2 if swept else j * TICKS + f // 2)
                + f % 2 for f in self._block_free(j, swept)]

    def _block_wires(self, G, j, swept):
        cells = range(j * TICKS, (j + 1) * TICKS)
        return [G[2 * k + h] for k in (reversed(cells) if swept else cells) for h in (0, 1)]

    def _encoder(self, j, swept):
        if swept and j == 0:
            if self._encode_block0 is None:
                self._encode_block0 = self.C.record_shape(
                    lambda W: block0_code.encode(self.C, W), [2 * TICKS])
            return self._encode_block0
        if self._encode_block is None:
            self._encode_block = self.C.record_shape(lambda W: tape_block.encode(self.C, W), [2 * TICKS])
        return self._encode_block

    def _encode(self, G, j, swept):
        """Code block j. Its free wires leave G, which keeps None in those places."""
        W = self._block_wires(G, j, swept)
        self.C.apply(self._encoder(j, swept), [W])
        wires = [W[f] for f in self._block_free(j, swept)]
        for slot, wire in zip(self._free_slots(j, swept), wires):
            assert G[slot] == wire
            G[slot] = None
        return wires

    def _block0_places(self, U, pool):
        """Where a merged pass leaves block 0's free wires: the last of its zero words, which fill
        U's top u_lent wires and then the pool's first lent_after_pass."""
        n = self.mod.n
        places = [(U, i) for i in range(n - self.u_lent, n)] + [(pool, i) for i in range(self.lent_after_pass)]
        return places[len(places) - len(self._block_free(0, self.swept)):]

    def _code_block0(self, U, pool):
        """Code block 0 of a pass that left it raw, with its free wires at _block0_places."""
        G = self.coded_tape(pool)
        for slot, (reg, i) in zip(self._free_slots(0, self.swept), self._block0_places(U, pool)):
            G[slot] = reg[i]
        self._encode(G, 0, self.swept)

    def _decode(self, G, j, wires, swept):
        assert len(wires) == len(self._block_free(j, swept))
        for slot, wire in zip(self._free_slots(j, swept), wires):
            G[slot] = wire
        self.C.apply(self._encoder(j, swept),
                     [self._block_wires(G, j, swept)], inverse=True)


    def _shared_difference(self, U, V, g0, g1, top=None, junk=None, high=0, held=None, scratch=None,
                           lazy=None):
        """Compare and form the even difference using one majority chain. With top set, bit w - 1
        is done by _top_bit on the domain that mode names. With junk, bits 1 .. len(junk) make
        their carries c_2 .. on those wires and leave them for the inverse pass, at 3 CCX a bit.
        With high, the top `high` bits are done by _high_word."""
        C = self.C
        (cin,) = C.alloc(1)
        chain = len(U) - high if high else len(U) if top is None else len(U) - 1
        m = len(junk) if junk else 0
        assert m <= max(chain - 1, 0)
        carries = [cin, U[0]] + list(junk or []) + U[m + 1:]
        for wire in V:
            C.x(wire)
        for i in range(chain):
            if i == 0:
                C.cx(U[i], V[i])
                C.cx(U[i], carries[i])
                C.cx(V[i], U[i])
            elif i > m:
                _maj(C, carries[i], V[i], U[i])
            else:
                # carries[i + 1] = c ^ (a ^ c)(b ^ c) = MAJ(a, b, c), leaving a ^ c and b ^ c.
                C.cx(carries[i], U[i])
                C.cx(carries[i], V[i])
                C.ccx(U[i], V[i], carries[i + 1])
                C.cx(carries[i], carries[i + 1])
        if high:
            self._high_word(U[chain:], V[chain:], carries[chain], g0, g1, held, scratch, lazy)
        elif top is None:
            C.ccx(g0, carries[chain], g1)
        else:
            self._top_bit(C, top, g0, U[-1], V[-1], carries[chain], g1)
        for i in reversed(range(chain)):
            if 1 <= i <= m:
                # a ^ c and b ^ c: select b ^ g0*(a ^ c); on g1 (so g0) b ^ c ^ (a ^ c) = NOT v ^ a.
                C.ccx(g0, U[i], V[i])
                C.x(V[i])
                C.ccx(g1, V[i], U[i])
                C.x(V[i])
                C.cx(carries[i], U[i])
                C.cx(carries[i], V[i])
                continue
            # Restore a_i and select S_i while unwinding the carry chain.
            if i == 0:
                C.cx(V[i], U[i])
                C.cx(g0, V[i])
            else:
                C.ccx(carries[i], V[i], U[i])
                C.ccx(g0, carries[i], V[i])
            C.cx(U[i], carries[i])
            C.cx(U[i], V[i])
            # S_i XOR c_i XOR 1 is a_i XOR v_i on the swap branch.
            C.cx(carries[i], V[i])
            C.x(V[i])
            if i != 0:
                C.ccx(g1, V[i], U[i])
            C.x(V[i])
            C.cx(carries[i], V[i])
        C.free([cin])
        for wire in V:
            C.x(wire)
        for wire in V[1:]:
            C.cx(g1, wire)

    def _high_word(self, U, V, cw, g0, g1, held, scratch, lazy=None):
        """The chain's bits above the tick's joint threshold, where U or V is zero. held is
        (z, zg): z is 1 when U is the nonzero one and may be either when both are zero, zg = g0*z.
        On entry V holds NOT v and cw the carry c into these bits; g1 is 0. The word X = U ^ V is
        the one that is there. A subtract leaves U = 0 and V = X - e, e = g0*(c ^ z) the borrow
        of the larger minus the smaller, and swaps exactly when z or the borrow leaves the top:
        with z the carries are NOT of the borrows and g1 = g0; without, they are the borrows and
        the top one is g1. scratch holds e and the borrows, and goes back clean.
        With lazy = (s, p) these bits are only the window of a longer word, and s says the word
        has nothing above them. Then the window's top borrow is the word's when s is set, and
        otherwise stays pending on p, the bits above being left as they are."""
        C = self.C
        z, zg = held
        E = scratch[:len(U) + 1]
        for u, v in zip(U, V):
            C.cx(u, v)
            C.x(v)
        for u, v in zip(U, V):
            C.ccx(zg, v, u)
        C.cx(z, cw)
        C.ccx(g0, cw, E[0])
        C.cx(z, cw)
        for i, v in enumerate(V):
            C.x(v)
            C.ccx(E[i], v, E[i + 1])
            C.x(v)
        C.cx(zg, g1)
        if lazy:
            # g1 ^ zg is now the top borrow E s, and p takes the rest of E: E ^ E s.
            s, p = lazy
            C.ccx(E[-1], s, g1)
            C.cx(E[-1], p)
            C.cx(g1, p)
            C.cx(zg, p)
        else:
            C.cx(E[-1], g1)
        for i, v in reversed(list(enumerate(V))):
            C.x(v)
            C.ccx(E[i], v, E[i + 1])
            C.x(v)
            C.cx(E[i], v)
        C.cx(z, cw)
        C.ccx(g0, cw, E[0])
        C.cx(z, cw)
        for u, v in zip(U, V):
            # Off g0 U is as found, and V takes v back. The caller's last X and CX from g1 then
            # leave X ^ borrow ^ top borrow, g1 ^ zg being the top borrow.
            C.cx(u, v)
            C.x(v)
            C.cx(zg, v)

    def _count(self, ctl, reg, scratch, down=False, twice=None):
        """reg += ctl, or reg -= ctl with down, on a two's complement register, low bit first.
        With twice the step is ctl + 2 * twice, for a twice that is never set beside ctl: it and
        the carry out of bit 0 are then never both set, and their XOR is the carry into bit 1."""
        C = self.C
        assert twice is None or len(reg) >= 2
        if down:
            for r in reg:
                C.x(r)
        a = [ctl] + list(scratch[:len(reg) - 1])
        for i in range(len(reg) - 1):
            C.ccx(a[i], reg[i], a[i + 1])
            if i == 0 and twice is not None:
                C.cx(twice, a[1])
        for i in range(len(reg) - 1, 0, -1):
            C.cx(a[i], reg[i])
            if i == 1 and twice is not None:
                C.cx(twice, a[1])
            C.ccx(a[i - 1], reg[i - 1], a[i])
        C.cx(ctl, reg[0])
        if down:
            for r in reg:
                C.x(r)

    def _or_into(self, bits, target, scratch):
        """target ^= OR of bits, through scratch wires that go back clean."""
        C = self.C

        def gate(a, b, t):
            C.x(a)
            C.x(b)
            C.ccx(a, b, t)
            C.x(a)
            C.x(b)
            C.x(t)

        if len(bits) == 1:
            C.cx(bits[0], target)
            return
        dst = list(scratch[:len(bits) - 2]) + [target]
        steps = [(bits[0] if j == 0 else dst[j - 1], bits[j + 1], dst[j]) for j in range(len(bits) - 1)]
        for a, b, t in steps:
            gate(a, b, t)
        for a, b, t in reversed(steps[:-1]):
            gate(a, b, t)

    @staticmethod
    def _top_bit(C, mode, g0, A, B, Cw, H):
        """Bit w - 1 of the shared difference in place of its MAJ, g1 and unwind (5 CCX). On entry
        A = U[w-1] holds a, B = V[w-1] holds b = NOT vhat_top, Cw = U[w-2] holds the carry c into
        bit w - 1, H = g1 is 0. On exit A = U'_top, B = b ^ g0*(a ^ c), Cw = c, H = g0*MAJ(a, b, c),
        as the five CCX leave them, on the domain the mode names:
        top_clear: b = 1 (V's top wire is 0 on entry), so g1 = g0*(a | c) and U'_top = a & ~g0.
        narrow: the step leaves U's top wire 0, so a implies g0 & b, g1 = a ^ g0*c*(a ^ b) and
            B = a ^ b ^ g0*c.
        both: both hold."""
        if mode == "top_clear":
            C.cx(Cw, A)
            C.ccx(g0, A, B)
            C.cx(Cw, A)
            C.ccx(g0, A, H)
            C.cx(H, A)
            C.x(B)
            C.ccx(Cw, B, H)
            C.x(B)
        elif mode == "narrow":
            C.ccx(g0, Cw, H)
            C.cx(A, B)
            C.ccx(H, B, A)
            C.cx(H, B)
            C.ccx(g0, Cw, H)
            C.cx(A, H)
            C.cx(H, A)
        elif mode == "both":
            # P1 = ~a & ~c on B, P2 = ~g0 & ~c on A; then g1 = g0 ^ P1 ^ P2, B = a ^ g0 ^ c ^ P2.
            C.x(B)
            C.x(A)
            C.x(Cw)
            C.ccx(A, Cw, B)
            C.x(A)
            C.x(Cw)
            C.cx(A, B)
            C.x(g0)
            C.x(Cw)
            C.ccx(g0, Cw, A)
            C.x(g0)
            C.x(Cw)
            C.cx(A, H)
            C.cx(B, H)
            C.cx(g0, H)
            C.cx(H, B)
            C.cx(Cw, B)
            C.cx(B, A)
            C.cx(g0, A)
            C.cx(Cw, A)
        else:
            raise ValueError(mode)

    @staticmethod
    def _carry_junk(w, top_clear, narrow, kappa=None):
        """The most carry wires a step at width w can leave: c_2 .. c_chain."""
        if kappa is not None:
            return max(kappa - 1, 0)
        top = w >= 2 and (top_clear or narrow)
        return max(w - (1 if top else 0) - 1, 0)

    def _euclid_step(self, w: int, top_clear: bool = False, narrow: bool = False, last: bool = False,
                     carry: int = 0, kappa=None, high: bool = False, fmode=None, init: bool = False,
                     lazy=None, far=None):
        """top_clear: the tick before ran at the same width, so V is below 2^(w - 1) on entry
        and the second strip's last swap, of two zero bits, never fires.
        narrow: the next width is w - 1, so by E7 this step leaves U below 2^(w - 1).
        last: the pass's last tick, whose input reaches (U, V, E) = (1, 0, 0) by E3.
        carry: the chain's lowest carries, this many, stay on junk wires ahead of the old offset's.
        kappa: the tick's flag threshold, at or above its joint threshold: U and V are never both
            2^kappa or more, and a subtract leaves U below it. The flag is [U >= 2^kappa].
        high: bits kappa and up run as one high word.
        fmode: how the flag goes on to the next tick. same: the next threshold is kappa, the flag
            wire takes flag AND NOT g0 and g0*flag stays as junk. drop: the next threshold is
            kappa - 1, a fresh wire takes that OR U'[kappa - 1] and the old flag stays as junk.
            final: no tick after this reads it.
        init: the flag is made here, from U.
        lazy: (raw, m_in, m_out, leave) for a high word longer than LAZY_WINDOW; see _lazy_plan.
            V then has one more wire on top, and lz holds the pending borrow p and the counter.
        far: (front, c, e, size) for tick e of a lag block of c; see _lag_plan. V then has the
            block's size wires, of which the step relabels only those under the front; the second
            strip and the move of U's word stop under the front; and a junk wire takes
            x = two AND NOT flag, which _realign counts."""
        if last:
            top_clear = narrow = False
            carry = 0
        assert not far or (lazy and not init and not lazy[3] and fmode in ("same", "drop") and far[0] <= w)
        key = (w, top_clear, narrow, last, carry, kappa, high, fmode, init, lazy, far)
        if key not in self._euclid_steps:
            C = self.C

            def step(U, V, E, cell, junk, flag, slot, scratch, lz, lzj):
                g0, b = cell
                if lazy:
                    raw, m_in, m_out, leave = lazy
                    cut = kappa + self.LAZY_WINDOW
                    cnt = lz[1:]
                    spare = scratch[len(cnt) - 1]
                offset = E[0]
                if last:
                    # The input is (1, 0, 0), (1, 1, 0) or (1, 0, 1), and the swept cell is
                    # (V0 ^ E, NOT g0): V[0] = g0 ^ E is cleared and E is left as the pass's junk.
                    C.cx(V[0], g0)
                    C.cx(offset, g0)
                    C.cx(g0, b)
                    C.x(b)
                    C.cx(g0, V[0])
                    C.cx(offset, V[0])
                    V[:] = V[1:] + V[:1]
                    return
                top = None
                if w >= 2 and (top_clear or narrow):
                    top = "both" if top_clear and narrow else "top_clear" if top_clear else "narrow"
                C.cx(offset, V[0])
                C.cx(V[0], g0)
                # A no-subtract tick strips a second bit when V + E is even. The chain leaves U and
                # V as found when g0 = 0, so the strip runs first, on V[1:], and two goes straight
                # onto b, which a subtract leaves 0 until the chain writes g1.
                low = V[1] if w >= 2 else V[0]
                C.cx(offset, low)
                C.x(g0)
                C.x(low)
                C.ccx(g0, low, b)
                C.x(low)
                C.x(g0)
                C.cx(offset, low)
                if lazy:
                    # rho = pi - 1 on max(m_in, m_out) bits, pi the ones of U's and V's bits from
                    # cut up. A bit that the tick before brought to the window's top leaves the
                    # count and takes the pending borrow; a threshold that dropped brings U[cut].
                    if m_in == 0:
                        for r in cnt:
                            C.x(r)
                    for r in cnt[m_in:] if m_in else []:
                        C.cx(cnt[m_in - 1], r)
                    # The second strip takes V[cut] out of the count as well: spare = b * V[cut].
                    C.ccx(b, V[cut], spare)
                    if raw:
                        bit, new, both = V[cut - 1], lzj[0], scratch[len(cnt)]
                        C.ccx(bit, spare, both)
                        C.cx(bit, spare)
                        self._count(spare, cnt, scratch, down=True, twice=both)
                        C.cx(bit, spare)
                        C.ccx(bit, spare, both)
                        C.x(bit)
                        C.ccx(lz[0], bit, new)
                        C.x(bit)
                        C.cx(lz[0], bit)
                        lz[0], lzj[0] = new, lz[0]
                    else:
                        # U[cut] set means U holds the word and V[cut] is 0, so the step is one
                        # up or one down: down is up between two flips of the counter.
                        for r in cnt:
                            C.cx(spare, r)
                        C.cx(U[cut], spare)
                        self._count(spare, cnt, scratch)
                        C.cx(U[cut], spare)
                        for r in cnt:
                            C.cx(spare, r)
                if w >= 2:
                    C.ccx(b, offset, V[1])
                # A lazy V is below 2^(w + 1), or 2^w with top_clear: one more bit than the true one.
                # In a lag block the strip leaves the far cells where they are: V[1], which is 0 on
                # b, goes to the cell under the front, a hole. Off the flag every cell the strip
                # moves lies in V's word; on it V is below 2^kappa and the cells it moves are 0.
                for i in range(1, far[0] - 1 if far else w - 2 + bool(lazy) if top_clear else w - 1 + bool(lazy)):
                    cswap(C, b, [V[i]], [V[i + 1]])
                if far:
                    C.x(flag[0])
                    C.ccx(b, flag[0], lzj[-2])
                    C.x(flag[0])
                if lazy:
                    # The second strip brought one more bit to the window's top, and the pending
                    # borrow, now one place lower, moves back up through it.
                    bit, t = V[cut - 1], lzj[-1]
                    C.ccx(b, bit, spare)
                    C.ccx(b, lz[0], t)
                    C.cx(t, bit)
                    C.x(bit)
                    C.ccx(t, bit, lz[0])
                    C.x(bit)
                    for r in cnt[m_out:]:
                        C.cx(cnt[m_out - 1], r)
                # g1 lands on the junk wire, then goes to b, where the swept tape wants a subtract's
                # own swap; the old offset moves onto the junk wire and E takes g1.
                (g1,) = junk[-1:]
                zg = None
                if init:
                    self._or_into(U[kappa:w], flag[0], scratch)
                if fmode:
                    zg = slot[0] if fmode == "same" else scratch[-1]
                    C.ccx(g0, flag[0], zg)
                if lazy:
                    # pi = 0 is rho's sign. Above the window the word only moves from U to V.
                    self._shared_difference(U[:cut], V[:cut], g0, g1, None, junk[:-1], self.LAZY_WINDOW,
                                            (flag[0], zg), scratch, (cnt[m_out - 1], lz[0]))
                    # In a lag block the far cells hold the word of whichever register has one, so
                    # only the cells under the front change sides.
                    cswap(C, zg, U[cut:far[0] if far else w], V[cut:far[0] if far else w])
                elif high:
                    self._shared_difference(U[:w], V[:w], g0, g1, None, junk[:-1], w - kappa,
                                            (flag[0], zg), scratch)
                else:
                    self._shared_difference(U[:w], V[:w], g0, g1, top, junk[:-1])
                C.cx(g1, b)
                cswap(C, g0, [offset], [g1])
                if fmode == "same":
                    C.cx(zg, flag[0])
                elif fmode == "drop":
                    new, bit = slot[0], U[kappa - 1]
                    C.cx(flag[0], zg)
                    C.cx(bit, new)
                    C.x(bit)
                    C.ccx(zg, bit, new)
                    C.x(bit)
                    C.cx(flag[0], zg)
                    C.ccx(g0, flag[0], zg)
                    flag[0], slot[0] = new, flag[0]
                elif fmode == "final":
                    C.ccx(g0, flag[0], zg)
                if lazy and leave:
                    # The next tick's V fits w - 1 bits. The last bit above the window takes the
                    # pending borrow, and V's extra wire, which held what p left there, p AND NOT
                    # bit, is cleared against it. The counter goes to 0 and p's wire stays as junk.
                    live = cnt[:m_out]
                    self._count(V[cut], live, scratch, down=True)
                    self._count(V[w], live, scratch, down=True)
                    for r in live:
                        C.x(r)
                    C.cx(lz[0], V[cut])
                    C.ccx(lz[0], V[cut], V[w])
                # The far cells stay: the 0 that leaves V's low end comes back in under the front.
                V[:] = V[1:far[0]] + V[:1] + V[far[0]:] if far else V[1:] + V[:1]

            sizes = [w, far[3] if far else w + bool(lazy), 1, 2, 0 if last else 1 + carry, 1 if fmode else 0,
                     1 if fmode in ("same", "drop") else 0, self._step_scratch(key),
                     1 + max(lazy[1:3]) if lazy else 0, lazy[0] + 1 + bool(far) if lazy else 0]
            self._euclid_steps[key] = C.record_shape(step, sizes)
        return self._euclid_steps[key]

    @staticmethod
    def _step_scratch(key):
        """Clean wires a step borrows: the high word's borrows, and g0*flag unless it stays."""
        w, kappa, high, fmode, init, lazy = key[0], key[5], key[6], key[7], key[8], key[9]
        word = w - kappa + 1 if high else 0
        if lazy:
            word = max(DialogMultiplier.LAZY_WINDOW + 1, max(lazy[1:3]) + lazy[0], w - kappa - 2 if init else 0)
        return word + (1 if fmode in ("drop", "final") else 0)

    def _euclid_tick0(self):
        """First Euclid step on U = q, E = 0. U is clean once q's bits flip, so it holds the add's
        carries, which then become the select's x without being cleared."""
        if self._euclid_tick0_block is None:
            C, q, n = self.C, self.mod.q, self.mod.n
            f1 = (1 << n) - q + 1
            j = (f1 & -f1).bit_length() - 1

            def step(U, V, E, cell):
                g0, b = cell
                C.cx(V[0], g0)
                for i in range(n):
                    if q >> i & 1:
                        C.x(U[i])
                # U[i] = g0*c_i, c_i the carry into bit i of x + f1; f1 has no bit below j.
                for i in range(j, n - 1):
                    if i == j:
                        C.ccx(g0, V[i], U[i + 1])
                    elif f1 >> i & 1:
                        # x_i OR c_i under g0 is d ^ x_i*(g0 ^ d), d = g0*c_i.
                        C.cx(U[i], U[i + 1])
                        C.cx(g0, U[i])
                        C.ccx(U[i], V[i], U[i + 1])
                        C.cx(g0, U[i])
                    else:
                        C.ccx(U[i], V[i], U[i + 1])
                for i in range(j, n):
                    if i > j:
                        C.cx(U[i], V[i])
                    if f1 >> i & 1:
                        C.cx(g0, V[i])
                # V = x ^ g0*(f1 ^ c) now, so U ^= V ^ g0*f1 leaves U = x on every lane.
                for i in range(n):
                    C.cx(V[i], U[i])
                    if f1 >> i & 1:
                        C.cx(g0, U[i])
                # Off g0, V still holds x, and U takes q back: U ^= NOT g0 * (x ^ q). Bit 0 of x is g0.
                C.x(g0)
                C.cx(g0, U[0])
                for i in range(1, n):
                    if q >> i & 1:
                        C.x(V[i])
                    C.ccx(g0, V[i], U[i])
                    if q >> i & 1:
                        C.x(V[i])
                C.x(g0)
                for wire in V:
                    C.cx(g0, wire)
                C.cx(g0, E[0])
                V[:] = V[1:] + V[:1]
                C.x(g0)
                C.x(V[0])
                C.ccx(g0, V[0], b)
                C.x(V[0])
                C.x(g0)
                for i in range(n - 2):
                    cswap(C, b, [V[i]], [V[i + 1]])
                # A subtract on U = q > x always swaps, and the swept cell keeps that swap.
                C.cx(g0, b)

            self._euclid_tick0_block = C.record_shape(step, [n, n, 1, 2])
        return self._euclid_tick0_block

    @property
    def tick0_fits(self):
        """Whether tick 0 runs as _euclid_tick0, which needs no scratch."""
        if self.widths[0] != self.mod.n:
            return False
        tick0 = self._euclid_tick0()
        return tick0.wires <= set(tick0.flat_in) and set(tick0.flat_out) == set(tick0.flat_in)

    def _clean_counts(self):
        """Per tick, the clean tape and vacated wires left once the tick has its cell."""
        n, T, W = self.mod.n, self.iters, self.widths
        clean, out = self.tape_external, []
        for i in range(T):
            clean -= 2
            out.append(clean)
            nxt = W[i + 1] if i + 1 < T else self.pass_tail_width
            clean += W[i] - nxt + sum(1 for j in range(nxt, W[i]) if j >= n - self.u_lent)
        return out

    def _tick_plan(self):
        """Per tick, None for tick 0 on U = q, else its step's key (w, top_clear, narrow, last,
        carry, kappa, high, fmode, init, lazy, far). The latest ticks keep their chain carries as junk while
        carry_junk wires last, the earliest of them perhaps only its lowest ones. A tick with
        HIGH_MIN or more bits above its joint threshold, and the clean wires to borrow, runs them
        as a high word; the flag is made at the first such tick and kept up to the last."""
        T, W = self.iters, self.widths
        joint, clean = tick_joint(self.mod.q, T), self._clean_counts()
        plan = []
        for i in range(T):
            w = W[i]
            if i == 0 and self.tick0_fits:
                plan.append(None)
                continue
            # A tick leaves V below 2^(w' - 1), w' its own width.
            top_clear = i >= 1 and w >= 3 and W[i - 1] == w
            # E7: the next tick's U fits its width. E3: the last tick ends at (1, 0, 0).
            last = i == T - 1
            narrow = not last and W[i + 1] == w - 1
            kappa = min(w, joint[i])
            high = not last and kappa >= 1 and w - kappa >= self.HIGH_MIN and w - kappa + 2 <= clean[i]
            plan.append([w, top_clear, narrow, last, 0, kappa, high, None, False, None, None])
        flagged = [i for i, key in enumerate(plan) if key and key[6]]
        if flagged:
            first, final = flagged[0], flagged[-1]
            plan[first][8] = True
            for i in range(first, final + 1):
                step = plan[i][5] - plan[i + 1][5] if i < final else None
                assert step in (None, 0, 1) and clean[i] >= 1
                plan[i][7] = "final" if i == final else "drop" if step else "same"
        self._lazy_plan(plan)
        self._lag_plan(plan)
        for i, key in enumerate(plan):
            if key and not key[7]:
                key[5] = None
        budget = self.carry_junk
        for key in reversed(plan):
            if key is None or key[3]:
                continue
            key[4] = min(budget, self._carry_junk(key[0], key[1], key[2], key[5] if key[6] else None))
            budget -= key[4]
            if not budget:
                break
        return [None if key is None else tuple(key) for key in plan]

    def _lazy_plan(self, plan):
        """Mark the lazy ticks: the high ticks with more than LAZY_WINDOW bits above the threshold.
        There V's physical word is the true one plus p * 2^cut, cut = kappa + LAZY_WINDOW, p the
        one pending borrow, so V has one more wire. A raw tick finds p one place lower, at
        2^(cut - 1), and the count taken over V from cut - 1 up, until its first crossing puts both
        right. Each gets (raw, m_in, m_out, leave):
        raw: the tick before kept its threshold, so the bit it brought to cut - 1 has not met p.
        m_in, m_out: the counter's width before and after. It holds rho = pi - 1, pi the ones of
            U and V from cut up, which number at most w - kappa - LAZY_WINDOW + 1.
        leave: the run's last tick.
        One p is enough because the window is all ones under a pending p at every subtract,
        which comes after the tick's crossing bits have met p. This holds along a run from the
        borrow that set p, not for an arbitrary p. A borrow leaves
        the window only when it is 0, and leaves it all ones over a low zone 2^kappa - d with
        d <= u + 2. Then v >= 3 * 2^kappa, so u < 2^kappa / 9 by the tick's potential bound, the
        low zone stays above u, and no subtract borrows from the window; each shift halves d,
        each subtract takes it to (d + u) / 2, and a bit that comes to the window's top under p
        is a 1 or clears p. The run starts where the word first has LAZY_WINDOW + 1 bits, after
        a drop at the same width, so the one bit above the window is new and V's is 0. It ends
        at that length before a narrower tick with the same threshold, where the last tick clears
        the counter and V's extra wire and leaves p's wire as junk; the plan is left alone
        otherwise."""
        R, W = self.LAZY_WINDOW, self.widths
        lazy = [i for i, key in enumerate(plan) if key and key[6] and key[0] - key[5] > R]
        if not lazy:
            return
        first, final = lazy[0], lazy[-1]
        if lazy != list(range(first, final + 1)) or plan[first - 1] is None or final + 1 >= len(plan):
            return
        before, after = plan[first - 1], plan[final + 1]
        if (before[0], before[5]) != (W[first], plan[first][5] + 1) or W[first] - plan[first][5] != R + 1:
            return
        if (after[0], after[5]) != (W[final] - 1, plan[final][5]) or W[final] - plan[final][5] != R + 1:
            return
        m_in = 0
        for i in lazy:
            w, kappa = plan[i][0], plan[i][5]
            m_out = (w - kappa - R).bit_length() + 1
            plan[i][9] = (i > first and plan[i - 1][5] == kappa, m_in, m_out, i == final)
            m_in = m_out

    @staticmethod
    def _lag_counts(c):
        """For a block of c ticks: the adders that count its strips, and per shift stage j the cells
        under the front whose swap waits on the count's lower bits."""
        J, sizes, adders = c.bit_length(), [c], 0
        for j in range(J):
            carries = 0
            while sizes[j] >= 2:
                sizes[j] -= 2 if sizes[j] >= 3 else 1
                carries += 1
            adders += carries
            sizes.append(carries)
        assert sizes[J] == 0 and all(sizes[:J])
        return adders, [min((1 << j) - 1, c - (1 << j)) if j == J - 1 else (1 << j) - 1 for j in range(J)]

    def _lag_cost(self, c, cells):
        """CCX of _realign on c count wires with this many cells from the front up."""
        adders, under = self._lag_counts(c)
        return 2 * (adders + sum(under)) + sum(under) + c.bit_length() * cells

    def _lag_plan(self, plan):
        """Mark the lag blocks: runs of lazy ticks, none the lazy run's first or last, of a length
        from LAG_BLOCK. Through a block of c ticks with front F the cells from F up are one
        register X for U's word and V's (the joint domain leaves at most one of them a bit there),
        and no gate of the block's ticks touches it. Under the front V shifts as ever: the relabel
        brings a 0 in under the front each tick, and so does each second strip that counts,
        x = two AND NOT flag. While U holds the word (the flag is set when the tick ends) X holds
        U's bits from F up, U's own wires being 0 there, and a 0 under the front is a 0 of V, which
        is below 2^kappa. A subtract that takes U's word to V moves only the cells under the
        front, and from then on X is V's and lags: after lag = sum of (NOT flag left by the tick)
        + sum of x, a cell i >= F holds V's bit i - lag and the cells [F - lag, F) are 0. _realign
        closes them after the last tick. F is the least front with F - 2 e >= cut_e + 2 for every
        tick e but the last, so that with at most 2 e holes before the tick, and one more after
        its strip, the cells it reads and the one its strip brings down, up to cut + 1, are V's
        own; the last tick needs only what it reads, F - 2 e >= cut_e + 1. That also
        keeps a U whose flag came on at tick e, which is below 2^kappa_(e-1), under the at most
        2 e holes counted by then, so the move lands under them. F is at most the width after the
        block. A run of blocks starts and ends by swapping the far cells of U and V on the flag;
        between two blocks the cells between their fronts are swapped on it. Each tick gets
        (front, c, e, size), size the wires V has when the block starts. The run is the choice
        that saves the most."""
        lengths = self.LAG_BLOCK if isinstance(self.LAG_BLOCK, tuple) else (self.LAG_BLOCK,)
        R, W, T = self.LAZY_WINDOW, self.widths, len(plan)

        def inner(i):
            return 0 <= i < T and plan[i] and plan[i][9] and plan[i][9][1] != 0 and not plan[i][9][3]

        def block(i, c):
            """(CCX saved, front) of a block of c ticks from tick i, or None."""
            ticks = range(i, i + c)
            if not c or not all(inner(t) for t in ticks):
                return None
            kappa = [plan[t][5] for t in ticks]
            front = max(kappa[e] + R + 2 * e + 2 - (e == c - 1) for e in range(c))
            if front > W[i + c] or front < 2 * c:
                return None
            eager = sum(2 * plan[t][0] - 1 - plan[t][1] - plan[t][5] - R for t in ticks)
            mine = sum(2 * front - 2 - kappa[e] - R for e in range(c)) + 3 * c + self._lag_cost(2 * c, W[i] + 1 - front)
            return eager - mine, front

        choice = [None] * T
        if self.LAG_MIN_SAVING is None:
            i = 0
            while i < T:
                found = next(((c, block(i, c)) for c in lengths if block(i, c)), None)
                if found:
                    choice[i] = (found[0], found[1][1])
                i += found[0] if found else 1
        else:
            # best[i]: the most a run that ends before tick i saves, its entry swap paid.
            best, back = [None] * (T + 1), [None] * (T + 1)
            for i in range(T):
                for c in lengths:
                    found = block(i, c)
                    if found is None:
                        continue
                    saved, front = found
                    start, prev = -(W[i] - front), None
                    if best[i] is not None and best[i] - abs(back[i][2] - front) > start:
                        start, prev = best[i] - abs(back[i][2] - front), i
                    if best[i + c] is None or start + saved > best[i + c]:
                        best[i + c], back[i + c] = start + saved, (i, c, front, prev)
            ends = [(best[i] - (W[i] - back[i][2]), i) for i in range(T + 1) if best[i] is not None]
            if ends and max(ends)[0] >= self.LAG_MIN_SAVING:
                i = max(ends)[1]
                while i is not None:
                    start, c, front, i = back[i]
                    choice[start] = (c, front)
        for i, found in enumerate(choice):
            if found:
                c, front = found
                for e in range(c):
                    plan[i + e][10] = (front, c, e, W[i] + 1)

    def _realign(self, c, cells):
        """Close a lag block's holes. V is the cells from c below the front up, xs the block's c
        count bits, k their sum: the cells [front - k, front) are 0, and every cell from the front up
        moves down k places, zeros coming in on top. The count is made by adders on scratch, and
        the shift goes low bit first: at stage j the cells from front - r up, r = k mod 2^j, move
        down 2^j, so the swaps of the 2^j - 1 cells under the front wait on k_j AND [r >= m], which
        are the upper half of the next stage's thresholds [k mod 2^(j+1) >= 2^j + m]."""
        key = (c, cells)
        if key not in self._realigns:
            C = self.C
            J = c.bit_length()
            adders, under = self._lag_counts(c)

            def body(V, xs, scratch):
                undo, free = [], list(scratch)

                def cx(a, b):
                    C.cx(a, b)
                    undo.append((a, b))

                def ccx(a, b, t):
                    C.ccx(a, b, t)
                    undo.append((a, b, t))

                buckets, k = [list(xs)] + [[] for _ in range(J)], []
                for j in range(J):
                    ws = buckets[j]
                    while len(ws) >= 2:
                        s = free.pop(0)
                        if len(ws) >= 3:
                            # s = MAJ(a, b, d) = a ^ (a ^ b)(a ^ d), and d takes the sum.
                            a, b, d = ws[:3]
                            del ws[:3]
                            cx(a, b)
                            cx(a, d)
                            ccx(b, d, s)
                            cx(a, s)
                            cx(b, d)
                            cx(a, d)
                            ws.append(d)
                        else:
                            a, b = ws[:2]
                            del ws[:2]
                            ccx(a, b, s)
                            cx(a, b)
                            ws.append(b)
                        buckets[j + 1].append(s)
                    k.append(ws[0])
                front, th = c, []
                for j in range(J):
                    prods = []
                    for m in range(under[j]):
                        prods.append(free.pop(0))
                        ccx(k[j], th[m], prods[m])
                    for i in range(front - under[j], len(V)):
                        cswap(C, prods[front - i - 1] if i < front else k[j], [V[i - (1 << j)]], [V[i]])
                    if j < J - 1:
                        # [r >= m] becomes k_j OR [r >= m] = k_j ^ [r >= m] ^ k_j [r >= m].
                        for m in range(under[j]):
                            cx(k[j], th[m])
                            cx(prods[m], th[m])
                        th = th + [k[j]] + prods
                for gate in reversed(undo):
                    (C.cx if len(gate) == 2 else C.ccx)(*gate)

            self._realigns[key] = C.record_shape(body, [c + cells, c, adders + sum(under)])
        return self._realigns[key]

    def _far_swap(self, flag, U, V, lo, hi):
        """On the flag, U's cells [lo, hi) change places with V's, which are 0 there."""
        for j in range(lo, hi):
            cswap(self.C, flag, [U[j]], [V[j]])

    def _lazy_width(self, plan):
        """The widest the lazy counter gets."""
        return max((key[9][2] for key in plan if key and key[9]), default=0)

    def _junk_sizes(self):
        """Junk wires each tick leaves for the inverse pass: its carries, the flag it makes, the
        wire its flag update leaves, a lazy tick's, and the retired offset. A lazy tick leaves
        the wire that held p if the tick is raw, its counted strip if it is in a lag block, and
        b * p; the first one brings V's extra wire, p and the counter, which
        the last one clears but for p."""
        plan = self._tick_plan()
        m = self._lazy_width(plan)

        def lazy(key):
            return (2 + m if key[9][1] == 0 else 0) + key[9][0] + 1 + bool(key[10]) if key[9] else 0

        return [0 if key is None or key[3] else 1 + key[4] + (key[7] in ("same", "drop")) + key[8] + lazy(key)
                for key in plan]

    @property
    def junk_width(self):
        """The pass's junk: E, then each tick's."""
        return 1 + sum(self._junk_sizes())

    def _euclid_pass(self, U: list[int], pool: list[int], junk: list[int], block0_raw: bool = False) -> None:
        # The pool starts as V plus clean external tape. It ends as the zero wires U did not
        # take back, then tape; U ends as its low wires and those zero wires. The boundary is a
        # permutation of U and the pool together. junk starts zero and ends holding each tick's
        # junk, which only the inverse pass clears.
        n = self.mod.n
        V, available, G = pool[:n], pool[n:], []
        lend = U[n - self.u_lent:]
        E = junk[:1]
        plan, sizes, used, flag = self._tick_plan(), self._junk_sizes(), 1, []
        counter, ext, pending = self._lazy_width(plan), None, []
        strips, after, whole = [], [], False
        assert len(junk) == 1 + sum(sizes)
        for i in range(self.iters):
            w = self.widths[i]
            assert (whole or len(V) == w) and len(available) >= 2
            cell, available = available[:2], available[2:]
            if plan[i] is None:
                self.C.apply(self._euclid_tick0(), [U, V, E, cell])
            else:
                key = plan[i]
                tick_junk = junk[used:used + sizes[i]]
                used += sizes[i]
                # The tick's junk: its carries, the flag it makes, its flag update's wire, g1's.
                extra = tick_junk[key[4]:-1]
                if key[8]:
                    flag = [extra.pop(0)]
                held = flag if key[7] else []
                slot = [extra.pop(0)] if key[7] in ("same", "drop") else []
                lz = []
                if key[9]:
                    if key[9][1] == 0:
                        ext, pending = junk.index(extra[0]), [extra[1]]
                        count = extra[2:2 + counter]
                        V.append(extra[0])
                        del extra[:2 + counter]
                    elif not whole:
                        V.append(spare)
                    lz = pending + count[:max(key[9][1:3])]
                if key[10]:
                    front, c, e, size = key[10]
                    if e == 0:
                        assert len(V) == size
                        if not (plan[i - 1] and plan[i - 1][10]):
                            self._far_swap(held[0], U, V, front, w)
                        # V keeps every wire it has to the block's end.
                        whole = True
                    strips.append(extra[-2])
                scratch = (available + junk[used:])[:self._step_scratch(key)]
                assert len(scratch) == self._step_scratch(key)
                self.C.apply(self._euclid_step(*key),
                             [U[:w], V, E, cell, tick_junk[:key[4]] + tick_junk[-1:], held, slot, scratch,
                              lz, extra if key[9] else []])
                if key[7]:
                    flag = held
                if key[10]:
                    # What the flag update left: g0 * flag on a same tick, the flag it found on a drop.
                    after.append((key[7], slot[0]))
                if key[9]:
                    pending = lz[:1]
            G.extend(cell)
            next_width = self.widths[i + 1] if i + 1 < self.iters else self.pass_tail_width
            assert w - 1 <= next_width <= w
            if plan[i] is not None and plan[i][9]:
                far = plan[i][10]
                if far and far[2] == far[1] - 1:
                    front, c = far[:2]
                    cells = V[front - 2 * c:]
                    adders, under = self._lag_counts(2 * c)
                    scratch = (available + junk[used:])[:adders + sum(under)]
                    assert len(scratch) == adders + sum(under)
                    # The lag: a counted strip, and a tick that leaves the flag off. The flag tick e
                    # leaves is the one tick e + 1 finds: a drop keeps it on its old wire, and a same
                    # tick's is the flag it leaves plus its g0 * flag, which are never both set.
                    left = list(self.C.alloc(c))
                    for e in range(c):
                        j = e + 1
                        while j < c and after[j][0] == "same":
                            self.C.cx(after[j][1], left[e])
                            j += 1
                        self.C.cx(after[j][1] if j < c else flag[0], left[e])
                        self.C.x(left[e])
                    self.C.apply(self._realign(2 * c, len(cells) - 2 * c), [cells, strips + left, scratch])
                    for e in range(c):
                        j = e + 1
                        self.C.x(left[e])
                        while j < c and after[j][0] == "same":
                            self.C.cx(after[j][1], left[e])
                            j += 1
                        self.C.cx(after[j][1] if j < c else flag[0], left[e])
                    self.C.free(left)
                    strips, after, whole = [], [], False
                if not whole:
                    # V's extra wire stays above the next width; the run's last tick leaves it clear,
                    # and it takes the place in the junk that the first one's came from.
                    spare = V.pop(next_width)
                    if plan[i][9][3]:
                        junk[ext] = spare
                    if far:
                        # U's word under the next front, or all of it after the run's last block,
                        # goes back to U's own cells.
                        nxt = plan[i + 1] and plan[i + 1][10]
                        lo, hi = sorted((front, nxt[0])) if nxt else (front, next_width)
                        self._far_swap(flag[0], U, V, lo, hi)
            if not whole:
                available.extend(V[next_width:])
                V[:] = V[:next_width]
            # By L6, u < 2**next_width from here on, and later steps touch only U[:w] for
            # w <= next_width, so U's wires from next_width up are zero for the rest of the pass.
            available.extend(reversed([u for u in U[next_width:w] if u in lend]))
        assert len(V) == self.pass_tail_width and not available
        # Each subtract cell already holds its own swap. On the merged path the B1 sweep runs
        # down the tape, and a block is coded only after the B1 step that reads its bottom cell.
        T, swept, clean = self.iters, self.swept, []
        both = None

        def mark(k):
            self.C.x(G[2 * k])
            self.C.x(G[2 * k - 2])
            self.C.ccx(G[2 * k], G[2 * k - 2], both)
            self.C.x(G[2 * k - 2])
            self.C.x(G[2 * k])

        # The pass reaches (1, 0) within T ticks, and each no-subtract cell after its last subtract
        # has b = 1. So when cells T - 1 and T - 2 are both no-subtract their b are equal: the B1
        # step between them never swaps, and mark(T - 2) is made at once.
        late = swept and T >= 3
        if swept and T >= 2:
            (both,) = self.C.alloc(1)
            mark(T - 2 if late else T - 1)
        for i in range(T - 1, -1, -1):
            if swept and i + 1 < T and not (late and i == T - 2):
                k = i + 1
                self._b1_step(G, k, both)
                if k >= 2:
                    # Turn mark(k,k-1) into mark(k-1,k-2) before any encode.
                    self.C.cx(G[2 * k - 4], G[2 * k])
                    self.C.x(G[2 * k - 2])
                    self.C.ccx(G[2 * k - 2], G[2 * k], both)
                    self.C.x(G[2 * k - 2])
                    self.C.cx(G[2 * k - 4], G[2 * k])
                else:
                    mark(1)
            done = i + 1 if swept else i
            if done % TICKS == 0 and done // TICKS < self.blocks:
                clean.extend(self._encode(G, done // TICKS, swept))
        if both is not None:
            self.C.free([both])
        # Block 0's free wires take their places among the zero words raw, and block 0 is coded
        # last, unless block0_raw leaves that to the caller.
        raw = []
        if swept and self.blocks:
            raw = [G[slot] for slot in self._free_slots(0, swept)]
            for slot in self._free_slots(0, swept):
                G[slot] = None
        zeros = V + clean + raw
        U[:] = U[:n - self.u_lent] + zeros[:self.u_lent]
        pool[:] = zeros[self.u_lent:] + [w for w in G if w is not None]
        if raw and not block0_raw:
            self._code_block0(U, pool)

    def _bezout(self):
        if self._bezout_step is None:
            C, mod = self.C, self.mod

            def step(g, R, S):
                self.narrow.step(C, g, R, S)

            self._bezout_step = C.record_shape(step, [2, mod.n, mod.n])
        return self._bezout_step

    def _bezout_tick(self, first, last, pebble=False, after=False):
        """after marks a first tick that follows the terminal cell, where S is 0 or R. It takes
        that cell's g0 as a last argument."""
        if (first, last, pebble, after) not in self._bezout_ticks:
            C, n = self.C, self.mod.n

            def tick(prev, g, nxt, R, S, flag, tau):
                # flag comes in holding [previous cell is B1] and H's two top bits, or 0 on a
                # first tick, which has no previous cell, and goes out holding this cell's, or
                # empty after the last.
                out = self.narrow.tick(C, prev or None, g, nxt[0] if nxt else None, R, S,
                                       pebble=pebble, held=flag, cell0=last,
                                       tau=tau[0] if after else None)
                flag[:] = [] if out is None else out

            sizes = [0 if first else 2, 2, 0 if last else 1, n, n, 3, int(after)]
            self._bezout_ticks[first, last, pebble, after] = C.record_shape(tick, sizes)
        return self._bezout_ticks[first, last, pebble, after]

    def _b1_step(self, G, i, both):
        """Use the pass-wide mark while both cells are still raw."""
        cswap(self.C, both, [G[2 * i + 1]], [G[2 * i - 1]])

    def _bezout_pass(self, G, Yr, S, block0_raw=False):
        """G is the coded tape. Each block is decoded on a clean wire before a tick reads it and
        coded again once no tick reads it, except block 0 on the merged path with block0_raw."""
        C, T, swept = self.C, self.iters, self.swept
        if not self.narrow.merged:
            bz = self._bezout()
            for i in reversed(range(T)):
                j = i // TICKS
                if j < self.blocks and i % TICKS == TICKS - 1:
                    self._decode(G, j, C.alloc(len(self._block_free(j, swept))), swept)
                C.apply(bz, [G[2 * i:2 * i + 2], Yr, S])
                if j < self.blocks and i % TICKS == 0:
                    C.free(self._encode(G, j, swept))
            return
        # A tick reads the previous cell and writes the next cell's g0, so a block is decoded
        # one tick before its first tick and coded one tick after its last. Across each edge
        # two blocks are decoded for two ticks.
        if self.blocks * TICKS == T:
            self._decode(G, self.blocks - 1,
                         C.alloc(len(self._block_free(self.blocks - 1, swept))), swept)
            if self.blocks == 1:
                block0_code.hold(C, self._block_wires(G, 0, swept))
                C.free([G[1]])
                G[1] = None
        # The final swept cell is (1, 0) or (0, 1): the terminal subtract
        # is followed only by no-subtract cells. Lend its b only when uncoded.
        lend_tail = T > 1 and self.blocks * TICKS < T
        tail = (not self.RAW_TAPE and lend_tail and T >= tail_code.CELLS + 1
                and T - self.blocks * TICKS >= tail_code.CELLS
                and tick_widths(self.mod.q, T)[T - tail_code.CELLS] <= 6)
        lo = 2 * (T - tail_code.CELLS)
        tail_encoded = False
        flag = []
        J = self._rank_j(lend_tail, tail)
        if J:
            self._rank_start(G, Yr, S, J)
        for i in reversed(range(T - J)):
            nxt_block = i // TICKS - 1
            if i % TICKS == 0 and 0 <= nxt_block < self.blocks:
                self._decode(G, nxt_block,
                             C.alloc(len(self._block_free(nxt_block, swept))), swept)
                if nxt_block == 0:
                    block0_code.hold(C, self._block_wires(G, 0, swept))
                    C.free([G[1]])
                    G[1] = None
            if i == 1 and self.blocks:
                G[1] = C.alloc(1)[0]
                block0_code.hold(C, self._block_wires(G, 0, swept), inverse=True)
            no_prev = lend_tail and i == T - 2 and not J
            if no_prev:
                # The terminal cell is never B1; release b before its last reader.
                C.cx(G[-2], G[-1])
                C.x(G[-1])
                C.free([G[-1]])
                G[-1] = None
            first, last = i == T - 1 or no_prev or (J and i == T - J - 1), i == 0
            edge = any(i in (j * TICKS, j * TICKS - 1) for j in range(1, self.blocks))
            pebble = edge or (tail and self.blocks > 0 and i == self.blocks * TICKS)
            pebble = pebble or (self.blocks > 0 and i in (0, 1))
            prev = [] if first else G[2 * i + 2:2 * i + 4]
            nxt = [] if last else G[2 * i - 2:2 * i - 1]
            if i == T - 1:
                # S starts at zero and the terminal cell is A0, A1 or J, never B1.
                for r, s in zip(Yr, S):
                    C.ccx(G[2 * i], r, s)
            else:
                if not flag:
                    # The pass's first tick call: no tick before it handed a flag over, and the
                    # cell before, if any, is the terminal cell, never B1, so H's bits are 0.
                    flag = C.alloc(3)
                    if not first:
                        self.narrow._mark_b1(C, prev, flag[0])
                after = first and not last and i == T - 2
                if last and self._cutting:
                    self._cell0_cut(prev, G[:2], Yr, S, flag, pebble)
                else:
                    C.apply(self._bezout_tick(first, last, pebble, after),
                            [prev, G[2 * i:2 * i + 2], nxt, Yr, S, flag,
                             [G[2 * T - 2]] if after else []])
            prev_block = (i + 1) // TICKS
            if (i + 1) % TICKS == 0 and i + 1 < T and prev_block < self.blocks:
                C.free(self._encode(G, prev_block, swept))
            if tail and i == T - tail_code.CELLS - 1:
                W = G[lo:lo + tail_code.LIVE]
                tail_code.apply(C, W)
                tail_encoded = True
                C.free([W[tail_code.FREE[0]]])
                G[lo + tail_code.FREE[0]] = None
        assert not flag
        if self.blocks and not block0_raw:
            C.free(self._encode(G, 0, swept))
        if tail:
            assert tail_encoded, "tail inverse requires its forward apply"
            G[lo + tail_code.FREE[0]] = C.alloc(1)[0]
            tail_code.apply(C, G[lo:lo + tail_code.LIVE], inverse=True)
        if lend_tail:
            G[-1] = C.alloc(1)[0]
            C.x(G[-1])
            C.cx(G[-2], G[-1])

    def _rank_j(self, lend_tail, tail):
        if not (self.RANK_J and self.RAW_TAPE and lend_tail and not tail and not self.blocks):
            return 0
        return rank_one.plan(self.mod.q, self.iters, self.widths, self.RANK_J)

    def _rank_start(self, G, Yr, S, J):
        """Cells T - 1 down to T - J on unreduced integers, then one conversion. The tick after it
        is a first tick: the last cell here is semantic, so no lane is held."""
        C, T, n = self.C, self.iters, self.mod.n
        ro = rank_one.RankOne(self.mod.q, J, self.widths[T - J])
        R, Sw = Yr + C.alloc(ro.L - n), S + C.alloc(ro.L - n)
        term = G[2 * T - 2]
        for r, s in zip(Yr, S):
            C.ccx(term, r, s)
        # The terminal cell is never B1; release b before its last reader.
        C.cx(G[-2], G[-1])
        C.x(G[-1])
        C.free([G[-1]])
        G[-1] = None
        p = None
        for i in range(T - 2, T - J, -1):
            p = ro.tick(C, G[2 * i:2 * i + 2], G[2 * i + 2:2 * i + 4], p, R, Sw)
        i = T - J
        ro.last_tick(C, G[2 * i:2 * i + 2], G[2 * i + 2:2 * i + 4], p, R, Sw)
        cells = [G[2 * k:2 * k + 2] for k in range(T - 2, T - J - 1, -1)]
        ro.convert(C, term, cells, R, Sw)
        Yr[:], S[:] = R, Sw

    def _body(self, Xr, Yr):
        C, q, n, T = self.C, self.mod.q, self.mod.n, self.iters
        U = C.alloc(n)
        load_const(C, U, q)
        pool = list(Xr) + C.alloc(self.tape_external)

        # On the merged path block 0 stays decoded from Bezout's end into the inverse pass, so
        # the pass is recorded without its last encode, which runs here instead.
        raw0 = self.swept and self.blocks > 0
        junk = C.alloc(self.junk_width)
        E = C.run(lambda U, pool, junk: self._euclid_pass(U, pool, junk, raw0), [U, pool, junk])
        if raw0:
            self._code_block0(U, pool)
        lent = self.lent_after_pass
        G = self.coded_tape(pool)
        # Terminal V and the wires the code freed are zero; lend them during Bezout.
        C.free(pool[:lent])
        # U is (1,0,...,0). Reuse its entire word for the reverse coefficient;
        # the surviving V wire stays zero for the inverse Euclid boundary.
        C.x(U[0])
        S = U
        self._bezout_pass(G, Yr, S, raw0)
        if raw0:
            raw = [G[slot] for slot in self._free_slots(0, True)]
            for slot in self._free_slots(0, True):
                G[slot] = None
        # Restore the tape boundary using the newly allocated tail wire.
        pool[lent:] = [w for w in G if w is not None]
        C.x(Yr[0])
        pool[:lent] = C.alloc(lent)
        if raw0:
            # Block 0's free wires take the places of the last zero words, which are released.
            places = self._block0_places(Yr, pool)
            C.free([reg[i] for reg, i in places])
            for (reg, i), w in zip(places, raw):
                reg[i] = w
        C.apply(E, [Yr, pool, junk], inverse=True)
        load_const(C, Yr, q)
        Xr[:] = pool[:n]
        C.free(Yr + pool[n:] + junk)
        Yr[:] = S

    def _cell0_cut(self, prev, g, R, S, flag, pebble):
        """Cell 0's tick, run in the body itself and cut at its zero test. That test is the tick's
        one gt_into call, t ^= z*[S >> 1 > 0] against R's clean wires on the pass's final S; here
        it emits nothing, and split()'s caller emits its own in the cut. two_bit.py stays as it is."""
        C = self.C

        def zero_test(C, a, b, t, ctrl=None):
            assert a == S[:-1] and self._test is None
            self._test = [S[-1:] + a, list(b), [t], [ctrl]]
            C.cut()

        saved, two_bit.gt_into = two_bit.gt_into, zero_test
        try:
            out = self.narrow.tick(C, prev or None, g, None, R, S, pebble=pebble, held=flag, cell0=True)
        finally:
            two_bit.gt_into = saved
        assert out is None
        flag[:] = []

    def split(self, Xr, Yr):
        """Whether the body splits into two blocks around cell 0's zero test, so that apply can
        take a test. It does not on a pass with no cell 0 tick."""
        if self._halves is None:
            self._cutting, self._test = True, None
            self._halves = self.C.record_cut(self._body, [Xr, Yr]) or ()
            self._cutting = False
            if self._halves:
                # Layout wires the second block never touches and ends on: clean at the cut and
                # idle after it.
                before, after = self._halves
                idle = set(after.flat_out[len(Xr) + len(Yr):]) - after.wires
                self._spare = [w for w in before.flat_in if w in idle][:4]
                if len(self._spare) < 4:
                    self._halves = ()
        return bool(self._halves)

    def apply(self, Xr, Yr, inverse=False, test=None, flag=None):
        """test(S, R, t, z, k, scratch), allowed once split() holds, is emitted between the two
        blocks in place of cell 0's zero test. It must do t ^= z*[s >= 2] for the pass's final S
        value s, which sits on S with bit 0 inverted, and leave S, the n - 1 clean wires R and the
        three clean scratch wires as found. k is the test's own flag wire. It rides a layout place
        the second block never touches, so it costs no wire: a forward pass hands it to the test
        clean and returns it; an inverse pass takes it as flag, and the test must clear it."""
        C = self.C
        if test is None:
            if self.block is None:
                self.block = C.record(self._body, [Xr, Yr])
            C.apply(self.block, [Xr, Yr], inverse)
            return None
        before, after = self._halves
        width = len(Xr) + len(Yr)
        if inverse:
            pad = after.flat_out[width:]
            rest, i = C.alloc(len(pad) - 1), pad.index(self._spare[0])
            regs = [Xr, Yr, rest[:i] + [flag] + rest[i:]]
        else:
            regs = [Xr, Yr, C.alloc(len(before.flat_in) - width)]
        C.apply(after if inverse else before, regs, inverse)
        at = dict(zip(before.flat_in, flatten(regs)))
        S, R, (t,), (z,), (k, *scratch) = ([at[w] for w in ws] for ws in self._test + [self._spare])
        assert inverse is False or k == flag
        test(S, R, t, z, k, scratch)
        C.apply(before if inverse else after, regs, inverse)
        if not inverse:
            regs[2].remove(k)
        C.free(regs[2])
        return None if inverse else k


def zero_test_less(C, S, R, t, z, k, scratch, value, generic):
    """DialogMultiplier.apply's test with a compare: t ^= z*[s >= 2] and k ^= g*[s < value], for
    1 <= value < 2^n. generic(w) XORs the flag g into w.

    A CX fan writes s = r ^ y*value, y on the wire of value's top bit and r on the rest. On y = 1,
    s ^ value = r; on y = 0 it is r ^ value, whose top bit r lacks. So r's zero-prefix chain, on
    R's wires, gives both the zero test and the compare."""
    n = len(S)
    C.x(S[0])
    top = value.bit_length() - 1
    y = S[top]
    fan = [S[j] for j in range(n) if j != top and value >> j & 1]
    for w in fan:
        C.cx(y, w)
    prefix = ZeroPrefix(C, S, [j for j in reversed(range(n)) if j != top], R)
    prefix.build()
    u, w, g = scratch
    # s < 2 when r is zero above bit 0 and, unless value is 1, y = 0.
    small = prefix.above(0)
    C.cx(z, t)
    prefix.xor(small, u)
    if top == 0:
        C.ccx(z, u, t)
    else:
        C.x(y)
        C.ccx(y, u, w)
        C.ccx(z, w, t)
        C.ccx(y, u, w)
        C.x(y)
    prefix.xor(small, u)
    on, off = prefix.less(value, 0), prefix.less(value, value)

    def less():
        prefix.xor(off, w)
        prefix.xor(on ^ off, u)
        C.ccx(y, u, w)
        prefix.xor(on ^ off, u)

    less()
    generic(g)
    C.ccx(g, w, k)
    generic(g)
    less()
    prefix.clear()
    for w in reversed(fan):
        C.cx(y, w)
    C.x(S[0])


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
        self.square_offset = SquareFold(self.mod).offset if SquareFold.applies(curve.q) else 0
        # Where the fold can take 3a, the square's J is (Y^2 - 3a) mod q and the kernel adds no 3a.
        # Its gates then depend on A, so the square is emitted inline, outside the definitions.
        self._absorb = SquareFold(self.mod).absorb(3 * A[0] % curve.q) if SquareFold.applies(curve.q) else None

    def _mulsub(self, gen, Yr, Xr):
        """X = X - gen * (Y^2 - square_offset) mod q, or X - gen * (Y^2 - 3a) with _absorb, Y untouched."""
        C, mod, n = self.C, self.mod, self.n
        if SquareFold.applies(mod.q):
            self._square_mulsub(gen, Yr, Xr)
            return
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

    def _square_mulsub(self, gen, Yr, Xr):
        C, mod, n = self.C, self.mod, self.n
        fold = SquareFold(mod)
        self._modadd_blk = C.record_shape(lambda k, a, b: fold.modadd(C, a, b, k[0]), [1, n, n])
        J = C.alloc(fold.junk)
        inline = self._absorb is not None
        if inline:
            cins, kappa = self._absorb
            blk = C.record(lambda y, j: fold.compute(C, y, j, kappa, cins), [Yr, J])
        else:
            blk = C.record_shape(lambda y, j: fold.compute(C, y, j), [n, fold.junk])
        C.apply(blk, [Yr, J], inline=inline)
        C.apply(self._modadd_blk, [[gen], J[:n], Xr], inverse=True)
        C.apply(blk, [Yr, J], inverse=True, inline=inline)
        C.free(J)

    def _kernel(self, c, gen, Xr, Yr, f0, f1, below):
        """below holds gen*[x < a] and is cleared. Returns gen and a wire holding gen*[x3 <= a]."""
        C, mod, q = self.C, self.mod, self.curve.q
        a, b = self.A

        def recode(inverse=False):
            if not inverse:
                C.cx(gen, c)
                C.cx(gen, f0)
            # gen ^= NOT c AND (f0 XOR f1): no c = 0 word has f0 = 0 and f1 = 1.
            C.cx(f1, f0)
            C.x(c)
            C.ccx(c, f0, gen)
            C.x(c)
            C.cx(f1, f0)
            if inverse:
                C.cx(gen, f0)
                C.cx(gen, c)

        # Cell 0's zero test in M reads the final Y of each pass, top down. Emitted here with a
        # compare, it clears the wrap flag of the first Y -= b and makes the last one's.
        split = self.M.split(Xr, Yr)

        def generic(w):
            # gen, as recode's inverse makes it.
            C.cx(f1, f0)
            C.x(c)
            C.ccx(c, f0, w)
            C.x(c)
            C.cx(f1, f0)

        def test(value):
            return (lambda *wires: zero_test_less(C, *wires, value, generic)) if split else None

        mod.cmodadd_const(C, gen, q - a, Xr, k=below)
        wrap = mod.cmodadd_const(C, gen, q - b, Yr, keep=split)
        C.ccx(f0, f1, Xr[0])
        recode()
        C.free([gen])
        self.M.apply(Xr, Yr, inverse=True, test=test(q - b), flag=wrap)
        (gen,) = C.alloc(1)
        recode(inverse=True)
        if self._absorb is None:
            mod.cmodadd_const(C, gen, 3 * a - self.square_offset, Xr)
        self._mulsub(gen, Yr, Xr)
        recode()
        C.free([gen])
        wrap = self.M.apply(Xr, Yr, test=test(b))
        (gen,) = C.alloc(1)
        recode(inverse=True)
        C.ccx(f0, f1, Xr[0])
        below = mod.creflect(C, gen, a, Xr, keep=True)
        mod.cmodadd_const(C, gen, q - b, Yr, k=wrap)
        return gen, below

    def _eq_flags(self, c, Xr, Yr, flags, points, below=None):
        """flags[i] ^= [the point is points[i]], c-controlled except on the x = 0 point.

        below = (k, late) also does k ^= g*[x < pair_x], where g = c*~[x = 0] ^ (the other three
        flags) is the generic lane's flag. Early (late false), the flags come in 0 and g reads them
        once set. Late, they come in set and g reads them before they clear. g is 0 at x = pair_x
        either way, so this is g*[x <= pair_x] too."""
        C, n = self.C, self.n
        groups = {}
        for f, p in zip(flags, points):
            x, y = (0, 0) if p is None else p
            groups.setdefault(x, []).append((f, y))
        assert len(groups) == 3 and len(groups.get(0, [])) == 1, "0/a-pair/2A groups"
        pair_x = next(x for x, leaves in groups.items() if len(leaves) == 2)
        single_x = next(x for x in groups if x not in (0, pair_x))
        # A CX network writes x = r ^ y1*pair_x ^ y2*single_x, with y1 and y2 on the two pivot
        # wires and r, zero on the pivots, on the rest. Then x is 0, pair_x or single_x exactly when
        # r = 0 and (y1, y2) is (0, 0), (1, 0) or (0, 1), so the three equalities share r's all-zero
        # product. The pivots are the top bits of pair_x, single_x and their XOR (two places), so
        # that the top bit of x ^ pair_x is r's top bit or a pivot on every (y1, y2) lane.
        p1 = pair_x.bit_length() - 1
        p2 = single_x.bit_length() - 1
        if p2 == p1:
            p2 = (pair_x ^ single_x).bit_length() - 1
        col = lambda v: (v >> p1 & 1) | (v >> p2 & 1) << 1
        pick = None
        for ops in ((), ((0, 1),), ((1, 0),)):
            rows = [1, 2]
            for src, dst in ops:
                rows[dst] ^= rows[src]
            for i1 in (0, 1):
                # Wire i1 must read 1 on pair_x and 0 on single_x, and the other wire the reverse.
                on = lambda i, v: (rows[i] & col(v)).bit_count() & 1
                if (on(i1, pair_x), on(i1, single_x), on(1 - i1, pair_x), on(1 - i1, single_x)) == (1, 0, 0, 1):
                    pick = pick or (ops, i1)
        assert pick, "the pivot bits separate the pair and the single"
        ops, i1 = pick
        pivots = [Xr[p1], Xr[p2]]
        x1, x2 = pivots[i1], pivots[1 - i1]
        order = [j for j in reversed(range(n)) if j not in (p1, p2)]

        def conjugate(inverse=False):
            steps = [(pivots[src], pivots[dst]) for src, dst in ops]
            steps += [(x1, Xr[j]) for j in order if pair_x >> j & 1]
            steps += [(x2, Xr[j]) for j in order if single_x >> j & 1]
            for src, dst in reversed(steps) if inverse else steps:
                C.cx(src, dst)

        undo = []

        def flip(wires):
            for w in wires:
                C.x(w)

        def chain(factors, target=None):
            """AND of the factors, each a (wire, negated) literal or None for 1. With a target, XOR
            it there and clear the chain; without, return it as a factor held until undo runs."""
            factors = [f for f in factors if f is not None]
            if target is None and len(factors) < 2:
                return factors[0] if factors else None
            flips = [w for w, neg in factors if neg]
            wires = [w for w, _ in factors]
            body = wires if target is None else wires[:-1]
            made = []
            flip(flips)
            acc = body[0] if body else None
            for w in body[1:]:
                (t,) = C.alloc(1)
                C.ccx(acc, w, t)
                made.append((acc, w, t))
                acc = t
            if target is not None:
                if not wires:
                    C.x(target)
                elif acc is None:
                    C.cx(wires[-1], target)
                else:
                    C.ccx(acc, wires[-1], target)
            flip(flips)

            def clear():
                flip(flips)
                for a, b, t in reversed(made):
                    C.ccx(a, b, t)
                C.free([t for _, _, t in made])
                flip(flips)

            if target is not None:
                clear()
                return None
            undo.append(clear)
            return acc, False

        conjugate()
        prefix = ZeroPrefix(C, Xr, order)
        prefix.build()
        undo.append(prefix.clear)
        zero = prefix.zero
        y0 = Yr[0]

        def less():
            k = below[0]
            # On lane (y1, y2), x ^ pair_x = r ^ d with d = pair_x*(1 ^ y1) ^ single_x*y2.
            f00, f10, f01, f11 = (prefix.less(pair_x, (0 if y1 else pair_x) ^ (single_x if y2 else 0))
                                  for y1, y2 in ((0, 0), (1, 0), (0, 1), (1, 1)))
            (w, g) = C.alloc(2)
            terms = [(x1, f00 ^ f10), (x2, f00 ^ f01), (None, f00 ^ f10 ^ f01 ^ f11)]

            def value():
                both = None
                prefix.xor(f00, w)
                for ctl, form in terms:
                    if not form:
                        continue
                    if ctl is None:
                        (both,) = C.alloc(1)
                        C.ccx(x1, x2, both)
                        ctl = both
                    (t,) = C.alloc(1)
                    prefix.xor(form, t)
                    C.ccx(ctl, t, w)
                    prefix.xor(form, t)
                    C.free([t])
                if both is not None:
                    C.ccx(x1, x2, both)
                    C.free([both])

            def generic():
                fz = groups[0][0][0]
                C.cx(c, g)
                C.ccx(c, fz, g)
                for f in flags:
                    if f != fz:
                        C.cx(f, g)

            value()
            generic()
            C.ccx(g, w, k)
            generic()
            value()
            C.free([w, g])

        if below and below[1]:
            less()
        # The zero leaf is uncontrolled; c joins only on the nonzero branches.
        chain([zero, (x1, True), (x2, True)], groups[0][0][0])
        common = chain([zero, (c, False)])
        leaf = chain([common, (x1, False), (x2, True)])
        for f, y in groups[pair_x]:
            chain([leaf, (y0, not y & 1)], f)
        f, y = groups[single_x][0]
        chain([common, (x1, True), (x2, False), (y0, not y & 1)], f)
        if below and not below[1]:
            less()
        for step in reversed(undo):
            step()
        conjugate(inverse=True)

    def build(self, hierarchical=False):
        curve, A = self.curve, self.A
        twoA = curve.add(A, A)
        if A[0] == 0 or twoA is None or twoA[0] == 0:
            raise ValueError("kernel zero flag requires nonzero x(A) and x(2A)")
        if pow(curve.b, (curve.q - 1) // 2, curve.q) != curve.q - 1:
            raise ValueError("kernel zero flag refuses curves with a point at x = 0")
        self.C = C = Circuit(hierarchical)
        n = self.n
        (c,) = C.alloc(1)
        Xr, Yr = C.alloc(n), C.alloc(n)
        declared = Xr + Yr
        self.M = DialogMultiplier(C, self.mod)
        self._modadd_blk = None

        flags = C.alloc(4)
        # The flag blocks' zero chain also compares x with a: the first block makes the wrap flag
        # of the kernel's X -= a, and the second clears the reflect's.
        (below,) = C.alloc(1)
        self._eq_flags(c, Xr, Yr, flags, [e for e, _ in self.exceptions], below=(below, False))
        f0, f1, f2, gen = flags
        for f, (e_in, e_out) in zip(flags[1:], self.exceptions[1:]):
            load_const(C, Xr + Yr, encode(e_in, n) ^ encode(e_out, n), f)
        # Move the uncontrolled x = 0 marker to f2 when c = 0.
        C.x(c)
        C.ccx(c, f0, f2)
        C.ccx(c, f2, f0)
        C.x(c)
        e_in, e_out = self.exceptions[0]
        load_const(C, Xr + Yr, encode(e_in, n) ^ encode(e_out, n), f0)
        # Only the c-controlled part of f2 contributes to the generic flag.
        for w in (c, f0, f1):
            C.cx(w, gen)
        C.ccx(c, f2, gen)
        # Two wires name the one set flag among three; the wire left over is zero.
        C.cx(f2, f0)
        C.cx(f2, f1)
        C.ccx(f0, f1, f2)
        C.free([f2])

        gen, below = self._kernel(c, gen, Xr, Yr, f0, f1, below)

        (f2,) = C.alloc(1)
        C.ccx(f0, f1, f2)
        C.cx(f2, f1)
        C.cx(f2, f0)
        for w in (c, f0, f1):
            C.cx(w, gen)
        C.ccx(c, f2, gen)
        flags = [f0, f1, f2, gen]
        self._eq_flags(c, Xr, Yr, flags, [e for _, e in self.exceptions], below=(below, True))
        C.free(flags + [below])

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
