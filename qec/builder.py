"""Reversible circuit builder over the v2 arithmetic IR's permutation subset.

An op is four int32 values (opcode, a, b, c). Unused operand slots repeat the
last used wire so wire remapping can treat every slot uniformly.

In hierarchical mode a body is its primitive op array plus a list of Calls,
each run just before primitive op `pos`. `wires` holds the body's own wire for
each of sorted(block.wires), so applying a block never remaps the calls inside
it.
"""
from array import array
from typing import NamedTuple

X, CX, CCX, SWAP = 0, 1, 2, 3
OPNAMES = {X: "X", CX: "CX", CCX: "CCX", SWAP: "SWAP"}


class Call(NamedTuple):
    pos: int
    block: "Block"
    wires: list
    inverse: bool


class Block:
    """A recorded op range that can be replayed forward or reversed onto any
    register layout of the same shape. Wires outside the layout are zero at the
    block's start and end and are mapped onto fresh ancillas."""

    def __init__(self, ops, flat_in, flat_out, calls=()):
        self.ops = ops
        self.flat_in = flat_in
        self.flat_out = flat_out
        self.calls = calls
        self.wires = set(ops[1::4]) | set(ops[2::4]) | set(ops[3::4])
        for call in calls:
            self.wires.update(call.wires)


def flatten(regs):
    return [w for r in regs for w in r]


def unflatten_into(regs, flat):
    i = 0
    for r in regs:
        r[:] = flat[i:i + len(r)]
        i += len(r)


class Circuit:
    def __init__(self, hierarchical=False):
        self.prims = array("i")
        self.calls = [] if hierarchical else None
        self.num_wires = 0
        self.free_stack = []
        self._depth = 0
        self._head = None

    # allocation

    def alloc(self, k=1):
        out = []
        for _ in range(k):
            if self.free_stack:
                out.append(self.free_stack.pop())
            else:
                out.append(self.num_wires)
                self.num_wires += 1
        return out

    def free(self, wires):
        self.free_stack.extend(reversed(wires))

    # gates

    def x(self, a):
        self.prims.extend((X, a, a, a))

    def cx(self, a, b):
        self.prims.extend((CX, a, b, b))

    def ccx(self, a, b, c):
        if a == b:
            self.prims.extend((CX, a, c, c))
        else:
            self.prims.extend((CCX, a, b, c))

    def swap(self, a, b):
        self.prims.extend((SWAP, a, b, b))

    @property
    def ops(self):
        if self.calls is not None:
            raise TypeError("a hierarchical Circuit holds only its top-level primitives; expand it with qec.flatten")
        return self.prims

    def __len__(self):
        return len(self.ops) // 4

    # blocks

    def _take(self, start, first_call):
        """Remove and return the ops and calls emitted since (start, first_call)."""
        ops = self.prims[start:]
        del self.prims[start:]
        calls = ()
        if self.calls is not None:
            calls = [call._replace(pos=call.pos - start // 4) for call in self.calls[first_call:]]
            del self.calls[first_call:]
        return ops, calls

    def _record(self, fn, regs):
        flat_in = flatten(regs)
        mark = len(self.prims), len(self.calls) if self.calls is not None else 0
        saved_stack, saved_hw = list(self.free_stack), self.num_wires
        outer, self._head = self._head, (mark, None)
        self._depth += 1
        fn(*regs)
        self._depth -= 1
        flat_out = flatten(regs)
        (_, head), self._head = self._head, outer
        tail = self._take(*mark)
        self.free_stack, self.num_wires = saved_stack, saved_hw
        unflatten_into(regs, flat_in)
        return flat_in, flat_out, head, tail

    def record(self, fn, regs):
        """Run fn(*regs) without emitting it. Returns a Block and leaves the
        circuit, the allocator, and regs exactly as they were."""
        flat_in, flat_out, head, tail = self._record(fn, regs)
        assert head is None, "cut() needs record_cut"
        return Block(tail[0], flat_in, flat_out, tail[1])

    def cut(self):
        """Inside record_cut's fn, outside any nested record: end the first block here."""
        mark, head = self._head
        assert self._depth == 1 and head is None
        self._head = mark, self._take(*mark)

    def record_cut(self, fn, regs):
        """record for an fn that calls cut() once: returns the blocks before and after the cut, or
        None when fn never cut. The first maps regs plus clean wires, one per other wire it
        touches, onto themselves, so every wire live at the cut is in its layout, in flat_in's
        order; the second takes that layout and ends on regs plus clean layout wires."""
        flat_in, flat_out, head, tail = self._record(fn, regs)
        if head is None:
            return None
        seen = Block(head[0], [], [], head[1]).wires
        live = flat_in + sorted(seen - set(flat_in))
        idle = [w for w in live if w not in set(flat_out)]
        pad = len(live) - len(flat_out)
        assert 0 <= pad <= len(idle)
        return Block(head[0], live, live, head[1]), Block(tail[0], live, flat_out + idle[:pad], tail[1])

    def apply(self, block, regs, inverse=False, inline=False):
        """Replay block onto regs. inline emits its ops here even in hierarchical mode, so a block
        whose ops depend on the addend stays out of the definitions; it must hold no calls."""
        src_live, dst_live = (block.flat_out, block.flat_in) if inverse else (block.flat_in, block.flat_out)
        cur = flatten(regs)
        assert len(cur) == len(src_live)
        if not block.wires:
            unflatten_into(regs, [dict(zip(src_live, cur))[w] for w in dst_live])
            return
        size = max(max(block.wires), max(src_live, default=0), max(dst_live, default=0)) + 1
        table = [-1] * size
        for w, v in zip(src_live, cur):
            table[w] = v
        scratch = sorted(w for w in block.wires if table[w] < 0)
        assert set(dst_live) <= set(src_live) | block.wires
        fresh = self.alloc(len(scratch))
        for w, v in zip(scratch, fresh):
            table[w] = v
        if self.calls is None or inline:
            assert not block.calls
            self._replay(block, table, inverse)
        else:
            wires = [table[w] for w in sorted(block.wires)]
            if len(set(wires)) != len(wires):
                raise ValueError(f"call wires are not distinct: {wires}")
            self.calls.append(Call(len(self.prims) // 4, block, wires, inverse))
        new_live = [table[w] for w in dst_live]
        live = set(new_live)
        self.free([w for w in list(cur) + fresh if w not in live])
        unflatten_into(regs, new_live)

    def _replay(self, block, table, inverse):
        src = block.ops
        out = array("i", bytes(len(src) * src.itemsize))
        if inverse:
            rev = src[::-1]
            out[0::4] = rev[3::4]
            slots = (rev[2::4], rev[1::4], rev[0::4])
        else:
            out[0::4] = src[0::4]
            slots = (src[1::4], src[2::4], src[3::4])
        get = table.__getitem__
        for k, slot in enumerate(slots, 1):
            out[k::4] = array("i", map(get, slot))
        self.prims.extend(out)

    def record_shape(self, fn, sizes):
        """Record fn on fresh placeholder registers of the given sizes."""
        saved_stack, saved_hw = list(self.free_stack), self.num_wires
        regs = [self.alloc(k) for k in sizes]
        block = self.record(fn, regs)
        self.free_stack, self.num_wires = saved_stack, saved_hw
        return block

    def run(self, fn, regs):
        block = self.record(fn, regs)
        self.apply(block, regs)
        return block

