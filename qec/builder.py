"""Reversible circuit builder over the v2 arithmetic IR's permutation subset.

An op is four int32 values (opcode, a, b, c). Unused operand slots repeat the
last used wire so wire remapping can treat every slot uniformly.
"""
from array import array

X, CX, CCX, SWAP = 0, 1, 2, 3
OPNAMES = {X: "X", CX: "CX", CCX: "CCX", SWAP: "SWAP"}


class Block:
    """A recorded op range that can be replayed forward or reversed onto any
    register layout of the same shape. Wires outside the layout are zero at the
    block's start and end and are mapped onto fresh ancillas."""

    def __init__(self, ops, flat_in, flat_out):
        self.ops = ops
        self.flat_in = flat_in
        self.flat_out = flat_out
        self.wires = set(ops[1::4]) | set(ops[2::4]) | set(ops[3::4])


def flatten(regs):
    return [w for r in regs for w in r]


def unflatten_into(regs, flat):
    i = 0
    for r in regs:
        r[:] = flat[i:i + len(r)]
        i += len(r)


class Circuit:
    def __init__(self):
        self.ops = array("i")
        self.num_wires = 0
        self.free_stack = []

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
        self.ops.extend((X, a, a, a))

    def cx(self, a, b):
        self.ops.extend((CX, a, b, b))

    def ccx(self, a, b, c):
        if a == b:
            self.ops.extend((CX, a, c, c))
        else:
            self.ops.extend((CCX, a, b, c))

    def swap(self, a, b):
        self.ops.extend((SWAP, a, b, b))

    def __len__(self):
        return len(self.ops) // 4

    # blocks

    def record(self, fn, regs):
        """Run fn(*regs) without emitting it. Returns a Block and leaves the
        circuit, the allocator, and regs exactly as they were."""
        flat_in = flatten(regs)
        start = len(self.ops)
        saved_stack, saved_hw = list(self.free_stack), self.num_wires
        fn(*regs)
        flat_out = flatten(regs)
        ops = self.ops[start:]
        del self.ops[start:]
        self.free_stack, self.num_wires = saved_stack, saved_hw
        unflatten_into(regs, flat_in)
        return Block(ops, flat_in, flat_out)

    def apply(self, block, regs, inverse=False):
        src_live, dst_live = (block.flat_out, block.flat_in) if inverse else (block.flat_in, block.flat_out)
        cur = flatten(regs)
        assert len(cur) == len(src_live)
        if len(block.ops) == 0:
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
        self.ops.extend(out)
        new_live = [table[w] for w in dst_live]
        live = set(new_live)
        self.free([w for w in list(cur) + fresh if w not in live])
        unflatten_into(regs, new_live)

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

