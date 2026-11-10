"""A top-down zero-prefix chain and the compares it gives for CX gates only.

Z_i = [r is zero on its top i places], i = 0 .. m, over r's places listed top first. Z_i ^ Z_(i+1)
says r's top set bit is place i of that list. Take x = r ^ e and a constant v, with d = e ^ v either
0 or with its top bit on a place r lacks. Then the top bit of x ^ v = r ^ d is r's top set bit when
that lies above d's top bit, and d's top bit otherwise, and x < v exactly when v has that bit. So
[x < v] is an XOR of chain wires: top_in(v, top of d), and above(top of d) when v has d's top bit."""


class ZeroPrefix:
    def __init__(self, C, wires, places, scratch=()):
        """wires[j] holds r's bit j; places lists the bits r lives on, top first. The chain's AND
        wires come from scratch, clean wires that clear() leaves clean, then from the allocator."""
        self.C, self.wires, self.places = C, wires, list(places)
        self.index = {j: i for i, j in enumerate(self.places)}
        self.scratch = list(scratch)
        self.ands = []

    def _flip(self):
        for j in self.places:
            self.C.x(self.wires[j])

    def build(self):
        C = self.C
        self._flip()
        acc = None
        for i, j in enumerate(self.places):
            if i == 0:
                acc = self.wires[j]
                continue
            (t,) = [self.scratch[i - 1]] if i - 1 < len(self.scratch) else C.alloc(1)
            C.ccx(acc, self.wires[j], t)
            self.ands.append((acc, self.wires[j], t))
            acc = t
        self._flip()

    def clear(self):
        C = self.C
        self._flip()
        for a, b, t in reversed(self.ands):
            C.ccx(a, b, t)
        self._flip()
        C.free([t for i, (_, _, t) in enumerate(self.ands) if i >= len(self.scratch)])
        self.ands = []

    @property
    def zero(self):
        """[r = 0] as a (wire, negated) literal, or None when r has no places."""
        if not self.places:
            return None
        if len(self.places) == 1:
            return self.wires[self.places[0]], True
        return self.ands[-1][2], False

    def xor(self, form, t):
        """t ^= the XOR of Z_i over i in form."""
        C = self.C
        for i in sorted(form):
            if i == 0:
                C.x(t)
            elif i == 1:
                C.cx(self.wires[self.places[0]], t)
                C.x(t)
            else:
                C.cx(self.ands[i - 2][2], t)

    def above(self, place):
        """{i} with Z_i = [r is zero above this place]."""
        return {sum(1 for j in self.places if j > place)}

    def top_in(self, mask, lo=-1):
        """The form of [r's top set bit is above lo and is a set bit of mask]."""
        form = set()
        for j in self.places:
            if j > lo and mask >> j & 1:
                form ^= {self.index[j], self.index[j] + 1}
        return form

    def less(self, v, d):
        """The form of [r ^ d ^ v < v], for d = 0 or with its top bit on a place r lacks."""
        if not d:
            return self.top_in(v)
        top = d.bit_length() - 1
        assert top not in self.index
        return self.top_in(v, top) ^ (self.above(top) if v >> top & 1 else set())
