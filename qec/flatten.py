"""Expand a hierarchical `ArithmeticCircuit` JSON into the bytes of its flat
`operations` array, by macro expansion. A forward call substitutes its wires
into the callee's body. An inverse call does the same over the body reversed,
with each nested call's inverse flag flipped. Shares no code with the builder,
so it can check the builder's hierarchy against the flat stream."""
import json
from operator import itemgetter
from typing import Callable, NamedTuple

CHUNK = 1 << 16
_templates = {}


class Op(NamedTuple):
    template: str
    pick: Callable


class Call(NamedTuple):
    callee: int
    wires: list
    inverse: bool


def _object(pairs):
    keys = tuple(k for k, _ in pairs)
    if keys[:1] == ("op",):
        name, wires = pairs[0][1], [v for _, v in pairs[1:]]
        assert wires and all(type(w) is int and w >= 0 for w in wires), pairs
        key = (name, keys)
        if key not in _templates:
            _templates[key] = "{" + ",".join([f'"op":{json.dumps(name)}'] + [f"{json.dumps(k)}:%d" for k in keys[1:]]) + "}"
        return Op(_templates[key], itemgetter(*wires))
    if keys == ("call", "wires", "inverse"):
        callee, wires, inverse = (v for _, v in pairs)
        assert type(callee) is int and type(inverse) is bool, pairs
        assert all(type(w) is int and w >= 0 for w in wires) and len(set(wires)) == len(wires), pairs
        return Call(callee, wires, inverse)
    return dict(pairs)


def load(fh):
    """Parse and validate the call graph: every call names a lower-numbered
    definition with the full wire map of its arity."""
    doc = json.load(fh, object_pairs_hook=_object)
    defs = doc.setdefault("definitions", [])
    for limit, body in [*enumerate(d["operations"] for d in defs), (len(defs), doc["operations"])]:
        for call in body:
            if type(call) is Call:
                assert 0 <= call.callee < limit and len(call.wires) == defs[call.callee]["arity"], (limit, call.callee)
    return doc


def operations_json(doc):
    """Yields the flat `operations` array text, `[` to `]`, in chunks."""
    defs = [d["operations"] for d in doc["definitions"]]
    out = []

    def expand(body, wires, inverse):
        for ins in reversed(body) if inverse else body:
            if type(ins) is Op:
                out.append(ins.template % ins.pick(wires))
            else:
                yield from expand(defs[ins.callee], [wires[w] for w in ins.wires], ins.inverse != inverse)
                if len(out) >= CHUNK:
                    yield ",".join(out)
                    out.clear()

    yield "["
    sep = ""
    for chunk in expand(doc["operations"], range(doc["num_wires"]), False):
        yield sep + chunk
        sep = ","
    if out:
        yield sep + ",".join(out)
    yield "]"
