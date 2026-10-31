"""Expand a hierarchical `ArithmeticCircuit` JSON into the bytes of its flat
`operations` array, by macro expansion. A forward call substitutes its wires
into the callee's body. An inverse call does the same over the body reversed,
with each nested call's inverse flag flipped. Shares no code with the builder,
so it can check the builder's hierarchy against the flat stream."""
import json
from operator import itemgetter
from typing import Callable, NamedTuple

CHUNK = 1 << 16
FIELDS = {"X": ("target",), "CX": ("control", "target"), "CCX": ("control1", "control2", "target"), "SWAP": ("wire1", "wire2")}
_templates = {name: "{" + ",".join([f'"op":"{name}"'] + [f'"{k}":%d' for k in keys]) + "}" for name, keys in FIELDS.items()}


class Op(NamedTuple):
    template: str
    wires: tuple
    pick: Callable


class Call(NamedTuple):
    callee: int
    wires: list
    inverse: bool


def _wires(ws):
    return all(type(w) is int and w >= 0 for w in ws)


def _object(pairs):
    d = dict(pairs)
    keys = FIELDS.get(d["op"]) if type(d.get("op")) is str else None
    if keys is not None and len(d) == len(keys) + 1 and all(k in d for k in keys):
        wires = tuple(d[k] for k in keys)
        if _wires(wires):
            return Op(_templates[d["op"]], wires, itemgetter(*wires))
    if d.keys() == {"call", "wires", "inverse"}:
        callee, wires, inverse = d["call"], d["wires"], d["inverse"]
        if type(callee) is int and type(inverse) is bool and type(wires) is list and _wires(wires):
            return Call(callee, wires, inverse)
    return d


def load(fh):
    """Parse and validate the call graph: every instruction is an X, CX, CCX,
    SWAP or call whose wires lie in its body's range, and every call names a
    lower-numbered definition with a distinct wire map of its arity."""
    doc = json.load(fh, object_pairs_hook=_object)
    defs = doc.setdefault("definitions", [])
    bodies = [(f"definition {i}", d["arity"], d["operations"]) for i, d in enumerate(defs)]
    for limit, (where, width, body) in enumerate(bodies + [("top level", doc["num_wires"], doc["operations"])]):
        for j, ins in enumerate(body):
            if type(ins) is Call:
                if not 0 <= ins.callee < limit:
                    raise ValueError(f"{where} instruction {j}: call {ins.callee} is not an earlier definition")
                if len(ins.wires) != defs[ins.callee]["arity"]:
                    raise ValueError(f"{where} instruction {j}: {len(ins.wires)} wires for arity {defs[ins.callee]['arity']}")
                if len(set(ins.wires)) != len(ins.wires):
                    raise ValueError(f"{where} instruction {j}: call wires are not distinct: {ins.wires}")
            elif type(ins) is not Op:
                raise ValueError(f"{where} instruction {j}: not an X, CX, CCX, SWAP or call: {json.dumps(ins)[:200]}")
            if max(ins.wires, default=-1) >= width:
                raise ValueError(f"{where} instruction {j}: wire {max(ins.wires)} out of range for {width} wires")
    return doc


def operations_json(doc):
    """Yields the flat `operations` array text, `[` to `]`, in chunks of CHUNK ops."""
    defs = [d["operations"] for d in doc["definitions"]]
    out = []

    def expand(body, wires, inverse):
        for ins in reversed(body) if inverse else body:
            if type(ins) is Op:
                out.append(ins.template % ins.pick(wires))
                if len(out) >= CHUNK:
                    yield ",".join(out)
                    out.clear()
            else:
                yield from expand(defs[ins.callee], [wires[w] for w in ins.wires], ins.inverse != inverse)

    yield "["
    sep = ""
    for chunk in expand(doc["operations"], range(doc["num_wires"]), False):
        yield sep + chunk
        sep = ","
    if out:
        yield sep + ",".join(out)
    yield "]"
