"""Stream a circuit as v2-ir `ArithmeticCircuit` JSON (serde field names from
contracts/ir/src/arithmetic.rs)."""
from .builder import CCX, CX, SWAP, X

CHUNK = 1 << 16


def _op(op, a, b, c):
    if op == CCX:
        return f'{{"op":"CCX","control1":{a},"control2":{b},"target":{c}}}'
    if op == CX:
        return f'{{"op":"CX","control":{a},"target":{b}}}'
    if op == X:
        return f'{{"op":"X","target":{a}}}'
    assert op == SWAP
    return f'{{"op":"SWAP","wire1":{a},"wire2":{b}}}'


def definitions(C):
    """Blocks the top level reaches, each after every block it calls."""
    order = []
    seen = set()

    def visit(calls):
        for call in calls:
            if call.block not in seen:
                seen.add(call.block)
                visit(call.block.calls)
                order.append(call.block)

    visit(C.calls)
    return order


def _instr_chunks(ops, calls, ids, local=None):
    def prims(lo, hi):
        for start in range(4 * lo, 4 * hi, 4 * CHUNK):
            it = iter(ops[start:min(start + 4 * CHUNK, 4 * hi)])
            if local is None:
                yield ",".join(_op(*t) for t in zip(it, it, it, it))
            else:
                yield ",".join(_op(op, local[a], local[b], local[c]) for op, a, b, c in zip(it, it, it, it))

    pos = 0
    for call in calls:
        yield from prims(pos, call.pos)
        args = call.wires if local is None else [local[w] for w in call.wires]
        yield f'{{"call":{ids[call.block]},"wires":[{",".join(map(str, args))}],"inverse":{"true" if call.inverse else "false"}}}'
        pos = call.pos
    yield from prims(pos, len(ops) // 4)


def _write_list(fh, pieces):
    total = fh.write("[")
    first = True
    for piece in pieces:
        if piece:
            total += fh.write(piece if first else "," + piece)
            first = False
    return total + fh.write("]")


def write_circuit(fh, C, layout):
    """Writes the JSON text to a text-mode handle and returns the counts it
    wrote: json_bytes, plus definitions, call_sites and primitive_ops for a
    hierarchical circuit, which also gets `definitions`, one per reachable block."""
    live = set(layout["control"]) | set(layout["x"]) | set(layout["y"])
    regs = [
        ("control", "CONTROL", layout["control"]),
        ("accumulator", "ACCUMULATOR", layout["x"] + layout["y"]),
        ("scratch", "SCRATCH", [w for w in range(C.num_wires) if w not in live]),
    ]
    reg_json = ",".join(f'{{"name":"{n}","role":"{r}","wires":[{",".join(map(str, ws))}]}}' for n, r, ws in regs)
    total = fh.write(f'{{"num_wires":{C.num_wires},"registers":[{reg_json}],')
    calls, ids, counts = [], {}, {}
    if C.calls is not None:
        calls = C.calls
        total += fh.write('"definitions":[')
        for i, block in enumerate(definitions(C)):
            order = sorted(block.wires)
            local = {w: k for k, w in enumerate(order)}
            total += fh.write(("," if i else "") + f'{{"arity":{len(order)},"operations":')
            total += _write_list(fh, _instr_chunks(block.ops, block.calls, ids, local))
            total += fh.write("}")
            ids[block] = i
        total += fh.write("],")
        counts = {"definitions": len(ids), "call_sites": len(calls) + sum(len(b.calls) for b in ids),
                  "primitive_ops": (len(C.prims) + sum(len(b.ops) for b in ids)) // 4}
    total += fh.write('"operations":')
    total += _write_list(fh, _instr_chunks(C.prims, calls, ids))
    total += fh.write("}")
    return {**counts, "json_bytes": total}
