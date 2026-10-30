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


def write_circuit(fh, C, layout):
    """Writes the JSON text to a text-mode handle and returns bytes written."""
    live = set(layout["control"]) | set(layout["x"]) | set(layout["y"])
    regs = [
        ("control", "CONTROL", layout["control"]),
        ("accumulator", "ACCUMULATOR", layout["x"] + layout["y"]),
        ("scratch", "SCRATCH", [w for w in range(C.num_wires) if w not in live]),
    ]
    reg_json = ",".join(f'{{"name":"{n}","role":"{r}","wires":[{",".join(map(str, ws))}]}}' for n, r, ws in regs)
    head = f'{{"num_wires":{C.num_wires},"registers":[{reg_json}],"operations":['
    total = fh.write(head)
    ops = C.ops
    first = True
    for start in range(0, len(ops), 4 * CHUNK):
        it = iter(ops[start:start + 4 * CHUNK])
        body = ",".join(_op(*t) for t in zip(it, it, it, it))
        if body:
            total += fh.write(body if first else "," + body)
            first = False
    total += fh.write("]}")
    return total
