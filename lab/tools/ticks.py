"""Independent ASAP tick estimate: CCX lowered to the textbook 15-gate Clifford+T sequence, every gate 1 tick."""
import json, sys
m = json.load(open(sys.argv[1])); c = m.get("circuit", m)
OPS = {"X": 0, "CX": 1, "CCX": 2, "SWAP": 3, "CZ": 4, "CCZ": 5, "Z": 6}
def comp(ops):
    out = []
    for r in ops:
        if "call" in r: out.append(("call", r["call"], r["wires"], r.get("inverse", False)))
        else:
            o = r["op"]
            if o in ("X", "Z"): out.append((o, r["target"]))
            elif o in ("CX", "CZ"): out.append((o, r["control"], r["target"]))
            elif o in ("CCX", "CCZ"): out.append((o, r["control1"], r["control2"], r["target"]))
            elif o == "SWAP": out.append((o, r["wire1"], r["wire2"]))
            else: raise SystemExit("unhandled op " + o)
    return out
defs = [comp(d["operations"]) for d in c.get("definitions", [])]
top = comp(c["operations"])
r = [0] * c["num_wires"]
ccx = 0
def run(body, wmap, inv):
    global ccx
    seq = reversed(body) if inv else body
    for g in seq:
        k = g[0]
        if k == "call":
            sub = [wmap[w] for w in g[2]] if wmap is not None else g[2]
            run(defs[g[1]], sub, inv ^ g[3]); continue
        w = [wmap[x] for x in g[1:]] if wmap is not None else list(g[1:])
        if k in ("X", "Z"): r[w[0]] += 1
        elif k in ("CX", "CZ"):
            t = max(r[w[0]], r[w[1]]) + 1; r[w[0]] = r[w[1]] = t
        elif k == "SWAP":
            t = max(r[w[0]], r[w[1]]) + 3; r[w[0]] = r[w[1]] = t
        else:
            ccx += 1
            a, b, t = w; ra, rb, rt = r[a], r[b], r[t]
            rt += 1                      # H t
            x = max(rb, rt) + 1; rb = rt = x   # CX b->t
            rt += 1                      # Tdg t
            x = max(ra, rt) + 1; ra = rt = x   # CX a->t
            rt += 1                      # T t
            x = max(rb, rt) + 1; rb = rt = x   # CX b->t
            rt += 1                      # Tdg t
            x = max(ra, rt) + 1; ra = rt = x   # CX a->t
            rb += 1; rt += 1             # T b ; T t
            rt += 1                      # H t
            x = max(ra, rb) + 1; ra = rb = x   # CX a->b
            ra += 1; rb += 1             # T a ; Tdg b
            x = max(ra, rb) + 1; ra = rb = x   # CX a->b
            r[a], r[b], r[t] = ra, rb, rt
sys.setrecursionlimit(100)
run(top, None, False)
print("CCX", ccx, "ticks", max(r), "wires", c["num_wires"], "qubit_time", max(r) * c["num_wires"])
