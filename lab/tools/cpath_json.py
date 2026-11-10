"""Critical path of a JSON circuit under the board tick model, attributed to definitions.

    python3 -I cpath.py <circuit.json> <labels.json> [top_k] [def_id_for_op_detail ...]

Each executed op starts at S: X, CX, SWAP at max ready of its wires; CCX at max(r_c1, r_c2 + 2, r_t + 3)
(its target is then free at S + 7 and its controls at S + 8). The wire achieving that max is the op's
binding wire, and its predecessor is the last op on that wire. Walking predecessors back from the op
that finishes last gives one critical path; each op on it is charged S(next) - S(op) ticks.
"""
import json, sys
from array import array
from collections import Counter

doc = json.load(open(sys.argv[1])); c = doc.get("circuit", doc)
labels = {d["id"]: d["label"] for d in json.load(open(sys.argv[2]))}
K = int(sys.argv[3]) if len(sys.argv) > 3 else 25
detail = [int(x) for x in sys.argv[4:]]

def comp(ops):
    out = []
    for r in ops:
        if "call" in r:
            out.append((7, r["call"], tuple(r["wires"]), bool(r.get("inverse", False))))
        else:
            o = r["op"]
            if o == "X": out.append((0, r["target"], 0, 0))
            elif o == "CX": out.append((1, r["control"], r["target"], 0))
            elif o == "SWAP": out.append((2, r["wire1"], r["wire2"], 0))
            else: out.append((3, r["control1"], r["control2"], r["target"]))
    return out

bodies = [comp(d["operations"]) for d in c["definitions"]]
top = comp(c["operations"])
nw = c["num_wires"]
del doc, c
ready = [0] * nw
last = [-1] * nw
S = array("i"); pred = array("i"); bid = array("i"); oid = array("i"); stk = array("i")
stacks = {(): 0}

def run(body_id, ops, wmap, inverse, sid):
    n = len(ops)
    for idx in (range(n - 1, -1, -1) if inverse else range(n)):
        k, a, b, t = ops[idx]
        if k == 7:
            sub = [wmap[w] for w in b] if wmap is not None else list(b)
            s2 = stacks.setdefault((sid, a), len(stacks))
            run(a, bodies[a], sub, inverse ^ t, s2)
            continue
        if wmap is not None:
            a = wmap[a]
            if k: b = wmap[b]
            if k == 3: t = wmap[t]
        i = len(S)
        if k == 0:
            s = ready[a]; pr = last[a]
            ready[a] = s + 1; last[a] = i
        elif k == 3:
            ra, rb, rt = ready[a], ready[b] + 2, ready[t] + 3
            if ra >= rb and ra >= rt: s, pr = ra, last[a]
            elif rb >= rt: s, pr = rb, last[b]
            else: s, pr = rt, last[t]
            ready[t] = s + 7; ready[a] = ready[b] = s + 8
            last[a] = last[b] = last[t] = i
        else:
            ra, rb = ready[a], ready[b]
            if ra >= rb: s, pr = ra, last[a]
            else: s, pr = rb, last[b]
            ready[a] = ready[b] = s + (1 if k == 1 else 3); last[a] = last[b] = i
        S.append(s); pred.append(pr); bid.append(-1 if body_id is None else body_id); oid.append(idx); stk.append(sid)

run(None, top, None, False, 0)
total = max(ready)
print("ops executed", len(S), "ticks", total)
wmax = max(range(nw), key=lambda w: ready[w])
i = last[wmax]
path = []
while i >= 0:
    path.append(i)
    i = pred[i]
path.reverse()
charge = Counter(); charge_op = Counter(); charge_stack = Counter()
inv = {v: k for k, v in stacks.items()}
for j, i in enumerate(path):
    nxt = S[path[j + 1]] if j + 1 < len(path) else total
    d = nxt - S[i]
    charge[bid[i]] += d
    charge_op[(bid[i], oid[i])] += d
    charge_stack[stk[i]] += d
print("path ops", len(path), "path start", S[path[0]])
print("by definition:")
for did, v in charge.most_common(K):
    print(f"  {v:9d} {100*v/total:6.2f}%  def {did} {labels.get(did, 'top')[:80]}")
def stackname(sid):
    names = []
    while sid:
        p, a = inv[sid]
        names.append(str(a))
        sid = p
    return "/".join(reversed(names))
print("by call stack (def ids, outermost first):")
for sid, v in charge_stack.most_common(K):
    print(f"  {v:9d} {100*v/total:6.2f}%  {stackname(sid)}")
print("hot ops:")
for (did, idx), v in charge_op.most_common(K):
    print(f"  {v:9d} def {did} op {idx} {bodies[did][idx] if did >= 0 else top[idx]}")
import os
phases = json.load(open(os.environ["PHASES"])) if os.environ.get("PHASES") else {}
execs = Counter()
for j, i in enumerate(path):
    pass
for d in detail:
    names = phases.get(str(d))
    agg = Counter()
    for (did, idx), v in charge_op.items():
        if did == d:
            kind = "X CX SWAP CCX".split()[(bodies[d] if d >= 0 else top)[idx][0]]
            nm = names[idx] if names else str(idx)
            if os.environ.get("COARSE"):
                nm = "/".join(x.split("#")[0] for x in nm.split("/")[:int(os.environ["COARSE"])])
            agg[nm + ":" + kind] += v
    tot = sum(agg.values())
    print(f"def {d} {labels.get(d, '')[:70]}: {tot} ticks on the path, by phase:")
    for k, v in sorted(agg.items(), key=lambda kv: -kv[1])[:K]:
        print(f"  {v:9d} {100*v/max(tot,1):6.2f}%  {k}")
if os.environ.get("AROUND"):
    tgt = int(os.environ["AROUND"]); shown = 0
    j = 0
    while j < len(path) and shown < int(os.environ.get("SHOW", "4")):
        if bid[path[j]] == tgt and (j == 0 or bid[path[j - 1]] != tgt):
            k = j
            while k < len(path) and bid[path[k]] == tgt:
                k += 1
            def desc(x):
                i = path[x]; b = bid[i]
                return f"[{stackname(stk[i])}] def {b} op {oid[i]} {bodies[b][oid[i]] if b >= 0 else top[oid[i]]} S={S[i]}"
            print("run of", k - j, "ops,", S[path[k]] - S[path[j]] if k < len(path) else "-", "ticks")
            for x in range(max(0, j - 3), j): print("   before", desc(x))
            print("   first ", desc(j))
            print("   last  ", desc(k - 1))
            for x in range(k, min(len(path), k + 3)): print("   after ", desc(x))
            shown += 1
            j = k
        else:
            j += 1
