import json, sys
m = json.load(open(sys.argv[1])); c = m.get("circuit", m)
done = []
def body(ops):
    n = 0
    for r in ops:
        n += done[r["call"]] if "call" in r else (r["op"] in ("CCX", "CCZ"))
    return n
for d in c.get("definitions", []): done.append(body(d["operations"]))
n = body(c["operations"]); print(sys.argv[1].split("/")[-1], "wires", c["num_wires"], "CCX", n, "T", 7 * n)
