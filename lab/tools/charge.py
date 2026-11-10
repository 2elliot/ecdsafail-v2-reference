import json, math, sys
c = json.load(open(sys.argv[1])); c = c.get("circuit", c)
done = []
def body(ops):
    t = {"ops": 0, "calls": 0, "call_wires": 0}
    for r in ops:
        if "call" in r:
            k = done[r["call"]]
            for key in t: t[key] += k[key]
            t["calls"] += 1; t["call_wires"] += len(r["wires"])
        else:
            t["ops"] += 1
    return t
for d in c.get("definitions", []): done.append(body(d["operations"]))
t = body(c["operations"])
print("expanded ops", t["ops"], "calls", t["calls"], "call wires", t["call_wires"],
      "charge", t["ops"] + math.ceil((1_329_430 * t["calls"] + 17_360 * t["call_wires"] + 918_090) / 20_710))
