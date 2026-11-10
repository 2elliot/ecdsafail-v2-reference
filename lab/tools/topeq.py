"""Compare two circuits' whole top-level maps on random full states (all wires, any values)."""
import json, random, sys
sys.path.insert(0, sys.argv[3])
import qlocal
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
ca, cb = a.get("circuit", a), b.get("circuit", b)
assert ca["num_wires"] == cb["num_wires"]
da = [qlocal.compile_body(d["operations"]) for d in ca["definitions"]]
db = [qlocal.compile_body(d["operations"]) for d in cb["definitions"]]
rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 9); L = 512; ALL = (1 << L) - 1
n = ca["num_wires"]; init = [rng.getrandbits(L) for _ in range(n)]
sa, pa = qlocal.simulate(da, qlocal.compile_body(ca["operations"]), n, init, ALL)
sb, pb = qlocal.simulate(db, qlocal.compile_body(cb["operations"]), n, init, ALL)
print("top level on 512 random full states:", "same" if (sa, pa) == (sb, pb) else "DIFFERS")
