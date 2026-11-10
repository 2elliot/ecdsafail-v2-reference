"""[integrity copy: adds --cases FILE, a fixed list of lanes from hard_cases.py, and per-lane failure ids]
Local stand-in for the board's check, for circuits we build or build on.

    python3 -I qlocal.py <circuit-or-method.json> [--draws N] [--seed S] [--reverse] [--no-sim]

Prints limits (charge as the member guide's qcheck.py computes it), T, ticks and
qubit-time (ASAP over the 15-gate Clifford+T Toffoli, which reproduced the board
exactly on reference-c4 and v2-tape), then a bit-sliced simulation of the 10
boundary cases, adversarial inputs and N random points. Exits 0 only if every
lane passes and every limit holds. This is a sampled test, like the board's.
"""
import json, math, os, random, sys

P = 2**256 - 2**32 - 977
GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
INF = (0, 0)
LIMITS = {"wires": 4096, "upload bytes": 104_857_600, "expanded ops charge": 21_696_718,
          "T": 3_250_000_000, "call depth": 32, "records": 10_000_000, "listed call wires": 10_000_000}


def ec_add(A, B):
    if A == INF: return B
    if B == INF: return A
    (x1, y1), (x2, y2) = A, B
    if x1 == x2:
        if (y1 + y2) % P == 0: return INF
        lam = 3 * x1 * x1 * pow(2 * y1, -1, P) % P
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, P) % P
    x3 = (lam * lam - x1 - x2) % P
    return (x3, (lam * (x1 - x3) - y1) % P)


def on_curve(x):
    rhs = (x * x * x + 7) % P
    y = pow(rhs, (P + 1) // 4, P)
    return y if y * y % P == rhs else None


G = (GX, GY)
NEG_G = (GX, P - GY)
G2 = ec_add(G, G)
NEG_G2 = (G2[0], P - G2[1])


# ---------------------------------------------------------------- structure and limits
def compile_body(ops):
    out = []
    for r in ops:
        if "call" in r:
            out.append((9, r["call"], tuple(r["wires"]), bool(r.get("inverse", False))))
            continue
        o = r["op"]
        if o == "X": out.append((0, r["target"]))
        elif o == "CX": out.append((1, r["control"], r["target"]))
        elif o == "CCX": out.append((2, r["control1"], r["control2"], r["target"]))
        elif o == "SWAP": out.append((3, r["wire1"], r["wire2"]))
        elif o == "Z": out.append((4, r["target"]))
        elif o == "CZ": out.append((5, r["control"], r["target"]))
        elif o == "CCZ": out.append((6, r["control1"], r["control2"], r["target"]))
        elif o == "NEG": out.append((7,))
        elif o == "CONTROLLED_NEG": out.append((8, r["control"]))
        else: raise SystemExit(f"unknown op {o}")
    return out


def validate(circuit):
    """Structural rules the board's parser enforces that we know of. Returns a list of problems."""
    problems = []
    defs = circuit.get("definitions", [])
    fields = ("target", "control", "control1", "control2", "wire1", "wire2")

    def check_body(ops, arity, where, idx):
        used = set()
        for r in ops:
            if "call" in r:
                j, ws = r["call"], r["wires"]
                if idx is not None and j >= idx:
                    problems.append(f"{where}: calls definition {j}, not an earlier one")
                elif j < len(defs) and len(ws) != defs[j]["arity"]:
                    problems.append(f"{where}: call to definition {j} lists {len(ws)} wires, arity {defs[j]['arity']}")
                if len(set(ws)) != len(ws):
                    problems.append(f"{where}: call to definition {j} repeats a wire")
                used.update(ws)
            else:
                ws = [r[f] for f in fields if f in r]
                if len(set(ws)) != len(ws):
                    problems.append(f"{where}: {r['op']} repeats a wire")
                used.update(ws)
        out = sorted(w for w in used if not 0 <= w < arity)
        if out:
            problems.append(f"{where}: wires out of range {out[:5]}")
        return used

    for i, d in enumerate(defs):
        used = check_body(d["operations"], d["arity"], f"definition {i}", i)
        missing = sorted(set(range(d["arity"])) - used)
        if missing:
            problems.append(f"definition {i}: definition wire {missing[0]} is never touched by its body"
                            + (f" ({len(missing)} such wires)" if len(missing) > 1 else ""))
    check_body(circuit["operations"], circuit["num_wires"], "top level", None)
    return problems


def structure(circuit):
    """Counts as qcheck.py's counts_of does: a definition may call only earlier ones."""
    done, records, listed = [], 0, 0

    def body(ops, idx):
        nonlocal records, listed
        t = {"ops": 0, "ccx": 0, "calls": 0, "call_wires": 0, "depth": 0}
        for r in ops:
            records += 1
            if r[0] == 9:
                if idx is not None and r[1] >= idx: raise SystemExit(f"definition {idx} calls {r[1]}, not an earlier one")
                d = done[r[1]]
                listed += len(r[2])
                for k in ("ops", "ccx", "calls", "call_wires"): t[k] += d[k]
                t["calls"] += 1; t["call_wires"] += len(r[2]); t["depth"] = max(t["depth"], d["depth"] + 1)
            else:
                t["ops"] += 1; t["ccx"] += r[0] in (2, 6)
        return t

    defs = [compile_body(d["operations"]) for d in circuit.get("definitions", [])]
    for i, d in enumerate(defs):
        done.append(body(d, i))
    top = compile_body(circuit["operations"])
    t = body(top, None)
    ps = t["calls"] * 1_329_430 + t["call_wires"] * 17_360 + 918_090
    t["charge"] = t["ops"] + math.ceil(ps / 20_710)
    t["records"], t["listed"] = records, listed
    return defs, top, t


# ---------------------------------------------------------------- ticks
def schedule(defs, top, nw):
    r = [0] * nw

    def run(body, wmap, inv):
        for g in (reversed(body) if inv else body):
            k = g[0]
            if k == 9:
                run(defs[g[1]], [wmap[w] for w in g[2]] if wmap is not None else g[2], inv ^ g[3]); continue
            # Calibrated on X, CX, CCX (and SWAP off the critical path). Z, CZ, CCZ, NEG and
            # CONTROLLED_NEG timings below are guesses: no calibrated circuit uses them.
            if k == 7: continue
            w = [wmap[x] for x in g[1:]] if wmap is not None else g[1:]
            if k in (0, 4, 8): r[w[0]] += 1
            elif k in (1, 5):
                t = max(r[w[0]], r[w[1]]) + 1; r[w[0]] = r[w[1]] = t
            elif k == 3:
                t = max(r[w[0]], r[w[1]]) + 3; r[w[0]] = r[w[1]] = t
            else:
                a, b, t = w; ra, rb, rt = r[a], r[b], r[t]
                rt += 1; x = max(rb, rt) + 1; rb = rt = x; rt += 1
                x = max(ra, rt) + 1; ra = rt = x; rt += 1
                x = max(rb, rt) + 1; rb = rt = x; rt += 1
                x = max(ra, rt) + 1; ra = rt = x
                rb += 1; rt += 2
                x = max(ra, rb) + 1; ra = rb = x; ra += 1; rb += 1
                x = max(ra, rb) + 1; ra = rb = x
                r[a], r[b], r[t] = ra, rb, rt

    run(top, None, False)
    return max(r)


# ---------------------------------------------------------------- simulation
def simulate(defs, top, nw, init, lanes_mask, inv=False):
    s = list(init)
    phase = 0
    ALL = lanes_mask

    def run(body, wmap, inv):
        nonlocal phase
        for g in (reversed(body) if inv else body):
            k = g[0]
            if k == 9:
                run(defs[g[1]], [wmap[w] for w in g[2]] if wmap is not None else g[2], inv ^ g[3]); continue
            if k == 7: phase ^= ALL; continue
            w = [wmap[x] for x in g[1:]] if wmap is not None else g[1:]
            if k == 2: s[w[2]] ^= s[w[0]] & s[w[1]]
            elif k == 1: s[w[1]] ^= s[w[0]]
            elif k == 0: s[w[0]] ^= ALL
            elif k == 3: s[w[0]], s[w[1]] = s[w[1]], s[w[0]]
            elif k == 4: phase ^= s[w[0]]
            elif k == 5: phase ^= s[w[0]] & s[w[1]]
            elif k == 6: phase ^= s[w[0]] & s[w[1]] & s[w[2]]
            elif k == 8: phase ^= s[w[0]]

    run(top, None, inv)
    return s, phase


def pack(values, nbits):
    """values[l] -> per-bit lane masks."""
    masks = [0] * nbits
    for l, v in enumerate(values):
        while v:
            low = v & -v
            masks[low.bit_length() - 1] |= 1 << l
            v ^= low
    return masks


def unpack(masks, nlanes):
    vals = [0] * nlanes
    for b, m in enumerate(masks):
        while m:
            low = m & -m
            vals[low.bit_length() - 1] |= 1 << b
            m ^= low
    return vals


def random_point(rng):
    while True:
        x = rng.randrange(P)
        y = on_curve(x)
        if y is not None:
            return (x, y if rng.getrandbits(1) else (P - y) % P)


def adversarial_points():
    """Points whose GCD inputs are long: x_R - G_x = 2^k (k odd near 255) and nearby shapes."""
    pts = []
    for k in (255, 253, 251, 254, 252):
        for delta in (1 << k, P - (1 << k)):
            x = (GX + delta) % P
            y = on_curve(x)
            if y is not None:
                pts.append((x, y))
    for x in (1, 2, 3, P - 1, P - 2, (P - 1) // 2, (P + 1) // 2, 1 << 255, (1 << 256) - (1 << 32) - 978):
        y = on_curve(x % P)
        if y is not None:
            pts.append((x % P, y))
    return pts


def build_cases(draws, seed, fixed=None):
    rng = random.Random(seed)
    cases = []
    if fixed is not None:
        for R in (INF, G, NEG_G, G2, NEG_G2):
            for c in (0, 1):
                cases.append(("B", c, R))
        for lane in fixed:
            R = (int(lane["x"], 16), int(lane["y"], 16))
            want = (int(lane["ex"], 16), int(lane["ey"], 16))
            assert R == INF or on_curve(R[0]) is not None and (R[1] * R[1] - R[0] ** 3 - 7) % P == 0, lane["id"]
            assert (ec_add(R, G) if lane["c"] else R) == want, ("expected output disagrees", lane["id"])
            cases.append(("H", lane["c"], R, lane["id"]))
        for _ in range(draws):
            cases.append(("R", rng.getrandbits(1), random_point(rng)))
        return cases, rng
    for R in (INF, G, NEG_G, G2, NEG_G2):
        for c in (0, 1):
            cases.append(("B", c, R))
    for R in adversarial_points():
        for c in (0, 1):
            cases.append(("A", c, R))
    for _ in range(draws):
        cases.append(("R", rng.getrandbits(1), random_point(rng)))
    return cases, rng


def main(argv):
    path = argv[0]
    draws = int(argv[argv.index("--draws") + 1]) if "--draws" in argv else 10_000
    seed = int(argv[argv.index("--seed") + 1]) if "--seed" in argv else int.from_bytes(os.urandom(8), "big")
    with open(path) as fh:
        doc = json.load(fh)
    circuit = doc.get("circuit", doc)
    size = os.path.getsize(path)
    nw = circuit["num_wires"]
    problems = validate(circuit)
    if problems:
        print(f"file {path}")
        for p in problems[:10]:
            print("  INVALID:", p)
        print(f"INVALID: {len(problems)} structural problems; the board would refuse this file")
        return 1
    defs, top, t = structure(circuit)
    T = 7 * t["ccx"]
    ticks = schedule(defs, top, nw)
    print(f"file {path}")
    print(f"wires {nw:,}  CCX {t['ccx']:,}  T {T:,}  ticks {ticks:,}  qubit-time {nw * ticks:,}")
    print(f"vs reference: T {T - 50_465_933:+,}  qubit-time {nw * ticks - 164_135_627_418:+,}")
    good = True
    for name, val in (("wires", nw), ("upload bytes", size), ("expanded ops charge", t["charge"]), ("T", T),
                      ("call depth", t["depth"]), ("records", t["records"]), ("listed call wires", t["listed"])):
        lim = LIMITS[name]
        flag = "OVER" if val > lim else ("near" if val > 0.9 * lim else "ok")
        good &= val <= lim
        print(f"  limit {name}: {val:,} of {lim:,} {flag}")
    if "--no-sim" in argv:
        print("simulation skipped"); return 0 if good else 1

    roles = {}
    for reg in circuit["registers"]:
        roles.setdefault(reg["role"], []).extend(reg["wires"])
    ctrl = roles["CONTROL"]; acc = roles["ACCUMULATOR"]
    assert len(ctrl) == 1 and len(acc) == 512, "need 1 control and 512 accumulator wires"
    borrowed = [w for role in ("BORROWED", "ADDEND", "ADDRESS") for w in roles.get(role, [])]
    fixed = None
    if "--cases" in argv:
        with open(argv[argv.index("--cases") + 1]) as fh:
            fixed = json.load(fh)["lanes"]
        if "--draws" not in argv:
            draws = 0
    cases, rng = build_cases(draws, seed, fixed)
    L = len(cases)
    ALL = (1 << L) - 1
    init = [0] * nw
    xs = pack([cs[2][0] for cs in cases], 256)
    ys = pack([cs[2][1] for cs in cases], 256)
    for i in range(256):
        init[acc[i]] = xs[i]; init[acc[256 + i]] = ys[i]
    init[ctrl[0]] = sum(cs[1] << l for l, cs in enumerate(cases))
    noise = {w: rng.getrandbits(L) for w in borrowed}
    for w, m in noise.items():
        init[w] = m
    out, phase = simulate(defs, top, nw, init, ALL)
    gotx = unpack([out[acc[i]] for i in range(256)], L)
    goty = unpack([out[acc[256 + i]] for i in range(256)], L)
    bad = 0
    coord = 0
    for l, cs in enumerate(cases):
        c, R = cs[1], cs[2]
        want = ec_add(R, G) if c else R
        if (gotx[l], goty[l]) != want:
            coord |= 1 << l
    keep = set(ctrl) | set(borrowed)
    preserved = 0
    for w in keep:
        preserved |= out[w] ^ init[w]
    garbage = 0
    accset = set(acc)
    for w in range(nw):
        if w not in accset and w not in keep:
            garbage |= out[w]
    phase_bad = phase ^ (ALL if phase & 1 else 0)
    fail = coord | preserved | garbage | phase_bad
    kinds = {"B": [0, 0], "A": [0, 0], "H": [0, 0], "R": [0, 0]}
    for l, cs in enumerate(cases):
        kinds[cs[0]][0] += 1; kinds[cs[0]][1] += fail >> l & 1
    failing_ids = []
    for l, cs in enumerate(cases):
        if fail >> l & 1:
            flags = [nm for nm, m_ in (("coord", coord), ("preserved", preserved), ("garbage", garbage), ("phase", phase_bad)) if m_ >> l & 1]
            failing_ids.append((cs[3] if len(cs) > 3 else f"{cs[0]}{l}", cs[1], flags))
    if failing_ids:
        print("  FAILING LANES:", failing_ids[:50])
    for k, label in (("B", "boundary"), ("A", "adversarial"), ("H", "hard fixed set"), ("R", "random")):
        print(f"  {label}: {kinds[k][1]} of {kinds[k][0]} failing")
    for name, m in (("coordinate_mismatch", coord), ("preserved_register_corruption", preserved),
                    ("workspace_garbage", garbage), ("logical_phase_mismatch", phase_bad)):
        if m: print(f"  {name}: {bin(m).count('1')} lanes")
    if "--reverse" in argv:
        back, _ = simulate(defs, top, nw, out, ALL, inv=True)
        rev = sum(1 for w in range(nw) if back[w] != init[w])
        print(f"  reverse run: {rev} wires differ from the input")
        good &= rev == 0
    good &= fail == 0
    print(f"seed {seed}: " + ("PASSED, every limit holds" if good else "FAILED"))
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
