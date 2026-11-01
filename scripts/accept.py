"""Acceptance gate for changes to the circuit generator.

Runs a registry of named invariant checks, prints one row per check, and writes
the full result to verify-out/accept.json. Exits 0 only when no check fails. A
pending check waits on the unpublished benchmark and fails nothing.

    python3 scripts/accept.py --quick    # toy exactness plus cheap 256-bit structure
    python3 scripts/accept.py            # adds the 256-bit sampled check and a wider addend sample
    python3 scripts/accept.py --only I5 I6
"""
import argparse
import functools
import hashlib
import io
import json
import random
import shutil
import subprocess
import sys
import time
import traceback
from collections import Counter
from pathlib import Path
from typing import Callable, NamedTuple

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from qec.curve import SECP256K1, SECP256K1_G, SECP256K1_R
from qec.emit import write_circuit
from qec.point_add import ControlledPointAdd

MAX_WIRES = 4096
MAX_UPLOAD_BYTES = 104_857_600
T_PER_CCX = 7
ALLOWED_GATES = {"X", "CX", "CCX", "SWAP"}
REFERENCE = {"commit": "8ca9b0e", "max_t_over_1G_2G_3G": 50_465_989, "wires": 2311, "mutations_caught": 163}
TOY_ARGS = ["--toy-bits", "4", "6", "7", "8", "9", "10", "--toy-per-size", "6"]
EMITS = {
    "toy127": ("toy127_circuit.json", ["--curve", "toy:127", "--addend-index", "5"]),
    "toy127.hier": ("toy127_circuit.hier.json", ["--curve", "toy:127", "--addend-index", "5", "--hier"]),
    "secp256k1_G.hier": ("secp256k1_G_circuit.hier.json", ["--hier"]),
}

PASS, FAIL, PENDING = "pass", "fail", "pending"


class Result(NamedTuple):
    status: str
    note: str
    evidence: dict


class Proc(NamedTuple):
    rc: int
    stdout: str
    seconds: float
    log: str


class Emitted(NamedTuple):
    path: Path
    meta: dict
    summary: dict


def summarize(text):
    """Counts of a v2-ir circuit JSON, with every call expanded, plus a digest per definition."""
    doc = json.loads(text)
    defs = doc.get("definitions", [])
    seen = set()

    def expand(body, done):
        kinds = Counter()
        for ins in body:
            if "call" in ins:
                kinds.update(done[ins["call"]])
            else:
                kinds[ins["op"]] += 1
                seen.add(ins["op"])
        return kinds

    expanded = []
    for d in defs:
        expanded.append(expand(d["operations"], expanded))
    top_own = Counter(ins["op"] for ins in doc["operations"] if "call" not in ins)
    total = expand(doc["operations"], expanded)
    digest = lambda obj: hashlib.sha256(json.dumps(obj, separators=(",", ":")).encode()).hexdigest()
    return {
        "num_wires": doc["num_wires"],
        "registers": {r["name"]: [r["role"], r["wires"]] for r in doc["registers"]},
        "bytes": len(text.encode()),
        "gates": dict(sorted(total.items())),
        "gate_names_seen": sorted(seen),
        "ops": sum(total.values()),
        "ccx": total["CCX"],
        "t": T_PER_CCX * total["CCX"],
        "top_ops": sum(top_own.values()),
        "top_ccx": top_own["CCX"],
        "definitions": len(defs),
        "definition_digests": [digest(d) for d in defs],
    }


def rel(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def label(k):
    return f"[{k}]G" if k < 4 else f"[0x{k:064x}]G"


def lean(summary):
    return {k: v for k, v in summary.items() if k not in ("registers", "definition_digests")}


class Run:
    """The sources checks read. Each is built once per run and shared."""

    def __init__(self, out, quick, seed, sample):
        self.out, self.quick, self.seed, self.sample = out, quick, seed, sample
        self.work = out / "accept"
        self.git_before = _git_status()

    def script(self, tag, name, *args):
        t0 = time.time()
        p = subprocess.run([sys.executable, "-I", str(ROOT / "scripts" / name), *map(str, args)], cwd=ROOT, capture_output=True, text=True)
        log = self.work / "logs" / f"{tag}.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(f"$ {name} {' '.join(map(str, args))}\nexit {p.returncode}\n--- stdout\n{p.stdout}--- stderr\n{p.stderr}")
        return Proc(p.returncode, p.stdout, round(time.time() - t0, 1), rel(log))

    @functools.cache
    def emitted(self, name):
        filename, args = EMITS[name]
        path = self.work / "emit" / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        proc = self.script(f"emit-{name}", "emit.py", *args, "--out", path)
        if proc.rc:
            raise RuntimeError(f"scripts/emit.py exited {proc.rc}, see {proc.log}")
        meta = json.loads(path.with_name(filename.split(".json")[0] + ".meta.json").read_text())
        return Emitted(path, meta, summarize(path.read_text()))

    @functools.cache
    def toy(self, hier):
        out = self.work / ("toy_exhaustive_hier.json" if hier else "toy_exhaustive.json")
        proc = self.script("check-toy-hier" if hier else "check-toy", "check.py", *TOY_ARGS, *(["--hier"] if hier else []), "--out", out)
        return proc, json.loads(out.read_text()) if out.exists() else []

    @functools.cache
    def secp_sampled(self):
        out = self.work / "secp256k1_sampled.json"
        proc = self.script("check-secp", "check.py", "--secp", "--secp-random", 120, "--seed", 2026, "--out", out)
        return proc, json.loads(out.read_text())[0] if out.exists() else None

    @functools.cache
    def addend(self, k):
        """[k]G built hierarchically in process with the code scripts/emit.py runs."""
        gen = ControlledPointAdd(SECP256K1, SECP256K1.mul(k, SECP256K1_G))
        buf = io.StringIO()
        write_circuit(buf, gen.build(hierarchical=True), gen.layout)
        return summarize(buf.getvalue())

    def addends(self):
        rng = random.Random(self.seed)
        return [1, 2, 3] + [rng.randrange(4, SECP256K1_R) for _ in range(self.sample)]


def i1_interface(run):
    problems, evidence = [], {}
    for name in ("toy127", "secp256k1_G.hier"):
        e = run.emitted(name)
        n, w = e.meta["n"], e.summary["num_wires"]
        want_regs = {"control": ["CONTROL", [0]], "accumulator": ["ACCUMULATOR", list(range(1, 2 * n + 1))],
                     "scratch": ["SCRATCH", list(range(2 * n + 1, w))]}
        want_layout = {"control": [0], "x": list(range(1, n + 1)), "y": list(range(n + 1, 2 * n + 1))}
        if e.summary["registers"] != want_regs:
            problems.append(f"{name}: registers are not control=[0], accumulator=1..{2 * n}, scratch={2 * n + 1}..{w - 1}")
        if e.meta["layout"] != want_layout:
            problems.append(f"{name}: layout is not control=[0], x=1..{n}, y={n + 1}..{2 * n}")
        if w > MAX_WIRES:
            problems.append(f"{name}: {w} wires > {MAX_WIRES}")
        evidence[name] = {"n": n, "num_wires": w, "file": rel(e.path)}
    w256 = evidence["secp256k1_G.hier"]["num_wires"]
    evidence.update(max_wires=MAX_WIRES, reference_wires=REFERENCE["wires"], problems=problems)
    note = "; ".join(problems) or f"control=[0], x=1..n, y=n+1..2n at n=7 and n=256; {w256} wires (reference {REFERENCE['wires']}, cap {MAX_WIRES})"
    return Result(FAIL if problems else PASS, note, evidence)


def i2_model(run):
    p = run.script("model-check", "model_check.py")
    rows = [json.loads(line) for line in p.stdout.splitlines() if line.strip()]
    small = [r for r in rows if r["q"] != "secp256k1"]
    bad = [r for r in small if not (r["worst_iterations"] <= r["bound"] and r["width_schedule_ok"] and r["multiplier_ok"])]
    secp = next(r for r in rows if r["q"] == "secp256k1")
    ok = p.rc == 0 and small and not bad and secp["multiplier_ok"]
    note = f"{len(small)} primes exhaustive, worst Euclid run within bound; secp256k1 bound {secp['bound']}, {secp['random_pairs']} random products"
    return Result(PASS if ok else FAIL, note if ok else f"exit {p.rc}, {len(bad)} failing primes", {"exit": p.rc, "rows": rows, "failing": bad, "log": p.log, "seconds": p.seconds})


def _toy_result(proc, rows, what):
    fails = sum(r["failures"] for r in rows)
    ok = proc.rc == 0 and rows and not fails
    lanes = sum(r["lanes"] for r in rows)
    note = f"{len(rows)} toy instances, {lanes} lanes, {fails} failures ({what})"
    evidence = {"exit": proc.rc, "instances": len(rows), "lanes": lanes, "failures": fails, "log": proc.log, "seconds": proc.seconds,
                "first_failures": [r for r in rows if r["failures"]][:3]}
    return Result(PASS if ok else FAIL, note, evidence)


def i2_toy(run):
    return _toy_result(*run.toy(False), "flat circuit, every point, c in {0,1}")


def i2_toy_hier(run):
    return _toy_result(*run.toy(True), "hierarchical emission expanded by qec/flatten.py")


def i2_components(run):
    p = run.script("component-check", "component_check.py")
    out = json.loads(p.stdout)
    parts = {k: v for k, v in out.items() if isinstance(v, dict)}
    bad = sorted(k for k, v in parts.items() if v["failures"])
    ok = p.rc == 0 and parts and not bad
    note = f"{len(parts)} 256-bit components on boundary inputs, {sum(v['lanes'] for v in parts.values())} lanes, failing: {bad or 'none'}"
    return Result(PASS if ok else FAIL, note, {"exit": p.rc, "components": len(parts), "failing": bad, "log": p.log, "seconds": p.seconds})


def i2_mutation(run):
    p = run.script("mutation-control", "mutation_control.py")
    out = json.loads(p.stdout) if p.rc == 0 else {}
    caught = out.get("caught", 0)
    note = f"{caught} of {out.get('trials')} CCX-to-CX demotions caught (reference {REFERENCE['mutations_caught']})"
    return Result(PASS if p.rc == 0 and caught else FAIL, note if p.rc == 0 else f"exit {p.rc}, see {p.log}", {"exit": p.rc, **out, "log": p.log})


def i2_secp(run):
    proc, row = run.secp_sampled()
    ok = proc.rc == 0 and row is not None and not row["failures"]
    note = f"{row['lanes']} lanes at A = G incl. inf, A, -A, -2A, 2A, 3A; {row['failures']} failures" if row else f"exit {proc.rc}, see {proc.log}"
    return Result(PASS if ok else FAIL, note, {"exit": proc.rc, "receipt": row, "log": proc.log, "seconds": proc.seconds})


def i3_gate_set(run):
    problems, evidence = [], {}
    for name in EMITS:
        seen = run.emitted(name).summary["gate_names_seen"]
        evidence[name] = seen
        if set(seen) - ALLOWED_GATES:
            problems.append(f"{name} emits {sorted(set(seen) - ALLOWED_GATES)}")
    proc, rows = run.toy(False)
    for r in rows:
        if sum(r["counts"].values()) != r["ops"]:
            problems.append(f"toy q={r['q']} A={r['addend']}: {r['ops'] - sum(r['counts'].values())} ops outside X, CX, CCX, SWAP")
    evidence["toy_instances_counted"] = len(rows)
    if not rows:
        problems.append(f"no toy rows to count, see {proc.log}")
    if not run.quick:
        proc, row = run.secp_sampled()
        if row is None or sum(row["counts"].values()) != row["ops"]:
            problems.append("256-bit flat build has ops outside X, CX, CCX, SWAP" if row else f"no 256-bit flat counts, see {proc.log}")
        evidence["secp256k1_flat_counts"] = row and row["counts"]
    evidence["problems"] = problems
    where = "emitted toy and 256-bit JSON, every toy instance" + ("" if run.quick else ", and the 256-bit flat build")
    return Result(FAIL if problems else PASS, "; ".join(problems) or f"only X, CX, CCX, SWAP, no measurement or reset, in {where}", evidence)


def i5_addend_independence(run):
    profiles = {"G via scripts/emit.py": run.emitted("secp256k1_G.hier").summary}
    ks = run.addends() if not run.quick else [1, 2, 3]
    profiles.update((label(k), run.addend(k)) for k in ks)
    ref = profiles["G via scripts/emit.py"]["definition_digests"]
    diffs = []
    for name, p in profiles.items():
        got = p["definition_digests"]
        if got != ref:
            first = next((i for i, (a, b) in enumerate(zip(got, ref)) if a != b), min(len(got), len(ref)))
            diffs.append({"addend": name, "definitions": len(got), "first_differing_definition": first})
    evidence = {"addends": list(profiles), "definitions": len(ref),
                "digest_of_definitions": hashlib.sha256("".join(ref).encode()).hexdigest(),
                "top_level_ccx": {name: p["top_ccx"] for name, p in profiles.items()}, "differences": diffs}
    if diffs:
        d = diffs[0]
        return Result(FAIL, f"{len(diffs)} addends change the definitions; {d['addend']} first differs at definition {d['first_differing_definition']}", evidence)
    return Result(PASS, f"all {len(ref)} definitions identical over {len(profiles)} builds; A reaches only the top level", evidence)


def i6_cost(run):
    rows = {label(k): lean(run.addend(k)) for k in run.addends()}
    ts = {name: r["t"] for name, r in rows.items()}
    worst, best = max(ts, key=ts.get), min(ts, key=ts.get)
    base = REFERENCE["max_t_over_1G_2G_3G"]
    problems = []
    g_file = run.emitted("secp256k1_G.hier").summary
    if (g_file["ops"], g_file["ccx"]) != (rows["[1]G"]["ops"], rows["[1]G"]["ccx"]):
        problems.append("in-process [1]G and scripts/emit.py G disagree")
    toy_hier, toy_flat = run.emitted("toy127.hier").summary, run.emitted("toy127").meta
    cross = {"toy127": {"hier_expanded": [toy_hier["ops"], toy_hier["ccx"]], "flat": [toy_flat["ops"], toy_flat["ccx"]]}}
    if cross["toy127"]["hier_expanded"] != cross["toy127"]["flat"]:
        problems.append("toy127 expanded hierarchical counts differ from the flat emission")
    if not run.quick:
        _, flat = run.secp_sampled()
        cross["secp256k1_G"] = {"hier_expanded": g_file["gates"], "flat": flat and flat["counts"]}
        if flat is None or g_file["gates"] != flat["counts"]:
            problems.append("256-bit G expanded hierarchical gate counts differ from the flat build")
    evidence = {"t_per_ccx": T_PER_CCX, "seed": run.seed, "sampled": run.sample, "per_addend": rows,
                "max": {"addend": worst, **rows[worst]}, "min": {"addend": best, **rows[best]},
                "reference_max_t": base, "max_t_delta_vs_reference": ts[worst] - base, "count_cross_checks": cross, "problems": problems}
    if problems:
        return Result(FAIL, "; ".join(problems), evidence)
    g = rows["[1]G"]
    note = (f"G: {g['ccx']:,} CCX, T {g['t']:,}, {g['ops']:,} ops, {g['num_wires']} wires; over {len(rows)} addends "
            f"max T {ts[worst]:,} ({ts[worst] - base:+,} vs reference {base:,}), min T {ts[best]:,}")
    return Result(PASS, note, evidence)


def i7_limits(run):
    g = run.emitted("secp256k1_G.hier")
    size = g.path.stat().st_size
    largest = max((run.addend(k)["bytes"], label(k)) for k in run.addends())
    wires = g.summary["num_wires"]
    problems = [m for ok, m in ((wires <= MAX_WIRES, f"{wires} wires > {MAX_WIRES}"),
                                (size < MAX_UPLOAD_BYTES, f"hierarchical JSON for G is {size:,} bytes >= {MAX_UPLOAD_BYTES:,}"),
                                (largest[0] < MAX_UPLOAD_BYTES, f"hierarchical JSON for {largest[1]} is {largest[0]:,} bytes")) if not ok]
    evidence = {"num_wires": wires, "max_wires": MAX_WIRES, "hier_json_bytes_G": size, "hier_json_file": rel(g.path),
                "largest_hier_json_over_addends": {"addend": largest[1], "bytes": largest[0]}, "max_upload_bytes": MAX_UPLOAD_BYTES, "problems": problems}
    note = f"{wires} of {MAX_WIRES} wires; hierarchical JSON {size:,} bytes for G, at most {largest[0]:,} over addends, cap {MAX_UPLOAD_BYTES:,}"
    return Result(FAIL if problems else PASS, "; ".join(problems) or note, evidence)


def i7_charge(run):
    return Result(PENDING, "needs the published benchmark: its charge formula is not public yet",
                  {"reason": "needs the published benchmark", "limit_as_reported": 33_610_176, "this_repo_as_reported": 19_492_771})


def r1_hygiene(run):
    after = _git_status()
    if run.git_before is None:
        return Result(PASS, "not a git checkout; nothing to compare", {})
    changed = sorted(set(after.splitlines()) ^ set(run.git_before.splitlines()))
    return Result(FAIL if changed else PASS, f"run changed {changed}" if changed else "git status is the same before and after the run",
                  {"changed": changed})


class Check(NamedTuple):
    id: str
    title: str
    run: Callable
    full_only: bool = False


CHECKS = (
    Check("I1", "interface: layout and wire count", i1_interface),
    Check("I2.model", "exactness: classical dialog model and iteration bound", i2_model),
    Check("I2.toy", "exactness: exhaustive gate-level check on toy curves", i2_toy),
    Check("I2.toy-hier", "exactness: same check on the hierarchical upload", i2_toy_hier),
    Check("I2.components", "exactness: 256-bit components on boundary inputs", i2_components),
    Check("I2.mutation", "exactness: the checker rejects mutated circuits", i2_mutation),
    Check("I2.secp", "exactness: 256-bit circuit on sampled and exceptional points", i2_secp, full_only=True),
    Check("I3", "gate set: X, CX, CCX, SWAP only", i3_gate_set),
    Check("I5", "addend independence of definitions", i5_addend_independence),
    Check("I6", "cost at G and over sampled addends (report)", i6_cost),
    Check("I7", "limits: wires and hierarchical upload", i7_limits),
    Check("I7.charge", "limits: benchmark charge", i7_charge),
    Check("R1", "run hygiene: no file git tracks changed", r1_hygiene),
)


def _git_status():
    p = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all"], cwd=ROOT, capture_output=True, text=True)
    return p.stdout if p.returncode == 0 else None


def _git_head():
    p = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True)
    return p.stdout.strip() if p.returncode == 0 else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="toy checks and cheap 256-bit structure checks only")
    ap.add_argument("--seed", type=int, default=2026, help="seed for the sampled addends")
    ap.add_argument("--sample", type=int, help="random addends besides [1]G, [2]G, [3]G (default 4 quick, 16 full)")
    ap.add_argument("--only", nargs="+", metavar="ID", choices=[c.id for c in CHECKS], help="run only these checks")
    ap.add_argument("--out", type=Path, default=ROOT / "verify-out")
    args = ap.parse_args()
    sample = args.sample if args.sample is not None else (4 if args.quick else 16)
    run = Run(args.out.resolve(), args.quick, args.seed, sample)
    shutil.rmtree(run.work, ignore_errors=True)
    run.work.mkdir(parents=True)
    t_start = time.time()

    results, skipped = [], []
    for check in CHECKS:
        if args.only and check.id not in args.only:
            continue
        if check.full_only and args.quick:
            skipped.append(check.id)
            continue
        t0 = time.time()
        try:
            r = check.run(run)
        except Exception as e:
            r = Result(FAIL, f"check crashed: {e!r}"[:300], {"traceback": traceback.format_exc()})
        secs = round(time.time() - t0, 1)
        results.append({"id": check.id, "title": check.title, "status": r.status, "seconds": secs, "note": r.note, "evidence": r.evidence})
        print(f"{check.id:<14}{r.status.upper():<9}{secs:>6.1f}s  {check.title}\n{'':<31}{r.note}", flush=True)

    tally = Counter(r["status"] for r in results)
    failed = [r["id"] for r in results if r["status"] == FAIL]
    pending = [r["id"] for r in results if r["status"] == PENDING]
    report = {"gate": FAIL if failed else PASS, "mode": "quick" if args.quick else "full", "seed": args.seed, "sample": sample,
              "git_head": _git_head(), "python": sys.version.split()[0], "seconds": round(time.time() - t_start, 1),
              "counts": dict(tally), "failed": failed, "pending": pending, "not_run_in_quick_mode": skipped, "results": results}
    (run.out / "accept.json").write_text(json.dumps(report, indent=1))
    print(f"\ngate {report['gate'].upper()}: {tally[PASS]} pass, {tally[FAIL]} fail, {tally[PENDING]} pending {pending or ''}"
          f"{f'; not run in --quick: {skipped}' if skipped else ''}; {report['seconds']}s; {rel(run.out / 'accept.json')}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
