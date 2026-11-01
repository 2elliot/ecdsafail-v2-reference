"""Acceptance gate for changes to the circuit generator.

Runs a registry of named invariant checks, prints one row per check, and writes
the full result to verify-out/accept.json. Exits 0 only when no check fails. A
pending check waits on the unpublished benchmark and fails nothing. The
checker scripts run as parallel subprocesses, all started before the first row.

    python3 scripts/accept.py --quick    # toy exactness, 256-bit structure, 256-bit upload simulated at 2 addends
    python3 scripts/accept.py            # simulates the upload at every sampled addend, and the flat circuit at G
    python3 scripts/accept.py --only I5 I6
"""
import argparse
import functools
import hashlib
import io
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
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
SECP_SPECIAL_POINTS = 6
G_MULTIPLES = (1, 2, 3)
SECP_RANDOM_POINTS = {True: 8, False: 120}
SECP_WORKERS = min(6, os.cpu_count() or 1)
SUMMARY_WORKERS = min(3, os.cpu_count() or 1)
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
        "sha256": hashlib.sha256(text.encode()).hexdigest(),
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
    return f"[{k}]G" if k <= G_MULTIPLES[-1] else f"[0x{k:064x}]G"


def lean(summary):
    return {k: v for k, v in summary.items() if k not in ("registers", "definition_digests", "layout", "sha256")}


def json_lines(text):
    rows = []
    for line in text.splitlines():
        try:
            rows.append(json.loads(line))
        except ValueError:
            pass
    return rows


class Job:
    """One checker script running in a subprocess, collected on first use."""

    def __init__(self, argv, log):
        self.argv, self.log = argv, log
        self.out, self.err = tempfile.TemporaryFile("w+"), tempfile.TemporaryFile("w+")
        self.p = subprocess.Popen([sys.executable, "-I", str(ROOT / "scripts" / argv[0]), *map(str, argv[1:])],
                                  cwd=ROOT, stdout=self.out, stderr=self.err, text=True)

    @functools.cached_property
    def proc(self):
        rc = self.p.wait()
        self.out.seek(0)
        self.err.seek(0)
        stdout, stderr = self.out.read(), self.err.read()
        self.log.write_text(f"$ {' '.join(map(str, self.argv))}\nexit {rc}\n--- stdout\n{stdout}--- stderr\n{stderr}")
        return Proc(rc, stdout, rel(self.log))


class Run:
    """The sources checks read. Each is built once per run and shared."""

    def __init__(self, out, quick, seed, sample):
        self.out, self.quick, self.seed, self.sample = out, quick, seed, sample
        self.work = out / "accept"
        self.git_before = _git_status()
        rng = random.Random(seed)
        self.addends = [*G_MULTIPLES, *(rng.randrange(G_MULTIPLES[-1] + 1, SECP256K1_R) for _ in range(sample))]
        self.simulated = [1, self.addends[-1]] if quick else self.addends
        self.secp_points = SECP_RANDOM_POINTS[quick]
        w = self.work
        chunks = [c for c in (self.simulated[i::SECP_WORKERS] for i in range(SECP_WORKERS)) if c]
        self.specs = {
            "model": [["model_check.py"]],
            "toy": [["check.py", *TOY_ARGS, "--out", w / "toy_exhaustive.json"]],
            "toy-hier": [["check.py", *TOY_ARGS, "--hier", "--out", w / "toy_exhaustive_hier.json"]],
            "components": [["component_check.py"]],
            "mutation": [["mutation_control.py"]],
            "secp-hier": [["check.py", "--secp", "--hier", "--secp-addend", *map(hex, c), "--secp-random", self.secp_points,
                           "--seed", seed, "--out", w / f"secp256k1_hier_{i}.json"] for i, c in enumerate(chunks)],
            **{f"emit:{name}": [["emit.py", *args, "--out", w / "emit" / filename]] for name, (filename, args) in EMITS.items()},
        }
        self.specs["addends"] = []
        if not quick:
            self.specs["secp-flat"] = [["check.py", "--secp", "--secp-random", 120, "--seed", 2026, "--out", w / "secp256k1_sampled.json"]]
        self.jobs = {}
        self.pool = None
        self.summaries = {}

    def start(self, group):
        if group == "addends":
            if self.pool is None:
                self.pool = ProcessPoolExecutor(SUMMARY_WORKERS)
                self.summaries = {k: self.pool.submit(addend_summary, k) for k in self.addends}
            return []
        if group not in self.jobs:
            (self.work / "logs").mkdir(parents=True, exist_ok=True)
            specs = self.specs[group]
            self.jobs[group] = [Job(argv, self.work / "logs" / f"{group.replace(':', '-')}{f'-{i}' if len(specs) > 1 else ''}.log")
                                for i, argv in enumerate(specs)]
        return self.jobs[group]

    def procs(self, group):
        return [job.proc for job in self.start(group)]

    def wait_all(self):
        for jobs in self.jobs.values():
            for job in jobs:
                job.proc
        if self.pool:
            self.pool.shutdown()

    @functools.cache
    def emitted(self, name):
        filename, _ = EMITS[name]
        path = self.work / "emit" / filename
        proc, = self.procs(f"emit:{name}")
        if proc.rc:
            raise RuntimeError(f"scripts/emit.py exited {proc.rc}, see {proc.log}")
        meta = json.loads(path.with_name(filename.split(".json")[0] + ".meta.json").read_text())
        return Emitted(path, meta, summarize(path.read_text()))

    @functools.cache
    def toy(self, hier):
        proc, = self.procs("toy-hier" if hier else "toy")
        out = self.work / ("toy_exhaustive_hier.json" if hier else "toy_exhaustive.json")
        return proc, json.loads(out.read_text()) if out.exists() else []

    @functools.cache
    def secp_flat(self):
        proc, = self.procs("secp-flat")
        out = self.work / "secp256k1_sampled.json"
        return proc, json.loads(out.read_text())[0] if out.exists() else None

    @functools.cache
    def secp_hier(self):
        """(procs, {k: receipt}) of the hierarchical 256-bit simulation."""
        procs = self.procs("secp-hier")
        return procs, {int(r["k"], 16): r for p in procs for r in json_lines(p.stdout)}

    @functools.cache
    def addend(self, k):
        self.start("addends")
        return self.summaries[k].result()


def addend_summary(k):
    """[k]G built hierarchically with the code scripts/emit.py runs, and summarized."""
    gen = ControlledPointAdd(SECP256K1, SECP256K1_G if k == 1 else SECP256K1.mul(k, SECP256K1_G))
    buf = io.StringIO()
    write_circuit(buf, gen.build(hierarchical=True), gen.layout)
    return {**summarize(buf.getvalue()), "layout": gen.layout}


def _interface_problems(name, n, registers, layout, wires):
    want_regs = {"control": ["CONTROL", [0]], "accumulator": ["ACCUMULATOR", list(range(1, 2 * n + 1))],
                 "scratch": ["SCRATCH", list(range(2 * n + 1, wires))]}
    want_layout = {"control": [0], "x": list(range(1, n + 1)), "y": list(range(n + 1, 2 * n + 1))}
    problems = []
    if registers != want_regs:
        problems.append(f"{name}: registers are not control=[0], accumulator=1..{2 * n}, scratch={2 * n + 1}..{wires - 1}")
    if layout != want_layout:
        problems.append(f"{name}: layout is not control=[0], x=1..{n}, y={n + 1}..{2 * n}")
    if wires > MAX_WIRES:
        problems.append(f"{name}: {wires} wires > {MAX_WIRES}")
    return problems


def i1_interface(run):
    problems, evidence = [], {}
    for name in ("toy127", "secp256k1_G.hier"):
        e = run.emitted(name)
        problems += _interface_problems(name, e.meta["n"], e.summary["registers"], e.meta["layout"], e.summary["num_wires"])
        evidence[name] = {"n": e.meta["n"], "num_wires": e.summary["num_wires"], "file": rel(e.path)}
    for k in run.addends:
        s = run.addend(k)
        problems += _interface_problems(label(k), 256, s["registers"], s["layout"], s["num_wires"])
    wires = {label(k): run.addend(k)["num_wires"] for k in run.addends}
    evidence.update(addend_wires=wires, max_wires=MAX_WIRES, reference_wires=REFERENCE["wires"], problems=problems)
    w256 = evidence["secp256k1_G.hier"]["num_wires"]
    note = "; ".join(problems) or (f"control=[0], x=1..n, y=n+1..2n at n=7, n=256 and all {len(wires)} built addends; "
                                   f"{w256} wires, at most {max(wires.values())} over addends (reference {REFERENCE['wires']}, cap {MAX_WIRES})")
    return Result(FAIL if problems else PASS, note, evidence)


def i2_model(run):
    p, = run.procs("model")
    rows = json_lines(p.stdout)
    small = [r for r in rows if r["q"] != "secp256k1"]
    bad = [r for r in small if not (r["worst_iterations"] <= r["bound"] and r["width_schedule_ok"] and r["multiplier_ok"])]
    secp = next(r for r in rows if r["q"] == "secp256k1")
    ok = p.rc == 0 and small and not bad and secp["multiplier_ok"]
    note = f"{len(small)} primes exhaustive, worst Euclid run within bound; secp256k1 bound {secp['bound']}, {secp['random_pairs']} random products"
    return Result(PASS if ok else FAIL, note if ok else f"exit {p.rc}, {len(bad)} failing primes", {"exit": p.rc, "rows": rows, "failing": bad, "log": p.log})


def _toy_result(proc, rows, what, hier):
    fails = sum(r["failures"] for r in rows)
    lanes = sum(r["lanes"] for r in rows)
    problems = [f"check.py exited {proc.rc} after {len(rows)} toy instances, see {proc.log}"] * bool(proc.rc and not fails)
    if not rows:
        problems.append("no toy instances")
    mismatched = [r for r in rows if hier and not r.get("expansion_matches_flatten")]
    if mismatched:
        problems.append(f"{len(mismatched)} expansions differ from flatten.operations_json")
    note = f"{len(rows)} toy instances, {lanes} lanes, {fails} failures ({what})" + "".join(f"; {m}" for m in problems)
    evidence = {"exit": proc.rc, "instances": len(rows), "lanes": lanes, "failures": fails, "log": proc.log, "problems": problems,
                "first_failures": [r for r in rows if r["failures"]][:3]}
    return Result(FAIL if fails or problems else PASS, note, evidence)


def i2_toy(run):
    return _toy_result(*run.toy(False), "flat circuit, every point, c in {0,1}", False)


def i2_toy_hier(run):
    return _toy_result(*run.toy(True), "hierarchical emission expanded by qec/flatten.py's rules", True)


def i2_components(run):
    p, = run.procs("components")
    out = json.loads(p.stdout)
    parts = {k: v for k, v in out.items() if isinstance(v, dict)}
    bad = sorted(k for k, v in parts.items() if v["failures"])
    ok = p.rc == 0 and parts and not bad
    note = f"{len(parts)} 256-bit components on boundary inputs, {sum(v['lanes'] for v in parts.values())} lanes, failing: {bad or 'none'}"
    return Result(PASS if ok else FAIL, note, {"exit": p.rc, "components": len(parts), "failing": bad, "log": p.log})


def i2_mutation(run):
    p, = run.procs("mutation")
    out = json.loads(p.stdout) if p.rc == 0 else {}
    caught = out.get("caught", 0)
    note = f"{caught} of {out.get('trials')} CCX-to-CX demotions caught (reference {REFERENCE['mutations_caught']})"
    return Result(PASS if p.rc == 0 and caught else FAIL, note if p.rc == 0 else f"exit {p.rc}, see {p.log}", {"exit": p.rc, **out, "log": p.log})


def i2_secp(run):
    procs, rows = run.secp_hier()
    lanes = 2 * (SECP_SPECIAL_POINTS + run.secp_points)
    failing = {label(k): f"{r['failures']} of {r['lanes']} lanes" for k, r in rows.items() if r["failures"]}
    problems = [f"check.py exited {p.rc}, see {p.log}" for p in procs if p.rc and not any(r["failures"] for r in json_lines(p.stdout))]
    problems += [f"{label(k)} not simulated" for k in run.simulated if k not in rows]
    problems += [f"{label(k)}: {r['lanes']} lanes, expected {lanes}" for k, r in rows.items() if r["lanes"] != lanes]
    g_file = run.emitted("secp256k1_G.hier").summary["sha256"]
    if 1 in rows and rows[1]["json_sha256"] != g_file:
        problems.append("the simulated G JSON is not the bytes scripts/emit.py wrote")
    evidence = {"addends": [label(k) for k in run.simulated], "lanes_per_addend": lanes, "failing": failing, "problems": problems,
                "logs": [p.log for p in procs], "first_failures": {label(k): r["first_failures"] for k, r in rows.items() if r["failures"]}}
    if failing or problems:
        return Result(FAIL, "; ".join([f"{a} fails {f}" for a, f in failing.items()] + problems), evidence)
    return Result(PASS, f"hierarchical upload expanded by flatten's rules at {len(rows)} addends incl. G, {lanes} lanes each "
                        f"(inf, A, -A, -2A, 2A, 3A and {run.secp_points} random R, c in {{0,1}}); 0 failures", evidence)


def i2_secp_flat(run):
    proc, row = run.secp_flat()
    ok = proc.rc == 0 and row is not None and not row["failures"]
    note = f"{row['lanes']} lanes at A = G incl. inf, A, -A, -2A, 2A, 3A; {row['failures']} failures" if row else f"exit {proc.rc}, see {proc.log}"
    return Result(PASS if ok else FAIL, note, {"exit": proc.rc, "receipt": row, "log": proc.log})


def i3_gate_set(run):
    problems, evidence = [], {}
    sources = [(name, run.emitted(name).summary["gate_names_seen"]) for name in EMITS]
    sources += [(label(k), run.addend(k)["gate_names_seen"]) for k in run.addends]
    for name, seen in sources:
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
    _, hier = run.secp_hier()
    problems += [f"{label(k)} expanded upload has ops outside X, CX, CCX, SWAP" for k, r in hier.items() if sum(r["counts"].values()) != r["ops"]]
    if not run.quick:
        proc, row = run.secp_flat()
        if row is None or sum(row["counts"].values()) != row["ops"]:
            problems.append("256-bit flat build has ops outside X, CX, CCX, SWAP" if row else f"no 256-bit flat counts, see {proc.log}")
        evidence["secp256k1_flat_counts"] = row and row["counts"]
    evidence["problems"] = problems
    where = f"emitted toy and 256-bit JSON, {len(run.addends)} built addends, every toy instance, {len(hier)} expanded 256-bit uploads"
    return Result(FAIL if problems else PASS, "; ".join(problems) or f"only X, CX, CCX, SWAP, no measurement or reset, in {where}"
                  + ("" if run.quick else ", and the 256-bit flat build"), evidence)


def i5_addend_independence(run):
    profiles = {"G via scripts/emit.py": run.emitted("secp256k1_G.hier").summary}
    profiles.update((label(k), run.addend(k)) for k in run.addends)
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
    rows = {label(k): lean(run.addend(k)) for k in run.addends}
    ts = {name: r["t"] for name, r in rows.items()}
    worst, best = max(ts, key=ts.get), min(ts, key=ts.get)
    base = REFERENCE["max_t_over_1G_2G_3G"]
    problems = []
    g_file = run.emitted("secp256k1_G.hier").summary
    if g_file["sha256"] != run.addend(1)["sha256"]:
        problems.append("in-process [1]G JSON differs from scripts/emit.py's G")
    _, hier = run.secp_hier()
    cross = {label(k): {"counted": run.addend(k)["gates"], "expanded": r["counts"]} for k, r in hier.items()}
    problems += [f"{label(k)}: counts differ from the expanded upload" for k, r in hier.items() if run.addend(k)["gates"] != r["counts"]]
    problems += [f"{label(k)}: counted JSON is not the simulated JSON" for k, r in hier.items() if run.addend(k)["sha256"] != r["json_sha256"]]
    toy_hier, toy_flat = run.emitted("toy127.hier").summary, run.emitted("toy127").meta
    cross["toy127"] = {"hier_expanded": [toy_hier["ops"], toy_hier["ccx"]], "flat": [toy_flat["ops"], toy_flat["ccx"]]}
    if cross["toy127"]["hier_expanded"] != cross["toy127"]["flat"]:
        problems.append("toy127 expanded hierarchical counts differ from the flat emission")
    if not run.quick:
        _, flat = run.secp_flat()
        cross["secp256k1_G_flat"] = {"hier_expanded": g_file["gates"], "flat": flat and flat["counts"]}
        if flat is None or g_file["gates"] != flat["counts"]:
            problems.append("256-bit G expanded hierarchical gate counts differ from the flat build")
    evidence = {"t_per_ccx": T_PER_CCX, "seed": run.seed, "random_addends": run.sample, "per_addend": rows,
                "max": {"addend": worst, **rows[worst]}, "min": {"addend": best, **rows[best]},
                "reference_max_t": base, "max_t_delta_vs_reference": ts[worst] - base, "count_cross_checks": cross, "problems": problems}
    if problems:
        return Result(FAIL, "; ".join(problems), evidence)
    g = rows["[1]G"]
    note = (f"G: {g['ccx']:,} CCX, T {g['t']:,}, {g['ops']:,} ops, {g['num_wires']} wires; over {len(rows)} addends "
            f"max T {ts[worst]:,} ({ts[worst] - base:+,} vs reference {base:,}), min T {ts[best]:,}; counts match {len(hier)} expanded uploads")
    return Result(PASS, note, evidence)


def i7_limits(run):
    g = run.emitted("secp256k1_G.hier")
    size = g.path.stat().st_size
    per = {label(k): (run.addend(k)["bytes"], run.addend(k)["num_wires"]) for k in run.addends}
    largest = max(per, key=lambda a: per[a][0])
    widest = max(per, key=lambda a: per[a][1])
    wires = g.summary["num_wires"]
    problems = [f"{wires} wires > {MAX_WIRES}"] * (wires > MAX_WIRES) + [f"hierarchical JSON for G is {size:,} bytes >= {MAX_UPLOAD_BYTES:,}"] * (size >= MAX_UPLOAD_BYTES)
    problems += [f"{a}: {w} wires > {MAX_WIRES}" for a, (_, w) in per.items() if w > MAX_WIRES]
    problems += [f"hierarchical JSON for {a} is {b:,} bytes" for a, (b, _) in per.items() if b >= MAX_UPLOAD_BYTES]
    evidence = {"num_wires": wires, "max_wires": MAX_WIRES, "hier_json_bytes_G": size, "hier_json_file": rel(g.path),
                "largest_hier_json_over_addends": {"addend": largest, "bytes": per[largest][0]},
                "most_wires_over_addends": {"addend": widest, "wires": per[widest][1]}, "max_upload_bytes": MAX_UPLOAD_BYTES, "problems": problems}
    note = (f"{wires} of {MAX_WIRES} wires, at most {per[widest][1]} over {len(per)} addends; hierarchical JSON {size:,} bytes for G, "
            f"at most {per[largest][0]:,} over addends, cap {MAX_UPLOAD_BYTES:,}")
    return Result(FAIL if problems else PASS, "; ".join(problems) or note, evidence)


def i7_charge(run):
    return Result(PENDING, "needs the published benchmark: its charge formula is not public yet",
                  {"reason": "needs the published benchmark", "limit_as_reported": 33_610_176, "this_repo_as_reported": 19_492_771})


def r1_hygiene(run):
    run.wait_all()
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
    jobs: tuple = ()
    full_only: bool = False


EMIT_JOBS = tuple(f"emit:{name}" for name in EMITS) + ("addends",)
CHECKS = (
    Check("I1", "interface: layout and wire count", i1_interface, EMIT_JOBS),
    Check("I2.model", "exactness: classical dialog model and iteration bound", i2_model, ("model",)),
    Check("I2.toy", "exactness: exhaustive gate-level check on toy curves", i2_toy, ("toy",)),
    Check("I2.toy-hier", "exactness: same check on the hierarchical upload", i2_toy_hier, ("toy-hier",)),
    Check("I2.components", "exactness: 256-bit components on boundary inputs", i2_components, ("components",)),
    Check("I2.mutation", "exactness: the checker rejects mutated circuits", i2_mutation, ("mutation",)),
    Check("I2.secp", "exactness: 256-bit hierarchical upload on exceptional and random points", i2_secp, ("secp-hier", "emit:secp256k1_G.hier")),
    Check("I2.secp-flat", "exactness: 256-bit flat circuit at G", i2_secp_flat, ("secp-flat",), full_only=True),
    Check("I3", "gate set: X, CX, CCX, SWAP only", i3_gate_set, ("toy", "secp-hier", "secp-flat") + EMIT_JOBS),
    Check("I5", "addend independence of definitions", i5_addend_independence, EMIT_JOBS),
    Check("I6", "cost at G and over sampled addends (report)", i6_cost, ("secp-hier", "secp-flat") + EMIT_JOBS),
    Check("I7", "limits: wires and hierarchical upload", i7_limits, EMIT_JOBS),
    Check("I7.charge", "limits: benchmark charge", i7_charge),
    Check("R1", "run hygiene: no file git tracks changed", r1_hygiene),
)


def _git(*args):
    p = subprocess.run(["git", *args], cwd=ROOT, capture_output=True)
    return p.stdout if p.returncode == 0 else None


def _git_status():
    p = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all"], cwd=ROOT, capture_output=True, text=True)
    return p.stdout if p.returncode == 0 else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="simulate the 256-bit upload at two addends only, and skip the flat 256-bit check")
    ap.add_argument("--seed", type=int, default=2026, help="seed for the random addends and the random points")
    ap.add_argument("--sample", type=int, help=f"random addends besides [1]G..[{G_MULTIPLES[-1]}]G (default 4 quick, 16 full)")
    ap.add_argument("--only", nargs="+", metavar="ID", choices=[c.id for c in CHECKS], help="run only these checks")
    ap.add_argument("--out", type=Path, default=ROOT / "verify-out")
    args = ap.parse_args()
    sample = args.sample if args.sample is not None else (4 if args.quick else 16)
    run = Run(args.out.resolve(), args.quick, args.seed, sample)
    shutil.rmtree(run.work, ignore_errors=True)
    (run.work / "emit").mkdir(parents=True)
    t_start = time.time()

    selected = [c for c in CHECKS if not (args.only and c.id not in args.only)]
    skipped = [c.id for c in selected if c.full_only and args.quick]
    selected = [c for c in selected if c.id not in skipped]
    for check in selected:
        for group in check.jobs:
            if group in run.specs:
                run.start(group)

    results = []
    for check in selected:
        t0 = time.time()
        try:
            r = check.run(run)
        except Exception as e:
            r = Result(FAIL, f"check crashed: {e!r}"[:300], {"traceback": traceback.format_exc()})
        secs = round(time.time() - t0, 1)
        results.append({"id": check.id, "title": check.title, "status": r.status, "seconds": secs, "note": r.note, "evidence": r.evidence})
        print(f"{check.id:<14}{r.status.upper():<9}{secs:>6.1f}s  {check.title}\n{'':<31}{r.note}", flush=True)
    run.wait_all()

    tally = Counter(r["status"] for r in results)
    failed = [r["id"] for r in results if r["status"] == FAIL]
    pending = [r["id"] for r in results if r["status"] == PENDING]
    report = {"gate": FAIL if failed else PASS, "mode": "quick" if args.quick else "full", "seed": args.seed, "sample": sample,
              "addends": [label(k) for k in run.addends], "simulated_addends": [label(k) for k in run.simulated],
              "git_head": (_git("rev-parse", "HEAD") or b"").decode().strip() or None, "python": sys.version.split()[0],
              "seconds": round(time.time() - t_start, 1), "counts": dict(tally), "failed": failed, "pending": pending,
              "not_run_in_quick_mode": skipped, "results": results}
    (run.out / "accept.json").write_text(json.dumps(report, indent=1))
    print(f"\ngate {report['gate'].upper()}: {tally[PASS]} pass, {tally[FAIL]} fail, {tally[PENDING]} pending{f' {pending}' if pending else ''}"
          f"{f'; not run in --quick: {skipped}' if skipped else ''}; {report['seconds']}s; {rel(run.out / 'accept.json')}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
