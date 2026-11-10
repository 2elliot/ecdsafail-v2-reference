# Lab brief: scouting new ideas for the qubit-time track

Shared brief for each scout task. Read it in full before you start. Your own task prompt names your idea, your ID and anything it changes here.

## The goal
QCB runs a friendly circuit competition. The board is jtaroreh/ecdsafail-v2-reference, scored at ecdsabreak.taroreh.com. Elliot (GitHub 2elliot) holds the best qubit-time entry, `cold-trees` (PR #12): 1,250 wires × 5,349,345 ticks = 6,686,681,250, with T 49,799,519.

Small optimizations of this tree are easy and are not what we want now. We want breakthroughs: new constructions that could cut qubit-time by a large factor. Take risks. A clear "this is dead, and here's why" with numbers is a useful result too.

This phase is scouting. Prototype the core of your idea in isolation, measure it with the exact tick model, project it to the full circuit, and report. Don't integrate it into the whole generator unless your prompt says so.

## The task the circuit must do
- **The map.** An exact controlled point addition |c⟩|x, y⟩ → |c⟩|(x, y) + c·G⟩ on secp256k1, with G the generator. One control wire and 512 accumulator wires (x then y, LSB first). Infinity is (0, 0).
- **Exactness.** It must be right on every input, including infinity, ±G and ±2G. Scratch starts and ends at 0. No measurement or reset (R and HMR are refused), so every uncompute is gates.
- **Gates.** X, CX, CCX, SWAP, plus CZ, CCZ, Z, NEG and CONTROLLED_NEG. The board's tick costs for the last five are unverified, so avoid them.
- **Costs.**
  - T = 7 per CCX. It must stay at or below 50,465,933, which is CCX ≤ 7,209,419.
  - Qubit-time = allocated wires × scheduled ticks, ranked on the qubit-time track. Every touched wire counts for the whole run.
- **Limits.**
  - Wires ≤ 4,096.
  - Expanded-ops charge ≤ 21,696,718, where charge = expanded ops + ⌈(1,329,430·calls + 17,360·call_wires + 918,090) / 20,710⌉. Use `lab/tools/charge.py`.
  - Upload ≤ 100 MiB.
  - Definitions may only call earlier ones, with call depth ≤ 32.

## The tick model (exact; it matches the board on every published circuit)
- **Scheduling.** ASAP. Every lowered gate takes 1 tick. X, CX and Z take 1 tick, and SWAP takes 3.
- **CCX(c1, c2, t)** lowers to H t; CX c2→t; T† t; CX c1→t; T t; CX c2→t; T† t; CX c1→t; T c2; T t; H t; CX c1→c2; T c1; T† c2; CX c1→c2.
  - It finishes at about max(r_c1 + 7, r_c2 + 9, r_t + 10), where r_w is when wire w is free.
  - It holds its controls for about 15 ticks, and frees its target after about 11.
- **Consequences** (all measured):
  - A ripple whose carry enters each CCX as a control costs about 7 ticks a bit going up. The undo costs about 11 a bit coming down.
  - Many CCX that share one control serialize on it, about 11 ticks each.
  - The builder hands out the most recently freed clean wire. So a stage can wait on a wire that something slow has only just released. Round 1's biggest gain came from choosing the clean wire that is free first.
- `lab/tools/ticks.py <circuit-or-method.json>` runs the model on a whole file. `lab/tools/tickmodel.py` runs it on a builder op list.

## The current design (tree at this branch's parent, depth/v5 @ 435d287)
The design is Jack's dialog-Euclid generator with Mike's PR8 changes and ours. `DESIGN.md` has one section per change; read the sections from "Depth pass" to the end.

The point addition calls a multiplier twice (M⁻¹ to get λ, then M to uncompute it). Between the calls it squares λ. Each multiplier call does three passes:
- **E, the Euclid pass.** It runs the binary-GCD dialog and writes a compact tape: a two-bit symbol per tick, with a held offset and a sweep.
- **B, the Bezout pass.** It applies the tape to the coefficients.
- **E⁻¹.** It cleans the tape.

Each pass runs tick_bound(q) = 384 ticks at secp256k1. 385 is proven; 384 rests on an unproved family lemma. Random inputs need about 302 ticks.

Where the cost is:
- **The Bezout tick** (`TwoBitMod.tick` in `qec/two_bit.py`) was about 75% of the critical path, measured before round 2. It ran about 6,000 ticks per tick, as 768 ticks in series over the two calls. Inside it are the add-or-shift (`add_or_shift_la`, about 40% of the tick), the compare (`compare_blocks`), the prefix P (`PrefixR`, 221 bits) and two carry chains. All of them are full-width carry propagations.
- **The Euclid passes** are about 14%, and the square about 7%.
- **Trees can be much faster than ripples.** Measured on round 2's arm D: one 256-bit constant add plus compare went from about 18,000 ticks as ripples to about 1,110 as trees, at 2.6× the CCX and about 513 scratch wires.
- **The binding limit is T.** The current circuit has 7,114,217 CCX, which leaves about 95,000 of headroom. Jack's `v2-lean` (PR #11) runs the same algorithm at 3,278,326 CCX and 19,861,888 ticks on 1,246 wires.
- **Logs.** `lab/logs/*-decision.tsv` is every attempt from the last two nights, with numbers. Read the one for your area so you don't repeat a dead end.

## Tools (in lab/tools; run Python as `python3 -I`, and the generator as `python3.12 -I -B`)
- **`sh lab/tools/measure.sh <tree> <tag> [MODULE.FLAG=value ...]`** is the frozen full-circuit cost harness.
  - It takes about 5 minutes and 1 GB, and is deterministic.
  - On this branch it prints `wires 1250 CCX 7114217 ... ticks 5349345 qubit_time 6686681250`.
  - `KEEP=<path>` keeps the output.
- **`emitvar.py`** emits under flags. **`buildonly.py`** gives CCX and wires without emitting, and is fast.
- **Screening and profiling:**
  - `bzpass.py` screens the Bezout pass in about 11 s and matches the harness's raw ticks.
  - `cpath_tick.py` charges the tick's critical path to source lines.
  - `cpath_json.py` with `phasemap.py` charges a whole JSON circuit's path to definitions and phases.
  - `prof.py` makes one part free. Its output is WRONG and only an upper bound.
  - `iso.py` times a component in isolation.
- **Exactness:**
  - `qlocal.py <file> --draws 10000 --seed 1 --reverse` is our bit-sliced stand-in for the board's test, and exits 0 on a pass.
  - `qhard.py <file> --cases lab/data/hard_cases_ext.json --draws 2000 --seed 2 --reverse` adds 2,313 hard inputs.
  - `topeq.py a.json b.json lab/tools [seed]` compares two circuits' maps on 512 random full states.
  - The tree's `scripts/accept.py --quick` is the generator's own gate. It takes about 15 minutes, needs a clean git tree, and its I2.pass check may be killed for memory (rerun it alone with `--only I2.pass`).
  - The tree's `scripts/*_check.py` are exhaustive component checks. Each has a `broken` mode that must fail.
- **`wrap.py`** is the board's upload wrapper, copied from its release qcb-8376098.

## Rules
- Exact circuits only. Every construction you propose must be exact on every input. Say how you would prove it, and test it exhaustively at small sizes. Never build or measure a circuit that is inexact on purpose.
- Local only. Never upload to the board, open pull requests, push, or post anywhere.
- One heavy process at a time (emit, harness, qlocal, accept). This machine has 2 cores and about 7 GB.
- Label every number: measured (you ran it), computed (arithmetic from measured numbers), source-read, inferred, or guess.
- Your report decides whether we spend days on this idea, so be honest about what blocks it.

## Deliverable
A report of at most 1,500 words, with these parts:
- **Verdict.** Promising, maybe, or dead, with one sentence why.
- **Prototype.** What you built and measured, with numbers.
- **Projection.** Projected full-circuit wires, ticks, CCX and qubit-time, with the method and its error bars.
- **Blockers and next steps.** What blocks it, and the next three steps if we go ahead.
- **Code.** Where your code is.

Write the report with the Projects tool to `claude/lab/<ID>-report.md`. Put your prototype code, as fenced blocks, in `claude/lab/<ID>-code.md`. If the Projects tool isn't available, publish both as one private Artifact titled "Lab <ID> report". Also give the report as your final message.
