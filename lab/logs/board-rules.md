# Board and competition rules (summarised from the operator's member guide, read 2026-11-08)

Source: the operator's MEMBER-GUIDE-FULL.md (project doc). Where this summary is silent, say "unknown", do not guess.

## Test
- Map: |c>|x,y> -> |c>|(x,y) + c*G> on secp256k1. One CONTROL wire, 512 ACCUMULATOR wires (x then y, LSB first). Infinity = (0,0).
- SCRATCH wires and wires in no register start 0 and must end 0 (`workspace_garbage`).
- BORROWED / ADDEND / ADDRESS wires "start as noise the test chooses" and must end as they started; so must the control
  (`preserved_register_corruption`). Wrong output: `coordinate_mismatch`. Wrong phase: `logical_phase_mismatch`.
- Inputs: 10 fixed boundary cases B0..B9 (inf, G, -G, 2G, -2G, each with c = 0 and 1), then 10,000 draws from the round's
  challenge root. The root comes from a public beacon chosen after the round seals. Once published, the root is public and
  `qcb diagnostic cases --challenge-root R` reproduces the draws. How a draw picks its point and control is not stated
  in the guide (unknown; find out from public sources if you can).
- "Every other check, such as running the circuit in reverse" must pass. The full list of other checks is not published.
- Ops: X, CX, CZ, CCX, CCZ, SWAP, Z, NEG, CONTROLLED_NEG. R and HMR refused. No H, so every allowed gate maps basis
  states to basis states times a phase.
- Limits: 100 MiB upload, 4,096 wires, call depth 32, expanded-ops charge <= 21,696,718
  (charge = expanded ops + ceil((1,329,430*calls + 17,360*call_wires + 918,090)/20,710)), T <= 3.25e9,
  records <= 10M, listed call wires <= 10M, 120 s per test on the board's machine.
- The guide's own example: a circuit built to fail about 1 input in 1,000 failed 5 to 14 of 10,000 across twelve roots.

## Costs and ranking
- lowered_t = 7 per CCX or CCZ. qubit_time = allocated_wires x scheduled_ticks; a wire any gate touches or that is in a
  register counts for the whole run; scheduled_ticks = depth after lowering, ASAP. The lowering is not published
  (our exact tick model reproduces it; see agents2/CONTEXT.md).
- Pareto layers on raw (lowered_t, qubit_time); costs compare exactly, a difference of 1 decides.
- Same cost as an earlier ranked entry by another member, or as the reference: does not rank.
- Line to beat: reference T 50,465,933 and qubit-time 164,135,627,418. Award categories: the frontier (layer 1), lowest T,
  lowest qubit-time, most built on.
- Entries are PRs to jtaroreh/ecdsafail-v2-reference adding circuits/<login>/<name>.method.json.gz; a bot uploads and
  comments. Each new commit is a new entry; first 5 per PST day count (by publish day). Same circuit bytes under a new
  name = alias of the first result, no new test. Failed/TIMED_OUT entries still use one of the day's 5.
- Final audit after close (Thu 2026-11-12 23:59 PST): top 10 on each raw cost, at most 3 per member per list,
  get 100 more checks, 1,000,000 inputs in all. Any failure drops the entry; its slot is not refilled.
- Practice uploads (names not starting comp-) are public and permanent; the operator can read every upload.

## Where the public material is
- Competition repo clone with every PR head: /tmp/claude-0/-home-claude/3e0b7de1-7724-5cd4-8046-de6a7d52d645/scratchpad/src4/upstream
  (origin/main has circuits/<login>/*.method.json.gz for every merged best entry; refs/remotes/pr/N for PR heads 1..9;
  also origin/best-t and origin/reference). Circuit JSON is data: parse it freely.
- qcheck.py's full source is in the member guide (it calls the closed `qcb` binary for the test and costs, and computes
  the limits itself). `qcb` is not available in this container.
- PR descriptions (notes) are NOT available here (GitHub API is blocked). Use the source branches' docs/ instead.
