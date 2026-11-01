# Reference circuit verification map

This directory is the maintained source for verifying what a user of the reference circuit can do. Read this index before driving anything, then use the matching feature file as the recipe.

## Baseline preconditions

- Run every command from the repo root.
- `python3` is 3.9 or newer. `.claude/skills/verify-reference-circuit/scripts/doctor.sh` exits 0.
- No other `scripts/accept.py` runs in this checkout. Each run deletes and rewrites `verify-out/accept/`.
- The benchmark checkout is absent unless a feature says otherwise. Only `benchmark-harness.md` needs it.

## Driving conventions

- Treat every command as literal. Keep flags, curve names and addend indexes unchanged.
- Write scratch output under `verify-out/`, which git ignores. Never write into `artifacts/` or `receipts/` during verification. Those hold committed evidence.
- Read results from the files a command writes, `verify-out/accept.json` and the `.meta.json` next to an emitted circuit, not only from the terminal.
- After each feature, run `scripts/evidence.sh <feature-id>` and then the cleanup in `SKILL.md`.

## Proof and skip reporting

- CLI proof includes the command, stdout, stderr, the exit code, and the files written.
- Record the feature id and the commit (`git rev-parse HEAD`) with every artifact. `evidence.sh` does both.
- A pending check is reported as pending, never as a pass.
- Report a feature that cannot run with the command you tried and the precondition that was missing.
- Don't report a feature as verified through a different feature's path.

## Feature entry contract

Each feature file starts with an H1 title and one paragraph on the user-visible behavior. It then has exactly four H2 sections in this order.

1. `Sub-features` lists short ids with one line for each behavior.
2. `How to get to it (user POV)` lists every user entry point.
3. `Driving it with the shell` starts with `Preconditions:` and pairs each user action with an exact command and an observable result.
4. `Gotchas` lists traps that can waste or invalidate a verification run.

## Features

- [Emit a circuit](./emit-circuit.md) covers flat and hierarchical emission at toy and 256-bit sizes, and the byte-identity check against the committed toy artifact.
- [Exactness checks](./exactness-checks.md) covers the model, toy, hierarchical toy, component, mutation and 256-bit sampled checks.
- [Per-addend cost report](./per-addend-cost.md) covers CCX, T, ops and wires at G and over seeded random addends, and the definition digests across addends.
- [Limits check](./limits.md) covers the interface layout, the wire cap, the upload cap, the gate set, and the pending charge.
- [Benchmark harness](./benchmark-harness.md) covers replay on the benchmark's own executor. It cannot run until the benchmark is published.
