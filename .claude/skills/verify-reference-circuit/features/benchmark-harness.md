# Benchmark harness

Not runnable until the benchmark is published. The harness replays an emitted circuit on the benchmark's own code. It parses and validates the circuit with `v2_ir`, executes test lanes on `v2_runtime::BitSlicedExecutor`, round-trips a `SubmittedMethod`, evaluates `QualificationGraph`, and lowers to Clifford+T with `v2_lowering`. It needs a checkout of the ecdsafail-v2 benchmark at `../ecdsafail-v2`, next to this repo, and that repository is not public yet.

## Sub-features

- `harness-build` builds `harness/` offline against the benchmark checkout.
- `harness-toy` replays the toy q = 127 circuit on all 254 lanes and writes a `SubmittedMethod`.
- `harness-secp` replays the 256-bit circuit at A = G on 510 lanes and lowers it.
- `harness-run-all` runs both through `scripts/run_all.sh`.

## How to get to it (user POV)

- Run `ECDSAFAIL_V2=../ecdsafail-v2 scripts/run_all.sh`.
- Build and run the harness binary directly, as in the commands below.

## Driving it with the shell

Preconditions:

- The benchmark is published and checked out at `../ecdsafail-v2`. `doctor.sh` reports `benchmark checkout at ../ecdsafail-v2`.
- `cargo` is on `PATH` and the crates the harness needs are in the local cargo cache, because the build is offline.
- Today neither holds outside the maintainers' machines. Report this feature as `could not run: benchmark not published` and stop.

- **Build.** Run `(cd harness && cargo build --release --offline --locked)`. Exit code 0 and `harness/target/release/qec-harness` exists.
- **Toy replay.** Run `mkdir -p verify-out/emit`, then `python3 -I scripts/emit.py --curve toy:127 --addend-index 5 --out verify-out/emit/toy127_circuit.json`, then `python3 -I scripts/make_cases.py --curve toy:127 --addend-index 5 --out verify-out/emit/toy127_cases.json`, then `harness/target/release/qec-harness --circuit verify-out/emit/toy127_circuit.json --cases verify-out/emit/toy127_cases.json --lower > verify-out/harness_toy127.json`. The JSON has `execution.lanes` 254, `execution.failures` 0, `lanes_with_dirty_scratch` 0 and `lowering.measurements` 0.
- **256-bit replay.** Run `mkdir -p verify-out/emit`, then `python3 -I scripts/emit.py --out verify-out/emit/secp256k1_G_circuit.json.gz`, then `python3 -I scripts/make_cases.py --random 250 --out verify-out/emit/secp256k1_cases.json`, then `gzip -dc verify-out/emit/secp256k1_G_circuit.json.gz | harness/target/release/qec-harness --circuit - --cases verify-out/emit/secp256k1_cases.json --lower > verify-out/harness_secp256k1.json`. The JSON has `execution.lanes` 510, `execution.failures` 0 and `operations` 18,995,461. On the published benchmark the validator rejects it on the operation cap, which README "Unresolved gates" describes.
- **Proof.** Run `.claude/skills/verify-reference-circuit/scripts/evidence.sh harness verify-out/harness_toy127.json verify-out/harness_secp256k1.json`.

## Gotchas

- `harness/Cargo.toml` reads the benchmark at `../../ecdsafail-v2` relative to `harness/`. A checkout anywhere else needs a symlink at `../ecdsafail-v2`. `run_all.sh` refuses to build when `ECDSAFAIL_V2` points somewhere else.
- `scripts/run_all.sh` rewrites the committed `receipts/` and `artifacts/`. Use it to regenerate receipts, not to verify a change. The direct commands above write only under `verify-out/`.
- Every replay ends `is_qualified: false`. The benchmark holds no sound theorem checker yet, so no artifact can qualify. That is not a harness failure.
- The 256-bit replay peaks at about 2.5 GB of memory (`receipts/harness_secp256k1.time`). Stream the gzip into it. Don't decompress the flat circuit to disk.
