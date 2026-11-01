# Emit a circuit

`scripts/emit.py` writes the controlled point-addition circuit as v2-ir `ArithmeticCircuit` JSON, flat or hierarchical, for a toy curve or for secp256k1 with A = G. It writes a `.meta.json` beside each circuit and prints the same metadata as one JSON line.

## Sub-features

- `emit-toy-flat` writes the flat toy circuit for q = 127 and the point at index 5.
- `emit-toy-hier` writes the same toy circuit as definitions and calls.
- `emit-secp-hier` writes the 256-bit circuit for A = G as definitions and calls.
- `emit-secp-flat` writes the flat 256-bit circuit, gzipped.
- `emit-identity` confirms the toy flat emission is byte-identical to the committed `artifacts/toy127_circuit.json`.

## How to get to it (user POV)

- Run `python3 scripts/emit.py --curve toy:<q> --addend-index <i> --out <file>` for a toy curve.
- Run `python3 scripts/emit.py --out <file>` for secp256k1 with A = G.
- Add `--hier` to either for the hierarchical form. End `--out` in `.gz` to gzip the output.

## Driving it with the shell

Preconditions:

- `doctor.sh` exits 0.
- `verify-out/emit/` does not exist or holds nothing you need.

- **Toy flat.** Run `python3 -I scripts/emit.py --curve toy:127 --addend-index 5 --out verify-out/emit/toy127_circuit.json`. Exit code 0. The printed line has `"num_wires": 70`, `"ops": 14803`, `"ccx": 5099` and `"json_bytes": 582590`. `verify-out/emit/toy127_circuit.meta.json` holds the same values and a `layout` of control `[0]`, x `1..7`, y `8..14`.
- **Byte identity.** Run `shasum -a 256 verify-out/emit/toy127_circuit.json artifacts/toy127_circuit.json`. Both digests are equal. A change under `qec/` that is meant to leave the circuit alone must keep them equal.
- **Toy hierarchical.** Run `python3 -I scripts/emit.py --curve toy:127 --addend-index 5 --hier --out verify-out/emit/toy127_circuit.hier.json`. The printed line has `"definitions": 11`, `"call_sites": 39`, `"primitive_ops": 1974` and `"json_bytes": 82002`.
- **256-bit hierarchical.** Run `python3 -I scripts/emit.py --hier --out verify-out/emit/secp256k1_G_circuit.hier.json`. It takes about 2 s. The printed line has `"num_wires": 2311`, `"definitions": 261`, `"call_sites": 1796`, `"primitive_ops": 972903` and `"json_bytes": 43422940`.
- **256-bit flat.** Run `python3 -I scripts/emit.py --out verify-out/emit/secp256k1_G_circuit.json.gz`. It takes about 40 s. The printed line has `"ops": 18995461`, `"ccx": 7209419` and `"num_wires": 2311`.
- **Proof.** Run `.claude/skills/verify-reference-circuit/scripts/evidence.sh emit-circuit verify-out/emit/toy127_circuit.meta.json verify-out/emit/secp256k1_G_circuit.hier.meta.json`. The evidence directory holds both meta files and `git.txt`.

## Gotchas

- `--addend-index` indexes `curve.points()`, where index 0 is infinity. Index 0 fails the generator's assertion, because A must be a finite point.
- `emit.py` supports only A = G at 256 bits. For other addends, `scripts/accept.py` builds `[k]G` in process with the same `write_circuit` call. See `per-addend-cost.md`.
- The flat 256-bit JSON is 822 MB uncompressed. Keep the `.gz` suffix, and remove it in cleanup.
- The meta file name drops everything from the first `.json`. `x.hier.json` gets `x.hier.meta.json`, and `x.json.gz` gets `x.meta.json`.
- The numbers above hold at commit `8ca9b0e` and every later commit that leaves the circuit alone. A change that lowers cost changes them on purpose. Then the byte-identity step is expected to differ, and the exactness checks carry the proof.
