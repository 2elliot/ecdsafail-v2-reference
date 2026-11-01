# Exactness checks

The exactness checks show that the circuit maps every promise input |c, R, 0> to |c, R + c*A, 0> with all scratch back at zero. `scripts/accept.py` runs the existing checkers and parses their output, under the ids `I2.model`, `I2.toy`, `I2.toy-hier`, `I2.components`, `I2.mutation` and `I2.secp`.

## Sub-features

- `exact-model` checks the classical dialog model exhaustively at 13 small primes and the iteration bound of 512 at secp256k1.
- `exact-toy` simulates every point of 36 toy curves with c = 0 and c = 1.
- `exact-toy-hier` runs the same simulation on the hierarchical emission after expansion.
- `exact-components` checks the 256-bit adders, doubler, negation and multiplier on boundary inputs.
- `exact-mutation` confirms the checker rejects CCX-to-CX demotions.
- `exact-secp` checks the 256-bit circuit at A = G on 120 random points plus infinity, A, -A, -2A, 2A and 3A. It runs only in full mode.

## How to get to it (user POV)

- Run `python3 scripts/accept.py --only I2.model I2.toy I2.toy-hier I2.components I2.mutation I2.secp`.
- Run `python3 scripts/accept.py --quick`, which includes every exactness check except `I2.secp`.
- Run a checker directly, for example `python3 -I scripts/check.py --toy-bits 4 6 7 8 9 10 --toy-per-size 6 --out verify-out/toy.json`.

## Driving it with the shell

Preconditions:

- `doctor.sh` exits 0.
- No other `accept.py` runs in this checkout.

- **Run the checks.** Run `python3 scripts/accept.py --only I2.model I2.toy I2.toy-hier I2.components I2.mutation I2.secp > verify-out/gate.log 2>&1; echo "exit $?" >> verify-out/gate.log`. It takes about 30 s. Its summary line starts with `gate PASS: 6 pass, 0 fail, 0 pending` and its last line is `exit 0`.
- **Read the rows.** Run `python3 -c "import json; [print(r['id'], r['status'], r['note']) for r in json.load(open('verify-out/accept.json'))['results']]"`. `I2.toy` and `I2.toy-hier` each report `36 toy instances, 25800 lanes, 0 failures`. `I2.components` reports `49 256-bit components` with `failing: none`. `I2.mutation` reports `163 of 200`. `I2.secp` reports `252 lanes` and `0 failures`.
- **Read the receipts.** Open `verify-out/accept/toy_exhaustive.json`. Every entry has `"failures": 0` and `"mode": "EXHAUSTIVE_PROMISE_DOMAIN"`. Open `verify-out/accept/secp256k1_sampled.json`. Its `counts` add up to its `ops`.
- **Proof.** Run `.claude/skills/verify-reference-circuit/scripts/evidence.sh exactness verify-out/gate.log`. The evidence directory holds `accept.json`, `toy_exhaustive.json`, `toy_exhaustive_hier.json`, `secp256k1_sampled.json` and every checker log.

## Gotchas

- `I2.mutation`'s count of 163 is for the reference circuit. A legitimate change can move it. The check fails only when the checker catches none, or when the unmutated circuit fails.
- A pass of `I2.secp` alone proves little. Random 256-bit points almost never reach the rare carries that separate an exact circuit from an approximate one.
- Exhaustive toy checks cover 256 bits only because the same code emits both. A code path that runs only at n = 256 passes every toy check. Read `INVARIANTS.md`, "What no script can check".
- The checkers write only where `accept.py` points them, under `verify-out/accept/`. Running `scripts/run_all.sh` instead rewrites the committed `receipts/`, which is not a verification step.
- `check.py` crashes rather than reporting failures when it meets a gate kind it does not know. `accept.py` reports that as `check.py exited 1`. Read the log it names.
