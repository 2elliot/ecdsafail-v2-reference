# Per-addend cost report

The cost report gives CCX, T = 7 x CCX, ops and wires for the 256-bit circuit at G, [2]G, [3]G and a seeded sample of random addends, then the maximum and the minimum. It also checks that every definition in the hierarchical circuit is byte-identical across those addends, so the addend reaches only the top level. `scripts/accept.py` runs these as `I6` and `I5`.

## Sub-features

- `cost-g` reports the cost at A = G.
- `cost-sample` reports the cost over [1]G, [2]G, [3]G and the seeded random addends.
- `cost-max` reports the maximum T and its delta against the reference, 50,465,989.
- `cost-crosscheck` confirms the expanded hierarchical counts equal the flat circuit, at toy size always and at 256 bits in full mode.
- `defs-identical` confirms all definition digests are equal across addends.

## How to get to it (user POV)

- Run `python3 scripts/accept.py --only I5 I6`.
- Change the sample with `--sample N` and `--seed S`. The default sample is 4 random addends with `--quick` and 16 without.

## Driving it with the shell

Preconditions:

- `doctor.sh` exits 0.
- No other `accept.py` runs in this checkout.

- **Run the report.** Run `python3 scripts/accept.py --quick --only I5 I6 > verify-out/gate.log 2>&1; echo "exit $?" >> verify-out/gate.log`. It takes about 11 s. The log shows `I5 PASS` with `all 261 definitions identical over 4 builds` and `I6 PASS` with `G: 7,209,419 CCX, T 50,465,933, 18,995,461 ops, 2311 wires`.
- **Read the maximum.** Run `python3 -c "import json; e=[r for r in json.load(open('verify-out/accept.json'))['results'] if r['id']=='I6'][0]['evidence']; print(e['max']['addend'], e['max']['t'], e['max_t_delta_vs_reference'])"`. At the reference it prints `[2]G 50465989 0`.
- **Read the per-addend rows.** In the same evidence, `per_addend` maps each addend label to `ccx`, `t`, `ops`, `num_wires`, `top_ccx` and `bytes`. At the reference, `[2]G` and `[3]G` have `ccx` 7,209,427 and `top_ccx` 18,131, and `[1]G` has 7,209,419 and 18,123.
- **Read the definition digests.** In the `I5` evidence, `differences` is an empty list and `definitions` is 261. `top_level_ccx` lists the top-level CCX per addend. Only these differ.
- **Proof.** Run `.claude/skills/verify-reference-circuit/scripts/evidence.sh per-addend-cost verify-out/gate.log`.

## Gotchas

- I6 never fails on cost. It is a report. It fails only when its own counting disagrees with the flat circuit.
- Quote the maximum over addends, not the cost at G. G is not the worst addend.
- The maximum is over the sample you ran. A different `--seed` can find a different minimum. [2]G and [3]G always run, and at the reference they reach the maximum seen over 1,024 addends.
- `I5` fails with the first addend and definition index that differ. An addend-dependent constant inside a recorded block is the usual cause.
- Each 256-bit addend costs about 3 s to build, emit and parse. `--sample 64` takes several minutes.
