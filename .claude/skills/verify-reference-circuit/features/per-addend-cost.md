# Per-addend cost report

The cost report gives CCX, T = 7 x CCX, ops and wires for the 256-bit circuit at G, at [2]G to [8]G and at a seeded sample of random addends, then the sample maximum and the sample minimum. Those two are statistics of the sample, not bounds over every addend. It also checks that every definition in the hierarchical circuit is identical across those addends, comparing the sha256 of each definition after parsing and compact re-serialization, so the addend reaches only the top level. `scripts/accept.py` runs these as `I6` and `I5`.

## Sub-features

- `cost-g` reports the cost at A = G.
- `cost-sample` reports the cost over [1]G to [8]G and the seeded random addends.
- `cost-max` reports the sample maximum T and its delta against the reference, 50,465,989.
- `cost-crosscheck` confirms the counts equal the expanded upload that `I2.secp` simulated, at every simulated addend, and equal the flat circuit at toy size and, in full mode, at G.
- `defs-identical` confirms all definition digests are equal across addends.

## How to get to it (user POV)

- Run `python3 scripts/accept.py --only I5 I6`.
- Change the sample with `--sample N` and `--seed S`. [1]G to [8]G always run. The default adds 1 random addend with `--quick` and 16 without.

## Driving it with the shell

Preconditions:

- `doctor.sh` exits 0.
- No other `accept.py` runs in this checkout.

- **Run the report.** Run `mkdir -p verify-out; python3 scripts/accept.py --quick --only I5 I6 > verify-out/gate.log 2>&1; echo "exit $?" >> verify-out/gate.log`. It takes about 15 s, because I6 waits for `I2.secp`'s expansion to check its counts. The log shows `I5 PASS` with `all 261 definitions identical over 10 builds` and `I6 PASS` with `G: 7,209,419 CCX, T 50,465,933, 18,995,461 ops, 2311 wires`.
- **Read the maximum.** Run `python3 -c "import json; e=[r for r in json.load(open('verify-out/accept.json'))['results'] if r['id']=='I6'][0]['evidence']; print(e['sample_max']['addend'], e['sample_max']['t'], e['sample_max_t_delta_vs_reference'])"`. At the reference it prints `[2]G 50465989 0`.
- **Read the per-addend rows.** In the same evidence, `per_addend` maps each addend label to `ccx`, `t`, `ops`, `num_wires`, `top_ccx` and `bytes`. At the reference, `[2]G` and `[3]G` have `ccx` 7,209,427 and `top_ccx` 18,131, and `[1]G` has 7,209,419 and 18,123.
- **Read the definition digests.** In the `I5` evidence, `differences` is an empty list and `definitions` is 261. `top_level_ccx` lists the top-level CCX per addend. Only these differ.
- **Proof.** Run `.claude/skills/verify-reference-circuit/scripts/evidence.sh per-addend-cost verify-out/gate.log`.

## Gotchas

- I6 never fails on cost. It is a report. It fails only when its counts disagree with the expanded upload or the flat circuit.
- Quote the sample maximum with its sample, not the cost at G. G is not the worst addend.
- The sample maximum does not bound the cost at addends outside the sample, and the true maximum over every A is unproven. A different `--seed` can find a different maximum. At the reference, [2]G, [3]G and [8]G reach the largest T seen over 1,024 addends. A change that adds cost only at an unsampled addend passes I6 unseen. Keep a change's cost independent of A, or argue its maximum. See `INVARIANTS.md`, "Report cost over a sample of addends".
- `I5` fails with the first addend and definition index that differ. An addend-dependent constant inside a recorded block is the usual cause.
- Each 256-bit addend costs about 1.3 s to build, emit and summarize, and about 10 s more to expand and simulate in full mode. Six simulations run at a time. `--sample 64` without `--quick` takes several minutes.
