# Invariants for changes to the circuit

This page is for you if you change `qec/` to make the circuit cheaper. It says what `scripts/accept.py` protects and why, what no script can check for you, and how to hand a change back.

The benchmark that will judge your change is not published yet. `accept.py` encodes the invariants that the circuit must keep under any scoring rule we expect. A change that passes it is likely to be accepted once the benchmark ships. Nothing here is a promise about the final scoring.

## Run the gate

```sh
python3 scripts/accept.py --quick   # 18 s, toy exactness and the cheap 256-bit structure checks
python3 scripts/accept.py           # 40 s, adds the 256-bit sampled check and 19 addends instead of 7
```

Both times are from an Apple M4 laptop with Python 3.14 under light load. The gate runs on one core, and under heavy load the same runs took 27 s and 79 s. Python is stdlib only. The gate prints one row per check and writes the full evidence to `verify-out/accept.json`. It exits 0 only when no check fails. A check marked pending waits on the published benchmark and fails nothing.

Run `--quick` while you iterate and the full mode before you hand a change back. `--only I5 I6` runs a subset. `--sample N` and `--seed S` change the random addends that I5, I6 and I7 build.

## What each check protects

The circuit implements |c>|R>|0> -> |c>|R + c*A>|0> for a fixed public point A on secp256k1, under the frozen profile `secp256k1-controlled-exact-v1`. Every check below guards one part of that contract or one limit of the benchmark.

**I1 interface.** The control is wire 0, x is wires 1..n, and y is wires n+1..2n, at n = 7 and at n = 256. Every other wire is scratch, and the wire count stays at most 4,096. The benchmark reads the declared registers from these positions, so a reordered layout computes the right sum in the wrong place. The reference uses 2,311 wires.

**I2 exactness.** The profile requires an exact unitary, so the circuit must be right on every input, including infinity and R in {A, -A, -2A}. I2 runs the existing checkers and parses what they print. It does not reimplement them.

- `I2.model` runs `scripts/model_check.py`. It checks the classical dialog model exhaustively at 13 small primes and checks the Euclid width schedule against the proven iteration bound.
- `I2.toy` runs `scripts/check.py` on 36 toy curves of 4 to 10 bits. It simulates every point with c = 0 and c = 1 at gate level and requires all scratch back at zero.
- `I2.toy-hier` runs the same check on the hierarchical emission after `qec/flatten.py` expands it. The upload is hierarchical, so a broken call map would pass `I2.toy` and still upload a wrong circuit.
- `I2.components` runs `scripts/component_check.py`. It checks the 256-bit adders, doubler, negation and multiplier on boundary inputs near q and 2^256 that random points almost never reach.
- `I2.mutation` runs `scripts/mutation_control.py`. It confirms that the checker still rejects mutated circuits. The reference catches 163 of 200 CCX-to-CX demotions. The rest are equivalent circuits.
- `I2.secp` runs only in full mode. It checks the 256-bit circuit at A = G on 120 random points plus the exceptional ones.

**I3 gate set.** Every emitted primitive is X, CX, CCX or SWAP, with no measurement and no reset. The profile disables measurement cleanup, and a circuit of these four gates is a permutation with no phases. That is why basis-state exactness with clean scratch gives the exact unitary on superpositions. A single Z, CZ or measurement breaks that argument even when every basis check passes. I3 reads the emitted JSON at toy and 256-bit sizes and counts the gate kinds in every toy instance.

**I5 addend independence of definitions.** The gate builds the 256-bit kernel hierarchically for [1]G to [8]G and the seeded random addends, 1 of them in quick mode and 16 in full mode. Every definition must be byte-identical across them. Today all 261 definitions are identical, and A enters only through the top level, in the `cmodadd_const` calls in `_kernel` and the four exception compares and loads in `build` (`qec/point_add.py:110-190`). Keep it that way. A whole attack runs this kernel for many addends, and shared definitions let one upload of the expensive parts serve every addend.

**I6 cost at G and sample statistics.** I6 reports CCX, T = 7 x CCX, ops and wires at G. Over the addends the gate built, [1]G to [8]G and a seeded random sample, it reports the sample maximum and the sample minimum. These are statistics of that sample, not bounds. I6 shows the sample maximum T against the reference, 50,465,989, the largest T any sampled addend reached at commit `8ca9b0e`. I6 never fails on cost. It fails only when its counts disagree with the expanded upload that `I2.secp` simulated, at every simulated addend, or with the flat circuit at toy size and, in full mode, at G.

**I7 limits.** Wires stay at most 4,096. The hierarchical JSON for G, which the gate writes to `verify-out/accept/emit/`, stays under 104,857,600 bytes, and so does the JSON for every sampled addend. The reference is 43,422,940 bytes.

**I7.charge.** The next benchmark release charges a hierarchical upload with a formula that is not public. This check stays pending until it is.

**R1 run hygiene.** The gate compares `git status` before and after it runs. Every checker writes under `verify-out/`, so a run that changes a tracked file is a bug in the gate.

## What no script can check

**One generator for every curve.** The exhaustive toy check says something about 256 bits only because the same code emits both sizes. Only the constants and the bound-derived width schedule differ. A code path that runs only at n = 256, or only for secp256k1's q, escapes every exhaustive check. One mutation used to test this gate added 2.8 million identity gates only when n = 256. It passed all of I2 and was caught only by the upload limit. Write size-generic code. If a trick needs a property of secp256k1's q, test that property on toy primes that share it, or prove it.

**The iteration bound is proven, not sampled.** The Euclid loop runs bitlen(q(q-1)) = 512 iterations at secp256k1 because u*v at least halves per step. Random inputs finish far earlier. The paper's approximate circuits use a shorter count and fail on rare inputs. Shrinking the bound needs a proof that covers every input, not a sample that passed.

**Random 256-bit sampling cannot tell exact from approximate.** Most CCX in the constant adders and the idle Euclid tails fire only on rare carries. A circuit that drops them passes random 256-bit lanes and is still wrong (README, exactness evidence item 5). Exactness rests on the exhaustive toy check, the 256-bit boundary components and the proven bound. A passing `I2.secp` alone proves little.

## Where the CCX are

At A = G the circuit has 7,209,419 CCX.

| Part | CCX | Share |
|---|---|---|
| Two dialog multipliers, 3,140,864 each | 6,281,728 | 87.1% |
| Lambda-squared subtract | 910,080 | 12.6% |
| Flags, exception loads, constant adds | 17,611 | 0.2% |

The multipliers and the lambda-squared subtract hold 99.8% of the CCX, and none of them reads A. Spend your effort there. That work also carries over to a future windowed kernel that takes the addend as a quantum input, because that kernel needs the same multiplier and subtract. Constant-adder tricks do not carry over, since a windowed kernel has no classical constant to exploit. They also act on the 0.2% that depends on A.

## Report cost over a sample of addends

Cost depends on A. [2]G and [3]G each cost 7,209,427 CCX, 8 more than G. Over 1,024 sampled addends at the reference, T ranged from 50,465,177 to 50,465,989. The benchmark may test any addend, and a whole attack runs hundreds of them.

The sample maximum that I6 prints is the largest cost among the addends the gate built. It does not bound the cost at other addends. The true maximum over every A is unproven. A change can raise the cost at an addend outside the sample, and the gate will not see it. So keep the cost of your change independent of A. A change inside the definitions does that, because I5 holds them identical across addends. If your change touches the top level, where A enters, argue its maximum over every A. A change that helps G and hurts another addend is a regression. Quote the I6 sample maximum with the sample it came from, not the cost at G.

## Don't trade ops for many small calls

The next benchmark release charges per call and per definition, and those terms may rise before it ships. Splitting the circuit into many small definitions or many calls can lower the op count and still raise the charge. Prefer changes that remove CCX and ops outright.

## Not settled yet

- The scoring rule and the caps of the next benchmark release.
- Which whole-attack route comes first. One candidate is the binary wrapper with 512 controlled additions of fixed points. The other is a windowed kernel with a coherent addend.

## Hand a change back

Open a pull request against this repository with:

1. The change, with nothing under `receipts/` edited by hand.
2. The output of `python3 scripts/accept.py` in full mode, pasted into the description, and its `verify-out/accept.json` attached.
3. The I6 sample maximum T, the sample it came from, and its delta against 50,465,989. For a change to the top level, add the argument for its maximum over every A.
4. For any change to the iteration bound, the width schedule or a modular reduction, the argument that it stays exact on every input.
