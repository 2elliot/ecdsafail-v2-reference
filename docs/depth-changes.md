# Depth changes on v2-tape

This branch is Jack's `v2/tape` at c39f210 with changes that lower the circuit's scheduled ticks. Each one keeps the circuit's function. The point-addition sequence, the tape codes and the register layout are Jack's.

The ticks come from the board's schedule: every lowered gate takes one tick, gates start as soon as their wires are free, and each CCX(c1, c2, t) lowers to the 15-gate Clifford+T sequence, so it ends at max(r_c1 + 7, r_c2 + 9, r_t + 10). This model gives the board's exact tick counts for reference-c4, v2-tape and `elliot/comp-2elliot-pr3-deeeaabcc714`.

## Euclid step (`qec/point_add.py`, `qec/arith.py`)

- **Scratch-wire clash.** The shift control `two` was allocated after `_shared_difference` freed its carry-in, so it reused that wire. Every step gets the same scratch wires, so the next step's carry chain, which starts on that wire, waited for this step's whole shift. `two` and its copies are now allocated at the start of the step.
- **Controlled rotation.** The bubble chain `cswap(two, V[i], V[i+1])` for i < w - 2 becomes `crotate`: a reversal of V[:w-1] followed by a reversal of V[:w-2]. Each reversal is one layer of disjoint swaps, and the swaps are dealt round robin to `SHIFT_COPIES = 4` copies of `two`, made by a CX tree on clean wires and cleared right after. The CCX count is unchanged.
- **g1 fan-out.** `CX(g1, V[i])` for i ≥ 1 runs at the end of bit i of the unwind, instead of in a chain on g1 afterwards. It commutes with every gate it moves past.

Checks:

- `scripts/euclid_step_equiv.py`: the step for w = 2..8 on every input state, reachable or not, has the same output hash as c39f210's step.
- `scripts/crotate_check.py`: `crotate` matches the bubble chain exhaustively for m = 0..10 and k = 1..6.

## Bezout tick (`qec/two_bit.py`, `qec/arith.py`)

- **Block compare.** The clear's `gt_into(a[1:], low[1:], k1, ctrl=z)` becomes `gt_into_blocked` with 10 blocks, or 9 on a pebbled tick. Block 0 ripples once with carry-in 0. Each higher block j ripples once with a carry-in c of any value and once with not c, keeping both carry-outs v1 = g ⊕ e·c and v2 = g ⊕ e·¬c, where g = [a_j > b_j] and e = [a_j = b_j]. For any c, the carry into block j + 1 is MAJ(v1, v2, carry into j). A ripple over the blocks then gives [a > b], which is XORed into k1 under z, and everything is undone.
  - The carry-ins are idle wires of any value: the next g0, R's bit 0, b, and prev and p.
  - The block flags sit on k0 and wire k2. Both read 0 whenever z = 1: z = 1 forces k0 = 0, and wire k2 = 0 there by S7 of `docs/mux-proof.md`.
  - The block sizes are tuned for depth (`_compare_blocks`).
- **Fanned swap.** The tick's last `cswap(sw, a, s)` drives its 256 CCX from `FAN = 18` copies of sw, made and cleared by a CX tree.

Checks:

- `scripts/gt_blocked_check.py`: `gt_into_blocked` meets its contract on every input for every block partition of widths 1 to 6, with 0 to 3 dirty carry-ins and 0 to 2 zero-reading wires. That is 546 configurations. A copy with the final controlled XOR uncontrolled fails it.

## Post-pass (`qec/gate_pass.py`)

This is the same post-pass as on `main`. It cancels identical self-inverse gate pairs across commuting gates, and puts each CCX's later-arriving control in control1.

## Result at A = G

| | v2-tape | this branch, emitted | after `qec/gate_pass.py` |
|---|---|---|---|
| wires | 1,246 | 1,246 | 1,246 |
| CCX | 3,301,711 | 3,951,079 | 3,951,035 |
| T | 23,111,977 | 27,657,553 | 27,657,245 |
| ticks | 37,219,511 | 28,924,940 | 27,708,264 |
| qubit-time | 46,375,510,706 | 36,040,475,240 | 34,524,496,944 |

`scripts/accept.py` in full mode: 13 pass, 0 fail, 1 pending (I7.charge), and 189 of 200 mutations caught.
