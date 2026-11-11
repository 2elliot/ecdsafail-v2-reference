# The rank-one start of a Bezout pass

The construction is in `qec/rank_one.py`, and `DialogMultiplier._rank_start` in `qec/point_add.py` runs it.
Each statement is numbered so that it can be tested as worded.

## What it does

The Bezout pass starts at (R, S) = (y, 0) and runs the swept tape's cells from the last to the first.
After its first j cells, (R, S) ≡ (U·y, V·y) mod q for the integer cofactor pair (U, V) of those cells, which starts at (1, 0).
`_rank_start` runs the first J cells on integers with no reduction.
Cell T − 1 is the terminal S ^= g0·R, as before.
Cells T − 2 to T − J + 1 are two_bit's add-or-shift and swap on words of n + w + 3 wires, with the held offset after a B1 as in `docs/mux-proof.md`.
Cell T − J is semantic: a B1 doubles S and adds nothing, so no lane is held after it.
One conversion then leaves (U·y mod q, V·y mod q) and every other wire 0.
The tick after it is a first tick (no previous cell), so it reads no held flag.

## Premises

**P1. The width.** After j reverse cells, U, V < 2^w with w = `tick_widths(q, T)[T − j]`, the bound the Euclid pass already uses for its own registers.
On the unswept tape the pair after j reverse cells is the Euclid pair (u, v) before Euclid tick T − j, because each Bézout cell inverts the Euclid tick that wrote it (`out/bezglue-w3/z3/rank_one_z3.py`, "Bezout cell inverts the Euclid tick").
The cells run on the swept tape, where that identity is false inside a shift run: the B1 sweep moves a run's single short tick, so at an interior boundary the swept pair is (u, 2v) where the Euclid pair is (u, v).
The two pairs agree at every run boundary.
At an interior boundary the forward tick before it was a double shift from (u, 4v), whose potential bounds max(u, 2v)² by the same B_t as the Euclid pair (`out/review-rankone/REVIEW.md`, "Exact premises and the width proof").
So the swept cofactors meet the same width; nothing beyond the stack's width premise is assumed.
At J = 12 the bound gives isqrt(B_372) = 5,792 < 2^13 after twelve cells, and the fold's premises hold at w = 13 (`out/bezglue-w3/toy/width-j12.txt`; the gate's `width_j.py` rechecks every k ≤ J from the generator).
On a held lane the stored word is H − R0 with H = S + R, so the words hold n + w + 1 bits, and two more wires are the add-or-shift's top.

**P2. The fold.** For X = hi·2ⁿ + lo < 2^w·q with (2^w + 2)·f ≤ 2ⁿ, m = hi + e with e = [lo + (hi + 1)f ≥ 2ⁿ] is ⌊X/q⌋, and lo + m·f mod 2ⁿ is X − m·q.
`rank_one.plan` turns the start off for a prime where this fails.

**P3. The last cell.** A cell after a B1 is a subtract (D1 of `docs/mux-proof.md`), so on a held lane of the last cell g0 = 1, and 2(X − (R >> 1)) + R0 = 2H − R for X = H − R0.

## The conversion

1. Fold R and S (P2): (r, mR) and (s, mS), mR = ⌊U·y/q⌋, mS = ⌊V·y/q⌋.
2. (U, V) from the cells on w-bit registers.
3. D = (U·s − V·r)/q is an integer, in (−V, U) when V > 0 since r, s < q, and 0 when V = 0 (then U = 1 by coprimality), so |D| < 2^w and w + 1 wires hold it in two's complement.
4. V·mR − U·mS = D, so mS becomes D in place on its own wires and one more: mS → U·mS mod 2^(w+1) (U is odd), complemented, plus V·mR, plus 1.
5. D mod U from |D| by restoring division, then U − rem on a negative D with rem ≠ 0.
6. The cells' tick mod U maps (mR, 0) in place to (U·mR mod U, V·mR mod U) = (0, D mod U), which the value of step 5 clears.
7. Step 5 is undone. D·q = U·s − V·r, so D → D·q mod 2^(w+1), minus U·s, plus V·r, on the low w + 1 bits of r and s, ends at 0. Step 2 is undone.

## Checks

- `out/bezglue-w3/tools/m_check.py` runs M(x, y) = (x, x·y mod q) on every x in [1, q) and every y: q = 61 (J = 3), 251 (J = 4), 4093 (J = 9) and 8191 (J = 10), with every other wire 0.
- The same script with a mutant: probing with hi in place of hi + 1, D with the wrong sign, a last cell that adds R on B1, no U − rem step, and words one bit short. Each fails on some lanes.
- `out/bezglue-w3/z3/rank_one_z3.py` proves steps 2 to 6, P2 at secp256k1 and P3, each with a wrong twin it refutes.
- `scripts/accept.py` runs the whole point addition at its toy primes; q = 61 takes the start with J = 3.
