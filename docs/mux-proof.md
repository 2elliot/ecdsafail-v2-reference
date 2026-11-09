# The mux merge on a held offset

The cell is in `qec/two_bit.py`, and the multiplier change is in `qec/point_add.py`.
Each statement is numbered so that it can be tested as worded.

## Credit

The representation T = S + R between a doubling and the following subtract is section VI of arXiv:2609.05625.
No novelty is claimed for the offset.
That paper's complement form leaves (c − 1)/p of inputs off the canonical range, so it is not used here.
The subtract below is a true integer map on the domain stated in D1 to D4.

## Notation

q is an odd prime of n bits and f = 2ⁿ − q, which is odd.
A prime takes the new path when `merged` holds: 4f < 2ⁿ and corner ≤ 4f, where corner = f ⊕ 2f ⊕ 3f.
Only 4f < 2ⁿ is used below; the second clause keeps the set of primes on the path where it was.
w = bitlen(4f).
Every other prime keeps today's `step`, with no sweep.

A cell is (g0, b), and the symbols are A0 = (1, 0), A1 = (1, 1), J = (0, 1), B1 = (0, 0).
A subtract cell is one with g0 = 1.
Ticks run in Bézout order, the reverse of the tape, so tick t reads cell T − 1 − t.

| Lane | Condition | Map |
|---|---|---|
| Plain add | g0 = 1, previous cell not B1 | S = 2T + R mod q, then swap on b |
| Held | g0 = 1, previous cell B1 | S = 4T + R mod q, then swap on b |
| J | J | S = 4T mod q |
| B1 | B1, not the last tick | no reduction, H = 2T + R is held |
| Last B1 | B1, last tick | S = 2T mod q |

## Domain

Every statement of exactness below names this domain.

**Top level.** R is on the curve or infinity, in the canonical encoding, with zero scratch.

**Tick level.** T and R are canonical, scratch is zero, and:

- **D1.** A cell after a B1 is a subtract cell.
- **D2.** A B1 that is not the last tick has next g0 = 1 on entry.
- **D3.** A B1 on the last tick has R = 0.
- **D4.** On a held tick, the input is the stored form of S4, not a canonical word.

Today's cell was exact on every symbol and canonical R, S.
This tick is not, so a per-cell all-input check does not apply to it.
D1 to D3 are properties of a swept Euclid tape, by S1 to S3.

## The tape

**S1. Run shape.**
In Bézout order, each maximal run of cells with g0 = 0 holds at most one B1, and it is the run's first cell.
`tick_tape` writes B1 exactly when v ≡ 2 mod 4, which leaves v odd, so the next tape cell is a subtract or the tape ends.
Reversing the tape gives the statement.

**S2. The sweep.**
`_b1_to_run_end` walks t = 0 to T − 2 and swaps b between cells t and t + 1 when both have g0 = 0.
It costs 3 CCX per adjacent pair and one clean wire.
It reads g0 and never writes it.
In each run the b values rotate left by one, so a B1 moves to the run's end and the Js shift down.

After the sweep every B1 is followed by a subtract cell or is the last cell, which gives D1 and D2 on the pass's own tape.
The number of ticks, the alphabet, the symbol counts, `tick_bound` and `tick_widths` are unchanged, and `model.py` is not edited.
B1 and J only double S, R is unchanged across a run, and doublings commute, so S and R after each run equal their unswept values.
The pass's result is therefore `tick_bezout`'s.
`_body` applies the sweep before the first tick and its reverse after the last, so the tape handed to the inverse Euclid block is the original.

**S3. The last tick.**
If the last cell is a B1 then R = 0 on it.
The pass ends with R = 0 by the existing proof of the recurrence, and a B1 does not change R.

## The integers

**S4. Stored form.**
Let H = 2T + R on a B1 that is not the last tick.
Then H < 3q < 2ⁿ⁺², and H ≡ R mod 2.
After the tick, in the rotated S register, wires 1 to n − 1 hold H's bits 1 to n − 1, wire 0 holds H's bit n, and the next cell's g0 wire holds H's bit n + 1.
H's bit 0 is not stored, because it equals R₀, which stands in R.
The previous cell's wires are unchanged, and they still read B1.
So a held tick can recover H from its own inputs.

**S5. The held tick's word.**
Let W = H − R₀ − (R >> 1), which is 2T + (R >> 1).
Then 0 ≤ W < 2ⁿ⁺².
H − R₀ is H with bit 0 cleared, so the subtract needs no carry in.
Let V = 2W + R₀.
Then V = 4T + R, and 0 ≤ V < 5q.
Nothing on this tick is negative, which is the reason for holding.

**S6. The reduction.**
Let V ≥ 0, V < 5q, low = V mod 2ⁿ and K = V >> n, so K is at most 4.
Let probe be f, 2f, 3f, 4f, f for K = 0, 1, 2, 3, 4, and e = [low + probe ≥ 2ⁿ] and k = K + e.
Then k = ⌊V/q⌋, k ≤ 4, [low + kf ≥ 2ⁿ] = e, and S = low + kf − e2ⁿ = V − kq.

For K ≤ 3, V ≥ (K + 1)q exactly when low + (K + 1)f ≥ 2ⁿ.
For K ≤ 2, V < (K + 1)2ⁿ ≤ (K + 2)q because (K + 2)f ≤ 4f < 2ⁿ.
For K = 3, V < 5q bounds k by 4.
For K = 4, low < 2ⁿ − 5f, so low + 4f < 2ⁿ, e = 0 and k = 4.
For K ≤ 3, if e = 0 the write has no carry, since kf ≤ probe, and if e = 1 the carry is the probe's, with kf = probe.

K = 3 needs a true 4f, where the earlier low-chain cell, which is not in this repository, used corner as a don't-care.
A separate check, also not in this repository, shows that corner fails on the held lane at secp256k1.

**S7. The clear.**
Write k = k₀ + 2k₁ + 4k₂, on wires k0, k1, k2.
S₀ = V₀ ⊕ k₀ on every lane, since f is odd, so the cell does not write bit 0 and `cx(l0, k0)` turns k0 into S₀.

- **Plain add.** V = 2T + R < 3q, so k ≤ 2, and k₂ = 0. Odd k is k = 1. For even k, k = 2 exactly when R >> 1 > S >> 1, because S and R agree on bit 0 and k = 0 gives S ≥ R while k = 2 gives S = R − 2(q − T) < R.
- **Held.** S ≡ R + kf mod 4, because 4T ≡ 0 and q ≡ −f. So k₁ = S₁ ⊕ R₁ ⊕ k₀R₀ ⊕ (f₁ · k₀), where f₁ is bit 1 of f. If k₀ = k₁ = 0 then k is 0 or 4, and k = 4 exactly when R >> 1 > S >> 1, by the same argument with S = R − 4(q − T) < R. There is no offset of f.
- **J.** S ≡ kf mod 4, as today.

On a held lane `cswap(p, k1, k2)` first exchanges the two wires, so the compare's target carries k₂, and wire k2 carries k₁.
`mark_even` makes z = (g0 ⊕ wire k2) · ¬k0, which on a held lane is 1 exactly when k is 0 or 4.
On other lanes wire k2 is 0, so z is today's g0 · ¬k0.

## The add and the shift

**S8. Complementary controls.**
Every tick adds or shifts: off J it adds R >> 1 to the word, and on J it shifts the word left.
On a held lane the word is complemented by `cx(p, x)` over all n + 2 wires before and after, and ~(~x + a) = x − a mod 2ⁿ⁺², which is S5's W.
A B1 adds like a plain add, with l0 left at 0 by g0 = 0.

This section describes `add_or_shift_cdkm`, the add the tick runs when `CDKM_ADD` is set (2elliot frontier/v1).

**S9. The add, with its sum writes off the carry chain.**
Let a_i be `ops[i]`, t_i be `tgt[i]` = word[i], c_i the carry into bit i (c_0 = 0), x_i = a_i ⊕ c_i, and u = ¬J.
The add is CDKM's ripple, and it runs on every lane. Only the sum writes depend on u.

- **Up.** MAJ_i leaves t_i ⊕ a_i on `tgt[i]`, x_i on `carries[i]` and c_(i+1) on `ops[i]`. MAJ_(i+1)'s `cx(a, c)` XORs into `ops[i]`, the target of MAJ_i's CCX, and MAJ_i reads neither a_(i+1) nor `ops[i]` as a control. So it commutes with that CCX and is issued first. The map is CDKM's.
- **Pre-write.** After MAJ_i, `tgt[i]` ⊕= u · ¬x_i, by X, CCX, X on `carries[i]`.
- **Down.** The restore CCX on bit i XORs x_i · `tgt[i]` into `ops[i]`. That is x_i · (t_i ⊕ a_i ⊕ u · ¬x_i) = x_i · (t_i ⊕ a_i), CDKM's value, because x_i · ¬x_i = 0. So `ops[i]` and `carries[i]` come back exactly as in CDKM, on every lane. The restore `cx(ops[i + 1], ops[i])` of `carries[i + 1]` XORs into the next CCX's target, so it is issued after that CCX.
- **Sum.** `tgt[i]` then takes u from a ¬J wire and a_i from the restored `ops[i]`. It ends as t_i ⊕ u · ¬x_i ⊕ u = t_i ⊕ u · x_i. Off J that is the sum bit, and on J it is t_i.

So off J, word[:n − 1] gains R >> 1 with carry c_(n−1) out, and on J it is unchanged.
Each bit costs 1 CCX up, 1 pre-write and 1 restore, as before.
The carry chain holds only MAJ's CCX going up and the restore's CCX going down. Going up, the chain enters each CCX as a control. Going down, it must wait for the CCX above to release a control. The measured costs are 7 ticks per bit up and 11 per bit down.

**S10. The top and the shift.**
The top is as before. `carry` = ¬J · c_(n−1). Then t = carry · word[n − 1], word[n + 1] ⊕= t · word[n], word[n] ⊕= t and word[n − 1] ⊕= carry. Then t is cleared from the new word[n − 1] under an X, and `carry` is cleared.
On J, the shift rotates word[:n + 1] left by one, as the old ladder did. word[n + 1] is untouched, and on the domain word[n] = 0 there.
That rotation is the reversal of word[0..n], with pairs (i, n − i), followed by the reversal of word[1..n], with pairs (i + 1, n − i). Each pair is a `cswap` under a copy of J.
Pair (i + 1, n − i) of the second reversal reads only the first reversal's pairs (i + 1, n − i − 1) and (i, n − i), and both are issued before it.
The down pass never touches word[i] after bit i. The top wires are final before it starts. So both pairs at bit i are issued right after that bit, and the swaps follow the down pass with no ladder.
The shift's swaps cost n CCX.
It spends 3(n − 1) + 5 + n CCX, which is 1,026 at n = 256. Over the whole circuit that is 579,096 fewer CCX than the blocked-lookahead add (`add_or_shift_la`), which stays in the tree behind `LOOKAHEAD_ADD`.

The add's ¬J copies serve the pre-writes and are cleared once the up pass ends. From then on J itself, flipped, carries ¬J.
The held lane's second complement keeps the tick's mirror-pair order. Those CXs share the control p and commute, so any order gives the same map.

## The lanes

**S11. B1 costs no CCX to read inside a bracket.**
Between `_begin` and `_end`, b holds J and g0 is true, so g0 ⊕ b is 0 on a B1 and 1 on the other three symbols.
The first bracket stays open through the first chain, so `one` = g0 ⊕ J is two CXs, up and down.

**S12. The previous cell.**
[previous cell is B1] needs a real AND, so `_after_b1` costs 1 CCX and `_drop_after_b1` 1 more, read twice per tick, 4 CCX.
Its wire does not stay through the reduction, because the reduction is the peak.

**S13. B1 is inert in the reduction.**
On a B1 that is not last, before either chain, k0 moves to l0 and k1 to the next cell's g0, whose 1 is first cleared.
This is sound because l0 = 0 there by g0 = 0, and next g0 = 1 by D2.
Then every selector is 0, so every logical carry of both chains is 0, and neither chain reads L[0] into anything.
The first chain's carries sit on H's wires, which equal the logical carries only when P = 1.
When P = 0 they may hold junk, which nothing reads except through P and which the first chain's reverse removes.
So e = 0, and every carry of the second chain is 0, both on its clean wires and on the selectors it perches on, which are 0 here.
`P.write(0)` and `M.write(0)` restore their words.
At the end `cx(l0, k0)` gives wire k0 the bit, and `ccx(L, k0, l0)` clears l0, which S11 keeps to 1 CCX.
That is S4's form.

**S14. Entering a held tick.**
`cswap(p, s[0], spill)` moves H's bit n to word[n] and leaves word[0] = 0, which is H − R₀ by S5.
`cswap(p, g0, high)` moves H's bit n + 1 to word[n + 1], and `cx(p, g0)` restores g0 = 1.
From there g0 is true, so the tape is whole again.

**S15. The last tick.**
With nxt = None, `one` = 1 and nothing is hidden.
By S3 a B1 there has R = 0, so its add is 0 and V = 2T < 2q, and the reduction runs as on any other lane.
z = 0 on it, and k₁ = 0, so the clear is S₀ alone.

**S16. The first tick.**
With prev = None there is no held lane.
Two gates are then inert on the whole domain: `ccx(t, word[n], word[n + 1])`, since word[n] = 0 on every lane, and `ccx(three, e, k2)`, since K = 3 with e = 1 needs V ≥ 4q, which a tick with no held lane cannot reach.
Those gates remain in the general tick, but M no longer calls that tick first.
M starts with S = 0 and its terminal tape cell is A0, A1 or J, never B1:
the final active cell is a subtract, and every padded cell after it is J;
the offset sweep changes b only on subtracts and the B1 sweep cannot introduce
B1 into that all-J tail. On A0 and A1 the first tick gives (R, R), and on J it
gives (R, 0). Thus M uses S XOR= g0*R, one CCX per bit, with R, cell and next
unchanged and no scratch. S keeps its wire list for the next tick. This uses
the existing tick bound, including the proved 384 bound at secp256k1, and reads
only the terminal cell's g0, so the tail b lend is unaffected. The inverse
uses the same XOR on the image of this boundary. Fallback steps are unchanged.

## The whole circuit

**S17. Both directions.**
`_body` is one recorded block, so M⁻¹ is its reverse, gate for gate.
Forward maps the domain (x, y) with canonical y and clean scratch to (x, x′y mod q), which is every canonical y again because x′ is invertible, so the reverse runs on exactly the forward image.
`_kernel` applies M⁻¹ first and then M.
With gen = 0 the two cancel at circuit level, as before.
The reverse of the pass applies the sweep, then the ticks backward from the last, then the sweep's reverse, so a reversed tick reads the swept tape and the stored form of S4.

**S18. Other moduli.**
The generator runs at every q in the gate.
Primes in the gate on the new path: 7, 13, 31, 61, 97, 127, 251, 397, 433, 509, 823, 907, 967, 1021, 4093, secp256k1.
Primes in the gate on today's path: 3, 5, 11, 43, 67, 79, 163, 211, 349, 547, 577, 733.
`f₁ = 1`, the `f >> 1 & 1` branches, is reached only at toy primes, because secp256k1 has f ≡ 1 mod 4.
Empty H and one-bit H are reached only at toy primes.
Both prefixes at their real lengths, H = 221 bits and the quiet run of 20, are reached only at secp256k1, and `Prefix` is byte-identical to the earlier low-chain cell's.

**S19. Peak.**
The cell needs 26 clean wires, the same as the earlier low-chain cell, at the top of the second chain.
`high` costs one more there, and clearing c₁ once c₂ exists pays it back at +2 CCX, by the same carry rebuilt on the way down.
c₁₂ still sits on `one`, because a perch can now read through other selectors, here k2.
The two held bits need no clean wire: S4 keeps them on wires the tick already owns.
The add's 10 scratch wires (`cin`, `carry`, `t`, 3 ¬J copies and 4 J copies) are live only during the add. The circuit keeps 1,249 wires with this add, the same as with the lookahead add.

**S20. Limits of this proof.**
It proves maps, and the wires by the arguments above and by test, not by a formal method.
`lift` on the new cell has not been tried, and its result belongs in the limit list either way.
Nothing here passes the checker, and `is_qualified` stays false.
