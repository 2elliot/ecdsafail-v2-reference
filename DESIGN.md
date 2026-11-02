# Exact controlled point addition for `secp256k1-controlled-exact-v1`

## What the challenge actually requires

Read from a checkout of the ecdsafail-v2 benchmark repository (read-only), which is not public until the benchmark is published.

- Profile `profiles/secp256k1-controlled-exact-v1.json` sets `arithmetic_contract: EXACT_UNITARY`, `approximate_arithmetic: DISABLED`, `measurement_cleanup: DISABLED`, target `CLIFFORD_PLUS_T`.
- SPEC (`ecdsa specs/synthesized-v2/SPEC.md`) DEF-01 requires complete controlled addition of fixed public points. QNT-03 requires identity on control 0 and exact addition on control 1. QNT-05 requires the default exact profile to handle every reachable case, including infinity and R in {A, -A, -2A}. QNT-04 requires correctness on superpositions and with entangled references. SEM-03 rejects `R` and `Hmr`.
- Encoding (QNT-01) is little-endian affine coordinates with `(0,0)` as infinity.
- A submission is a `SubmittedMethod` JSON (`assurance/qualification/src/bundle.rs`) whose `circuit` is an `ArithmeticCircuit` of `X`, `CX`, `CCX`, `SWAP`, `Z`, `CZ`, `CCZ`, `NEG`, `CONTROLLED_NEG` (`contracts/ir/src/arithmetic.rs`).
- Service caps (`service/src/controller.rs`, `capacity/`, `contracts/api/src/quota.rs`) are 4096 wires, 10,000,000 operations, and a 100 MiB upload. The upload is parsed as plain JSON with no compression.
- `limits.json` caps compiled T at 3.25e9 for the whole attack. `compiler/lowering` lowers each CCX to the 7-T Toffoli.
- Qualification (`assurance/qualification/src/graph.rs`) needs a universal certificate checked by a `SoundTheoremChecker`. The production registry holds none, and `qcb diagnostic check` executes nothing. No artifact can be `is_qualified` in this repository today.

So Schrottenloher's published circuits (arXiv 2606.02235) cannot pass as published. They use MSB-only comparisons, pseudo-Mersenne shortcuts that skip rare carries, a fixed iteration count tuned for random inputs (failure up to 2^-13.3), measurement-based uncomputation, and they ignore doubling and infinity. Each of these violates one of the clauses above. What carries over is the exact structure underneath: the dialog Euclid algorithm and the out-of-band Bezout reconstruction that turns an inverse into an in-place multiplication.

## Designs compared

Costs are Toffoli estimates at n = 256 before implementation. Every design here is exact and uses only X, CX, CCX, SWAP.

**A. Exact dialog multiplier (chosen).** Point addition follows the paper's Algorithm 1. It uses an in-place multiplier M(x, y) = (x, xy) and its reverse, plus one out-of-place lambda-squared subtraction. M runs Euclid forward to turn x into a dialog bit-vector, then runs Bezout on (y, 0), then runs Euclid backward. Making it exact means three changes. The iteration count becomes the proven bound bitlen(q(q-1)) = 2n. Comparisons and modular reductions become full width. Uncomputation uses reverse circuits instead of measurement. Euclid widths still shrink, because u*v at least halves each iteration, so the width can be min(n, bitlen(q(q-1) >> i)). Estimated cost is about 58n^2 per M and about 136n^2, roughly 9M Toffolis, per addition. Qubits are about 9n.

**B. Kaliski or binary-EEA inversion plus schoolbook multiplications (Roetteler et al. 2017, Haner et al. 2020).** This design computes lambda with an explicit inverse, then uncomputes the inverse. That costs two inversion passes plus at least three multiply-and-uncompute pairs. Published exact counts are about 2.5e8 Toffolis per addition for Roetteler. Haner's improvements remain several times above design A. Rejected because it costs more on every axis that binds.

**C. Projective complete formulas (Renes, Costello, Batina).** These need no inversion and are complete. They are rejected because projective coordinates are not canonical. An in-place R -> R + A on the declared affine encoding still needs a normalizing inversion to uncompute the old point. That adds design A's cost on top of the formulas.

## Completeness on the fixed addend

For a fixed public A = (a, b), the generic kernel fails only on R in {inf, A, -A, -2A}. On those inputs x - a = 0 or a - x' = 0, and M is not a bijection there. All four inputs and their sums {A, 2A, inf, -A} are public constants. The circuit therefore:

1. computes flags f_k = c AND [R == E_k] for the four exceptional inputs, and gen = c AND NOT(f1 xor f2 xor f3 xor f4);
2. XORs E_k xor O_k into the accumulator under f_k;
3. runs the generic kernel controlled only through gen, using M^-1 ... M. With gen = 0, every gen-controlled step is the identity, and M^-1 followed by M cancels exactly at the circuit level;
4. uncomputes gen, then uncomputes each f_k from the output as c AND [R' == O_k]. This works because R -> R + A is a bijection.

M itself is total. For x = 0, Euclid emits an all-zero dialog and leaves u = q instead of 1. Bit 1 of u then distinguishes the two cases for free, and a controlled swap puts y back in place.

## Exact merged Bézout reduction

For every odd q >= 3 and canonical R, S, `Mod.dbl_add` computes (2S + g0 R) mod q with clean scratch. The Bézout step then applies the existing g1-controlled swap; its contract and [2, n, n] layout are unchanged.

q is odd with n = bitlen(q), and f = 2^n − q, so f is odd and 2^n < 2q.
m is any integer with bitlen(2f) ≤ m ≤ n, so 2f < 2^m, and the generator takes the smallest.
That range is never empty, because 2f < 2^n.
At secp256k1, m = 34 and h = n − m = 222.

R, S < q and g0 is a bit.
s = 2S + g0·R, so s ≤ 3q − 3.
t = s >> n, low = s mod 2^n, L = low mod 2^m, H = low >> m, and A = [H = 2^h − 1], with A = 1 when h = 0.
For j in {0, 1, 2}, c_j is the carry out of L + j·f on m bits, so c_0 = 0 and c_1 ≤ c_2 ≤ 1.
e = A·c_(t+1) for t ≤ 1, e = 0 for t = 2, and k = t + e.

### Gate order

Every wire below is a named position in `s`'s bits: bit 0 is a fresh wire `l0`, bits 1 to n are S's wires, and bit n + 1 is a fresh wire `top`.
So `low` is `[l0] + S[:n-1]`, t0 is `S[n-1]`, and t1 is `top`.

| # | Action | CCX at n = 256 |
|---:|---|---:|
| 1 | `add(R[1:] + [l0], S, top, ctrl=g0)`, with `l0` still zero as a's top bit | 3n + 1 = 769 |
| 2 | l0 ^= g0·R[0] | 1 |
| 3 | Prefix-AND chain over H, the carries of H + 1, giving A_i and A_h = A | h − 1 = 221 |
| 4 | Set selectors in place: X t1, CX t0 → t1, so t1's wire is [t = 0] and t0's wire is [t = 1] | 0 |
| 5 | Carry chain up on L + [t = 0]·f + [t = 1]·2f, giving c at the top | m = 34 |
| 6 | Fresh `e` ^= A·c, a CX when h = 0 | 1 |
| 7 | The same chain down, then undo step 4 | m = 34 |
| 8 | t += e in place: CCX(t0, e, t1), CX e → t0. The wires now hold k0 and k1 | 1 |
| 9 | Carry chain up on L + k0·f + k1·2f, giving c_k at the top | m = 34 |
| 10 | e ^= A·c_k, which clears `e`, a CX when h = 0 | 1 |
| 11 | For i from h − 1 down to 1: clear A_(i+1) from A_i and the old H_i, then H_i ^= c_k·A_i. Then H_0 ^= c_k, a CX | 2(h − 1) = 442 |
| 12 | For i from m − 1 down to 1: clear c_(i+1) from the old L_i and c_i, then L_i ^= addend_i xor c_i. Then clear c_1 from `l0` and the bit-0 addend, leaving `l0` alone | m = 34 |
| 13 | Fresh z = g0·¬k0 | 1 |
| 14 | `gt_into(R[1:], low[1:], t1, ctrl=z)`, which clears k1 | 2(n − 1) + 1 = 511 |
| 15 | Clear z | 1 |
| 16 | CX l0 → t0, so t0's wire holds S'[0], then l0 ^= g0·R[0], which clears `l0` | 1 |
| 17 | Free `l0`, `top`, `z`, then rotate `S` by one wire, as `moddbl` does | 0 |
| 18 | `cswap(g1, R, S)`, in `_bezout` | n = 256 |

A low chain has one AND per position, with the addend bit at position i being one of 0, the first selector, the second selector, or their XOR.
That is MAJ(L_i, addend, carry), and CX conjugation gives it with one CCX.
`_carries(C, H, 1)` already builds step 3, up to A_(h−1), and A_h is one more carry on top of it, so reuse what fits.
`_reduce`'s down loop has step 11's order: it uncarries position i + 1 before it writes position i.

**Cleanup order in steps 11 and 12.** A stored wire reads the old bit below it, so each one clears before that bit is written.
A_h goes first, right after step 10, because it is used only by steps 6 and 10 and writing H_(h−1) first would leave it dirty. With h = 2, H = 11 and c_k = 1, A_2 would stay 1.
A_1 is an alias of H_0 and has no wire to clear.
With h = 1, A is H_0, step 3 is empty, steps 6 and 10 read H_0, and step 11 is the single CX.
With h = 0, A is the constant 1, step 3 is empty, steps 6 and 10 are CX from c and c_k, and step 11 is empty.
Step 7 runs step 5 in reverse with nothing written, so it has no such order to get wrong.

The cell is inline, with no `C.apply` or `C.record` inside, because `lift` refuses a leaf with calls.

### Correctness and cleanup

**C1.** t is 0, 1 or 2, because s < 3q < 3·2^n.

**C2.** For j ≤ 2, low + j·f ≥ 2^n exactly when A·c_j = 1.
This holds because low + j·f = (H + c_j)·2^m + ((L + j·f) mod 2^m), and H + c_j reaches 2^h only from H = 2^h − 1 with c_j = 1.
The carry is one bit because 2f < 2^m.

**C3.** k = floor(s/q), and S' = s − kq = low + k·f − e·2^n.
With e = 0 and t ≤ 1, C2 gives low + (t + 1)f < 2^n, so 0 ≤ low + t·f < q.
With e = 1, C2 gives low + (t + 1)f ≥ 2^n, so S' ≥ 0, and S' < q is the same as s < (t + 2)q, which holds for t = 0 because s < 2^n < 2q, and for t = 1 because s < 3q.
With t = 2, S' = low + 2f < q is the same as s < 3q.

**C4.** Step 5's addend is (t + 1)·f below 2^m for t ≤ 1 and 0 for t = 2, so step 6 computes the e above.
Step 7 finds the selectors as step 5 left them, because step 6 touches only `e`.

**C5.** e = A·c_k, on every reachable state.
If e = 1 then c_k = c_(t+1) and A = 1.
If e = 0 and t ≤ 1 then c_k = c_t ≤ c_(t+1), so A·c_k ≤ A·c_(t+1) = 0.
If t = 2 then low + 2f < q < 2^n, so A·c_2 = 0 by C2.
So step 10 clears `e` while A and c_k are still live.
This is the sketch's trap: `e` is the overflow out of H as well as k − t, so it needs no clear of its own from the final state.

**C6.** S' = ((H + c_k) mod 2^h)·2^m + ((L + k·f) mod 2^m), which follows from C3, and S' < 2^n.
Steps 11 and 12 write exactly that, and k0, k1 are exclusive because k ≤ 2.
Every stored A_(i+1) and c_(i+1) reads 0 after its clear, because each clear runs on the same inputs that set the wire: the old bit, and a lower chain wire that has not been cleared yet.

**C7.** k1 = 1 needs g0 = 1, because g0 = 0 gives s = 2S < 2q.

**C8.** Step 14 clears k1.
If z = 1 then k is even, so S' and R have the same parity and [R > S'] = [R >> 1 > S' >> 1].
With k = 0, S' = 2S + R ≥ R.
With k = 2, S' = R − 2(q − S) ≤ R − 2.
So the compare equals k1.
If z = 0 then k1 = 0, by C7 or because k0 = 1.
k = 1 can have R > S', so the compare needs z and not g0 alone.

**C9.** S'[0] = k0 xor g0·R[0], because 2S is even and q is odd.
So step 16 leaves S' on the rotated wires and clears `l0`.

**C10.** For each g, S → 2S + g0·R mod q is a bijection on [0, q), because q is odd.
The circuit is a permutation that takes canonical states with zero scratch onto canonical states with zero scratch, so the reversed cell that M⁻¹ runs is exact on canonical inputs too.

Nothing here uses primality, n = 256, or a run of ones.
The run of ones only makes m small.

## Data shape

- **Wire.** An integer index. **Register.** A Python list of wires, least significant bit first. Shifts are relabels, implemented as list rotations with no gates.
- **Circuit.** A flat `array('i')` of (opcode, a, b, c), a free-list wire allocator that requires ancillas to be returned clean, and a high-water mark that becomes `num_wires`.
- **Block.** A recorded op range with input layout, output layout and internal wires. It can be instantiated forward or reversed onto any matching layout. That is how M^-1, Euclid-backward and every compute/uncompute pair are emitted without hand-written inverse code.
- **Simulator.** A bit-sliced simulator in which each wire is a Python int and bit j is test lane j. It shares no code with the generator's arithmetic.

## Why the permutation is enough for QNT-04

The circuit contains only X, CX, CCX and SWAP. It is therefore a permutation matrix with no phases. If every promise basis state |c, R, 0> maps to |c, R + cA, 0>, the same matrix acts exactly on any superposition of promise states and on any extension with a reference system. Basis exactness plus clean scratch is the whole unitary contract. Outside the promise, the circuit is still a permutation, so QNT-06's valid-channel requirement holds by construction.
