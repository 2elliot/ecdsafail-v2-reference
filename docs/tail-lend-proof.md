# Last-cell tape lend

Named limit: **`I2.model`'s zero-caller row narrowed from the reference's gate**.

The pass entry is U=q and 1<=V<q with clean scratch; M remaps x=0 to1. At the terminal subtract, coprimality and equality imply U=V=1, so that swept cell is (g0,b)=(1,0). Subsequent cells are (0,1). The offset sweep restores the reference tape; the B1 sweep exchanges b only between adjacent no-subtract cells and cannot alter this terminal suffix. The proven tick bound therefore implies g0 XOR b=1 at the final cell. At secp256k1 this uses the existing 384-tick bound, not a sampled shorter bound.

On the merged path, the last cell is read by the first two reverse-order Bezout ticks. After the second, CX(g0,b) then X(b) makes b zero. When the final cell is in the uncoded tail, no later tick or codec reads it, so the allocator may reuse its physical wire. The code sets its tape slot to None while it is lent. After Bezout, a newly allocated zero wire receives X then CX(g0,b), restoring the logical final b. The final g0 remains live and untouched. The pool's tape mapping is rebuilt from G before inverse Euclid, so restoration need not use the old wire number. This boundary remapping is part of the circuit and its inverse.

The original physical wire is not promised zero while allocated to other work: it may hold a decoded block or cell scratch. Zero is required at release, and on the fresh allocation before restoration. Any probe while the wire is lent must distinguish the logical absent b from the reused physical storage.

The lend is enabled only on the merged path, with at least two ticks and an uncoded final cell. Fallback passes and fully coded tails retain their prior behavior. No tick bound, codec, cell, gate row or arithmetic operation changes. The added logical operations per enabled M are two CX and two X; physical emitted counts also reflect the restored register mapping.

An independent checker checks the actual generated and decoded tape on every x at all small gate primes (model, component and toy domains), with x=0 remapped; it also checks M both ways from independently expected clean states. The data-cell CX/X control yields nonzero b; deleting restoration makes M fail. Secp inputs are directed/seeded, including known long runs, not exhaustive. Counts and allocator peaks are measured separately. The wider all-input claim rests on the suffix lemma and bound.

Inherited compose.rs refusal for vacated-U/V tape, held mux lift refusal, charge and qualification limits remain. No checker admission is implied.
