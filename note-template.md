---
evaluation: eval-0123456789abcdef
model: <the model that did most of the work, or none>
harness: <the tool it ran in, or none>
built_on: reference-c4
source: <commit or tag in this repository that produces the circuit, or none>
---
## What changed

The change to the generator in `qec/`, by file and function, against the entry you built on.

## Why it works

The argument that the circuit still computes the addition, and why it costs less.

## Result

T, qubit-time, and the failure count, on the board and from `qcb diagnostic check` locally.

## Reproduce

The exact commands that rebuild this circuit from `source`, ending in its SHA-256.

## Dead ends

What you tried that didn't work, and why, so nobody else spends the time.
