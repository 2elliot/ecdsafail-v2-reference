#!/bin/bash
# Reproduces every receipt in receipts/ outside baseline/, and the tracked artifacts/.
# Needs the ecdsafail-v2 sibling checkout (read-only) and ~/.cargo/bin for the harness.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-~/ecdsa_venv/bin/python}
export PATH="$HOME/.cargo/bin:$PATH"
exec > >(tee receipts/run_all.log) 2>&1

$PY -I scripts/model_check.py > receipts/model_check.jsonl
$PY -I scripts/check.py --toy-bits 4 6 7 8 9 10 --toy-per-size 6 --out receipts/toy_exhaustive.json > /dev/null
$PY -I scripts/mutation_control.py > receipts/mutation_control.json
$PY -I scripts/component_check.py > receipts/secp256k1_components.json
/usr/bin/time -l $PY -I scripts/check.py --secp --secp-random 120 --seed 2026 --out receipts/secp256k1_sampled.json > receipts/secp256k1_sampled.log 2>&1
$PY -I scripts/breakdown.py > receipts/secp256k1_breakdown.json

$PY -I scripts/emit.py --curve toy:127 --addend-index 5 --out artifacts/toy127_circuit.json
$PY -I scripts/emit.py --out artifacts/secp256k1_G_circuit.json.gz
$PY -I scripts/emit.py --curve toy:127 --addend-index 5 --hier --out artifacts/toy127_circuit.hier.json
$PY -I scripts/emit.py --hier --out artifacts/secp256k1_G_circuit.hier.json
$PY -I scripts/make_cases.py --curve toy:127 --addend-index 5 --out artifacts/toy127_cases.json
$PY -I scripts/make_cases.py --random 250 --out artifacts/secp256k1_cases.json

(cd harness && CARGO_TARGET_DIR="$PWD/target" cargo build --release --offline --quiet)
H=harness/target/release/qec-harness
$H --circuit artifacts/toy127_circuit.json --cases artifacts/toy127_cases.json --lower --write-submission artifacts/toy127_submission.json > receipts/harness_toy127.json
gzip -dc artifacts/secp256k1_G_circuit.json.gz | /usr/bin/time -l -o receipts/harness_secp256k1.time $H --circuit - --cases artifacts/secp256k1_cases.json --lower > receipts/harness_secp256k1.json
echo "all receipts written in ${SECONDS}s"
