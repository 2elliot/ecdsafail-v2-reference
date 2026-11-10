#!/bin/sh
# Frozen cost harness: emits a tree's secp256k1 circuit (python3.12, with the tree's depth pass),
# wraps it, runs the tree's post-pass qec/gate_pass.py, and prints one line:
#   tag wires W CCX C T t raw_ticks R ticks K qubit_time Q T_headroom H
# Usage: sh lab/tools/measure.sh <tree> <tag> [MODULE.FLAG=value ...]
# KEEP=<path> keeps the post-passed method file. Deterministic: the tick model is exact.
set -e
here=$(cd "$(dirname "$0")" && pwd)
tree=$(cd "$1" && pwd); tag=$2; shift 2
d=$(mktemp -d "${TMPDIR:-/tmp}/measure.XXXXXX")
python3.12 -I -B "$here/emitvar.py" "$tree" "$d/c.json" "$@" > /dev/null
python3 -I "$here/wrap.py" "$d/c.json" "$d/raw.json" > /dev/null
gp=$(python3 -I "$tree/qec/gate_pass.py" "$d/raw.json" "$d/out.json")
cc=$(python3 -I "$here/ccx.py" "$d/out.json")
python3 -I - "$tag" "$gp" "$cc" <<'PY'
import json, sys
tag, gp, cc = sys.argv[1], json.loads(sys.argv[2]), sys.argv[3].split()
ccx = int(cc[cc.index("CCX") + 1]); w = gp["wires"]; t = gp["ticks_out"]
print(f"{tag} wires {w} CCX {ccx} T {7*ccx} raw_ticks {gp['ticks_in']} ticks {t} qubit_time {w*t} T_headroom {50465933-7*ccx}")
PY
if [ -n "$KEEP" ]; then mv "$d/out.json" "$KEEP"; fi
rm -rf "$d"
