#!/bin/sh
# Read-only. Exits 1 when this checkout cannot run the gate.
cd "$(dirname "$0")/../../../.." || exit 1
PY=${PY:-python3}
bad=0
say() { printf '%-5s %s\n' "$1" "$2"; }

if ver=$($PY -c 'import sys; print("%d.%d" % sys.version_info[:2]); sys.exit(sys.version_info < (3, 9))' 2>/dev/null); then
  say ok "python $ver at $(command -v "$PY")"
else
  say FAIL "$PY is missing or older than 3.9"; bad=1
fi
if $PY -I -c 'import sys; sys.path.insert(0, "."); import qec.point_add, qec.emit, qec.flatten' 2>/dev/null; then
  say ok "qec imports with the standard library only"
else
  say FAIL "qec does not import; run: $PY -I -c 'import sys; sys.path.insert(0, \".\"); import qec.point_add'"; bad=1
fi
[ -f scripts/accept.py ] && say ok "scripts/accept.py present" || { say FAIL "scripts/accept.py missing"; bad=1; }
if [ -e verify-out ] && [ ! -w verify-out ]; then say FAIL "verify-out is not writable"; bad=1; else say ok "verify-out writable at $PWD/verify-out"; fi
if head=$(git rev-parse --short HEAD 2>/dev/null); then
  say ok "git $head on $(git rev-parse --abbrev-ref HEAD), $(git status --porcelain -- qec | wc -l | tr -d ' ') changed files under qec/"
else
  say info "not a git checkout; accept.py check R1 will have nothing to compare"
fi
running=$(pgrep -f 'scripts/accept.py' | tr '\n' ' ')
[ -n "$running" ] && say warn "accept.py already running (pids $running); two gates in one checkout share verify-out/accept" || say ok "no other accept.py running"
if [ -n "${ECDSAFAIL_V2:-}" ] && [ -f "$ECDSAFAIL_V2/contracts/ir/Cargo.toml" ] && [ -d ../ecdsafail-v2 ]; then
  say ok "benchmark checkout at ../ecdsafail-v2; the harness feature can run"
else
  say info "no benchmark checkout (ECDSAFAIL_V2 unset or ../ecdsafail-v2 missing); the harness feature cannot run"
fi
exit $bad
