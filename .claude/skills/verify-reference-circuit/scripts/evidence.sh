#!/bin/sh
# Copies the last gate run into verify-out/evidence/<label>-<UTC time>/ and prints that path.
# Usage: evidence.sh <label> [extra file ...]
set -eu
cd "$(dirname "$0")/../../../.."
[ $# -ge 1 ] || { echo "usage: $0 <label> [extra file ...]" >&2; exit 2; }
if command -v sha256sum > /dev/null; then sha256() { sha256sum "$1"; }
elif command -v shasum > /dev/null; then sha256() { shasum -a 256 "$1"; }
else echo "evidence.sh needs sha256sum or shasum on PATH" >&2; exit 1; fi
label=$1; shift
dest=verify-out/evidence/$label-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$dest"
[ -f verify-out/accept.json ] && cp verify-out/accept.json "$dest/"
[ -d verify-out/accept/logs ] && cp -R verify-out/accept/logs "$dest/logs"
for f in verify-out/accept/*.json; do [ -f "$f" ] && cp "$f" "$dest/"; done
if [ -d verify-out/accept/emit ]; then
  for f in verify-out/accept/emit/*; do
    printf '%s  %s bytes  %s\n' "$(sha256 "$f" | cut -d' ' -f1)" "$(wc -c < "$f" | tr -d ' ')" "$f"
  done > "$dest/emitted.sha256"
fi
for f in "$@"; do cp "$f" "$dest/"; done
{ git rev-parse HEAD; git status --porcelain; } > "$dest/git.txt" 2>/dev/null || true
echo "$dest"
