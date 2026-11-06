#!/usr/bin/env bash
# Build the staphmix and staphmix-db conda packages from this checkout.
# The database directory is never copied into the source tree of the code
# package, so builds stay fast.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${ROOT}/dist/conda"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# --- code package -----------------------------------------------------------
rsync -a \
  --exclude 'db' \
  --exclude '.git' \
  --exclude '__pycache__' \
  --exclude '*.pkl' \
  --exclude 'dist' \
  "$ROOT/" "$WORK/src/"

conda build "$WORK/src/conda/recipe" --output-folder "$OUT"

# --- database package -------------------------------------------------------
# Requires conda/recipe-db/meta.yaml to point at the DB tarball or directory.
if [[ "${BUILD_DB:-0}" == "1" ]]; then
  conda build "$ROOT/conda/recipe-db" --output-folder "$OUT"
fi

echo
echo "Packages written to $OUT"
echo "Attach them to a GitHub release, or upload them to an anaconda.org channel."
