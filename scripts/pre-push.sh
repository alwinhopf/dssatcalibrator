#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

PYTHON_BIN="$(command -v python || command -v python3 || true)"
if [ -n "$PYTHON_BIN" ]; then
  "$PYTHON_BIN" -m pytest -m "not slow" -q
else
  python -m pytest -m "not slow" -q
fi
