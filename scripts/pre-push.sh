#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="$(command -v python || command -v python3)"
"$PYTHON_BIN" -m pytest -q -m "not slow" --junit-xml=test-results.xml 2>&1 | tee fast-tests.log
