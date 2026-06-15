#!/usr/bin/env bash
# Run the offline end-to-end smoke test (no internet, no external bioinformatics
# tools required). Exits non-zero if any check fails.
set -euo pipefail
cd "$(dirname "$0")"
PYTHONPATH=. python tests/smoke_test.py
