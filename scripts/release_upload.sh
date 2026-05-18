#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
if [[ -z "${TWINE_USERNAME:-}" || -z "${TWINE_PASSWORD:-}" ]]; then
  echo "TWINE_USERNAME and TWINE_PASSWORD must be set in the environment." >&2
  exit 1
fi
TARGET="${1:-pypi}"
if [[ "$TARGET" == "testpypi" ]]; then
  python3 -m twine upload --repository testpypi dist/*
else
  python3 -m twine upload dist/*
fi
