#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
rm -rf dist build *.egg-info
python3 -m build
python3 -m twine check dist/*
echo "Build and distribution checks completed:"
ls -lh dist
