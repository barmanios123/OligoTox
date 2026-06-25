#!/usr/bin/env bash
set -euo pipefail
cd ~/oligotox

if [ -f scripts/train_v3.py ]; then
  python scripts/train_v3.py
else
  echo "scripts/train_v3.py not found"
  exit 1
fi
