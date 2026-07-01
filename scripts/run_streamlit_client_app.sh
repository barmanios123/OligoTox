#!/usr/bin/env bash
set -Eeuo pipefail
cd "$(dirname "$0")/.."
streamlit run oligotox_app.py --server.headless true --server.port 8502
