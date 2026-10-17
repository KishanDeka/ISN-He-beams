#!/usr/bin/env bash
# Regenerates everything in data/ (takes ~15 min on my laptop).
set -euo pipefail
cd "$(dirname "$0")/.."

python scripts/compute_annual_beams.py     --out data/annual_beams.npz
python scripts/compute_skymaps.py          --out data/skymaps.npz
python scripts/compute_density_profile.py  --out data/density_profile.npz
python scripts/check_convergence.py        | tee data/convergence.txt
