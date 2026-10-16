#!/usr/bin/env bash
# Reproduce all data products used by notebooks/figures.ipynb.
# Total run time on a laptop: roughly 15-20 minutes.
set -euo pipefail
cd "$(dirname "$0")/.."

python scripts/compute_annual_beams.py     --out data/annual_beams.npz
python scripts/compute_skymaps.py          --out data/skymaps.npz
python scripts/compute_systematics.py      --out data/systematics.npz
python scripts/compute_density_profile.py  --out data/density_profile.npz
python scripts/check_convergence.py        | tee data/convergence.txt
