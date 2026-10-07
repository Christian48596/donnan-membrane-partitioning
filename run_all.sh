#!/usr/bin/env bash
# Reproduce the numerical results, scientific checks and figures.
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONPATH="$PWD/code"
export MPLBACKEND=Agg
mkdir -p results

reproduce() {
  echo "[1/7] Scientific tests"
  python -m unittest discover -s tests -v
  echo "[2/7] Equilibrium sorption and association"
  python code/analysis.py
  echo "[3/7] Fixed-charge calibration and transfer"
  python code/dle.py
  echo "[4/7] Lithium/sodium exchange"
  python code/ion_specific.py
  echo "[5/7] Spatial-average approximation checks"
  python code/checks.py
  echo "[6/7] Figure data and artwork"
  python code/make_figures.py
  echo "[7/7] Numerical and export verification"
  python code/verify_outputs.py
  echo "Reproduction complete. Results: results/. Figures: figures/."
}

reproduce 2>&1 | tee results/reproduction.log
