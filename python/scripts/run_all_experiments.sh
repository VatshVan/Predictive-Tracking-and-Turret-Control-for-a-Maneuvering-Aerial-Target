#!/usr/bin/env bash
# Regenerate every result, table and figure used in the report and paper.
# Usage (from python/):  bash scripts/run_all_experiments.sh
# Afterwards, in MATLAB (matlab/): parity_check, time_controller, run_scenarios
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
step() { echo "[$(date +%H:%M:%S)] $*"; }

step "unit tests";               python -m pytest -q tests
step "design analysis";          (cd scripts && python design_analysis.py > ../results/design_analysis.txt)
step "scenario figures";         python scripts/plot_scenarios.py > results/scenarios.txt
step "parity logs";              python scripts/export_parity.py > /dev/null
step "Monte Carlo ablations";    FORCE=1 bash scripts/run_ablations.sh
step "logged ensemble";          python scripts/mc_ensemble.py --n 300 --f-ctrl 100
step "stress noise ensembles";   python scripts/stress_ensemble.py --reps 50 > /dev/null
step "rate study";               python scripts/rate_study.py > results/rate_study.txt
step "statistics";               python scripts/mc_statistics.py > /dev/null
step "experiment figures";       python scripts/plot_experiments.py > /dev/null
step "interface robustness";     python scripts/interface_robustness.py --n 300 > /dev/null
step "report figures";           python scripts/plot_report_extras.py > /dev/null
step "paper numbers";            python scripts/export_paper_numbers.py
step "done"
