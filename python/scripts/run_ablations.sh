#!/usr/bin/env bash
# Re-run the final baseline and every one-factor ablation on identical seeds.
# Usage (from python/):  bash scripts/run_ablations.sh
# Existing results are skipped; set FORCE=1 to recompute everything.
set -e
# one BLAS thread per worker: the process pool already uses every core
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

mc() {
  local tag=$1
  shift
  if [ -z "$FORCE" ] && [ -f "results/mc_$tag.json" ]; then echo "skip: $tag"; return; fi
  python scripts/monte_carlo.py --tag "$tag" "$@" > /dev/null
  echo "done: $tag"
}

mc S_final    --n 1000 --f-ctrl 100
mc S_order3   --n 500 --f-ctrl 100 --set order=3 sigma_top=8.0 tau_top=1.5
mc S_order2   --n 500 --f-ctrl 100 --set order=2
mc S_noreid   --n 500 --f-ctrl 100 --set use_reid=False
mc S_ambush   --n 500 --f-ctrl 100 --set use_ambush=True
mc S_pdev10   --n 500 --f-ctrl 100 --set p_dev=10.0
mc S_pdev30   --n 500 --f-ctrl 100 --set p_dev=30.0
mc S_sigma40  --n 500 --f-ctrl 100 --set sigma_top=40.0
mc S_pmin075  --n 500 --f-ctrl 100 --set p_min=0.75
mc S_psure0999 --n 500 --f-ctrl 100 --set p_sure=0.999
mc S_rate20   --n 500 --f-ctrl 20
mc S_theta60  --n 500 --f-ctrl 100 --theta0-deg 60
