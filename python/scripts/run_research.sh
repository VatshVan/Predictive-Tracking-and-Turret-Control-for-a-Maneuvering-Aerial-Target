#!/bin/sh
# Research campaign (docs/research_lemma_homotopy.md): 1000 independent seeds per family.
# MC0-999 give the failure bound; HG0-999 are worst-case evidence, reported separately.
set -e
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=results/research
EX4="ref_mode=exact_taylor ref_n=4 order=5 sigma_top=150.0 tau_top=0.2"
run() { [ -f $OUT/mc_$2.json ] || python scripts/monte_carlo.py --n 1000 --prefix $1 --out $OUT --tag $2 $3 $4 $5 $6 $7 $8; }
python scripts/build_hg_cache.py --n 1000
for P in MC HG; do
  run $P ${P}_final
  run $P ${P}_exact3 --set ref_mode=exact_taylor ref_n=3
  run $P ${P}_exact4 --set $EX4
  run $P ${P}_mm --set ref_mode=moment_matched
  for A in 0.05 0.1 0.2 0.3 0.5; do
    run $P ${P}_blend$A --set ref_mode=blend blend_a=$A
  done
  run $P ${P}_homotopy --set phit_mode=homotopy
  run $P ${P}_noreid --set use_reid=False
done
run MC MC_final_pert --perturb-plant
run MC MC_mm_pert --perturb-plant --set ref_mode=moment_matched
