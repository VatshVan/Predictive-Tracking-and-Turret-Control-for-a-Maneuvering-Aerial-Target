#!/bin/sh
# Second campaign stage: contact-order pole-placed reference modes (research_pole_shift.py).
set -e
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=results/research
run() { [ -f $OUT/mc_$2.json ] || python scripts/monte_carlo.py --n 1000 --prefix $1 --out $OUT --tag $2 $3 $4 $5 $6; }
for P in MC HG; do
  run $P ${P}_contact1 --set ref_mode=contact
  run $P ${P}_contact2 --set ref_mode=contact "ref_poly=(1,0.5,0.125,0.015625,0.0009765625)"
done
