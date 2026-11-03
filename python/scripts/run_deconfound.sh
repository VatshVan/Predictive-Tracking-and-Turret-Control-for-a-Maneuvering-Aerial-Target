#!/bin/sh
# De-confounding check: the HG adversary re-tuned against each candidate design (200 seeds).
set -e
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=results/research/deconfound
N=200
run() { [ -f $OUT/mc_$2.json ] || python scripts/monte_carlo.py --n $N --prefix $1 --out $OUT --tag $2 $3 $4 $5; }
for T in contact1 exact3; do
  python scripts/build_hg_cache.py --n $N --target $T
done
for T in contact1 exact3; do
  run HX${T}_ HX${T}_final
  run HX${T}_ HX${T}_contact1 --set ref_mode=contact
  run HX${T}_ HX${T}_exact3 --set ref_mode=exact_taylor ref_n=3
done
