"""Generate the HG<seed> adversarial engagements once (closed-loop search, guardian/homotopy.py).

python scripts/build_hg_cache.py --n 1000
Writes results/hg_cache/HG<seed>.npz and results/research/hg_cache_summary.json.
"""
import argparse
import json
import os
import sys
from collections import Counter
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from guardian.homotopy import CLASSES, hg_def  # noqa: E402


def _one(job):
    seed, target = job
    h = hg_def(seed, target)['hg']
    h['seed'] = seed
    return h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=1000)
    ap.add_argument('--target', default='final')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    with Pool(a.workers) as pool:
        rows = pool.map(_one, [(s, a.target) for s in range(a.n)], chunksize=2)
    out = os.path.join(os.path.dirname(__file__), '..', 'results', 'research')
    os.makedirs(out, exist_ok=True)
    summ = dict(n=len(rows),
                realised=Counter(CLASSES[r['realised']] if r['realised'] >= 0 else 'all_through' for r in rows),
                intended_realised=sum(r['realised'] == r['cls'] for r in rows),
                lead=Counter(r['lead'] for r in rows),
                miss_mean=sum(r['miss_mean'] for r in rows) / len(rows))
    name = 'hg_cache_summary.json' if a.target == 'final' else f'hg_cache_summary_{a.target}.json'
    with open(os.path.join(out, name), 'w') as f:
        json.dump(dict(summary=summ, rows=rows), f, indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))


if __name__ == '__main__':
    main()
