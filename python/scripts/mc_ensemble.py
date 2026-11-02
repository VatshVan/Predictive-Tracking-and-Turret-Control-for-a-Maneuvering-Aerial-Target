"""Logged Monte Carlo ensemble for time-series statistics.

python scripts/mc_ensemble.py --n 300 --f-ctrl 100
Saves results/ensemble.npz with per-run histories on the controller grid:
  err   cross-range pointing error, turret vs true target [m]
  est   cross-range estimation error, nowcast vs true target [m]
  sig   filter-predicted cross-range std of the nowcast [m]
  phit  predicted hit probability
  mode  0 acquire, 1 track, 2 coast, 3 ambush
plus blackout intervals, pulse records and path features.
"""
import argparse
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from guardian.sim import run  # noqa: E402
from guardian.target_motion_generator import blackout_intervals  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), '..', 'results')
MODES = {'acquire': 0, 'track': 1, 'coast': 2, 'ambush': 3}


def _one(job):
    seed, f_ctrl = job
    s = f'MC{seed}'
    o = run(s, f_ctrl=f_ctrl)
    err = o['r_true'] * np.sin(o['th'] - o['th_true'])
    est = o['r_true'] * np.sin(o['th_hat'] - o['th_true'])
    mode = np.array([MODES.get(m, 0) for m in o['mode']], dtype=np.int8)
    pulses = [(p['t'], p['miss']) for p in o['pulses']]
    return (err.astype(np.float32), est.astype(np.float32), o['sig'].astype(np.float32),
            o['p_hit'].astype(np.float32), mode, blackout_intervals(s), pulses,
            o['metrics']['t_first_hit'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=300)
    ap.add_argument('--f-ctrl', type=float, default=100.0)
    a = ap.parse_args()
    with Pool(max(1, (os.cpu_count() or 2) - 1)) as pool:
        out = pool.map(_one, [(s, a.f_ctrl) for s in range(a.n)], chunksize=4)
    t = np.arange(len(out[0][0])) / a.f_ctrl
    bo = np.full((a.n, 2, 2), np.nan)
    for i, o in enumerate(out):
        for j, (t_on, dur) in enumerate(o[5][:2]):
            bo[i, j] = (t_on, dur)
    pulses = np.array([(i, t_p, miss) for i, o in enumerate(out) for (t_p, miss) in o[6]])
    os.makedirs(RESULTS, exist_ok=True)
    np.savez_compressed(os.path.join(RESULTS, 'ensemble.npz'), t=t,
                        err=np.stack([o[0] for o in out]), est=np.stack([o[1] for o in out]),
                        sig=np.stack([o[2] for o in out]), phit=np.stack([o[3] for o in out]),
                        mode=np.stack([o[4] for o in out]), blackouts=bo, pulses=pulses,
                        t_first_hit=np.array([o[7] for o in out]), f_ctrl=a.f_ctrl)
    print(f'saved ensemble of {a.n} runs, {len(t)} samples each')


if __name__ == '__main__':
    main()
