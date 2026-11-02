"""Controller-rate study: performance and compute time versus controller call frequency.

python scripts/rate_study.py
The plant is always integrated at <= 1 ms; only the controller call rate changes.
"""
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from guardian import turret_controller as tc  # noqa: E402
from guardian.sim import run  # noqa: E402

SCENARIOS = ['A', 'B', 'C', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6'] + [f'MC{i}' for i in range(15)]
RATES = [20, 40, 60, 100, 200, 500, 1000, 2000, 5000]


def _one(job):
    scenario, f = job
    stats = {'n': 0, 'sec': 0.0, 'max': 0.0}

    def ctrl(*args):
        t0 = time.perf_counter()
        out = tc.turret_controller(*args)
        dt = time.perf_counter() - t0
        stats['n'] += 1
        stats['sec'] += dt
        stats['max'] = max(stats['max'], dt)
        return out

    tc.reset()
    m = run(scenario, f_ctrl=f, dt_plant=min(1e-3, 1.0 / f), controller=ctrl)['metrics']
    return f, m, stats


def main():
    jobs = [(s, f) for f in RATES for s in SCENARIOS]
    with Pool(max(1, (os.cpu_count() or 2) - 1)) as pool:
        results = pool.map(_one, jobs, chunksize=1)
    print('%6s %9s %10s %9s %9s %9s %10s %10s' % ('f [Hz]', 'success', 'hits', 'fh mean', 'fh max',
                                                   'rms med', 'call mean', 'call max'))
    for f in RATES:
        rows = [(m, st) for (ff, m, st) in results if ff == f]
        succ = sum(m['success'] for m, _ in rows)
        hits = sum(m['n_hits'] for m, _ in rows)
        pulses = sum(m['n_pulses'] for m, _ in rows)
        fh = [m['t_first_hit'] for m, _ in rows]
        rms = [m['point_rms_m'] for m, _ in rows]
        n_calls = sum(st['n'] for _, st in rows)
        mean_ms = 1e3 * sum(st['sec'] for _, st in rows) / n_calls
        max_ms = 1e3 * max(st['max'] for _, st in rows)
        print('%6d %5d/%-3d %5d/%-4d %8.2fs %8.2fs %8.2fm %8.3fms %8.2fms' % (
            f, succ, len(rows), hits, pulses, np.nanmean(fh), np.nanmax(fh), np.median(rms), mean_ms, max_ms))
    tc.reset()
    m = run('B', f_ctrl=100, integrator='euler', dt_plant=1e-3)['metrics']
    print('\nforward-Euler plant at 1 ms, controller 100 Hz, scenario B: success=%s, first hit %.2f s'
          % (m['success'], m['t_first_hit']))


if __name__ == '__main__':
    main()
