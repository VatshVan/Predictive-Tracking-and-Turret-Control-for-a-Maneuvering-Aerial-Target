"""Controller-rate study: performance and compute time versus controller call frequency.

python scripts/rate_study.py
The plant is always integrated at <= 1 ms; only the controller call rate changes.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from guardian import turret_controller as tc  # noqa: E402
from guardian.sim import run  # noqa: E402

SCENARIOS = ['A', 'B', 'C', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6'] + [f'MC{i}' for i in range(15)]
RATES = [20, 40, 60, 100, 200, 500, 1000, 2000, 5000]


def timed_controller():
    stats = {'n': 0, 'sec': 0.0, 'max': 0.0}

    def ctrl(*args):
        t0 = time.perf_counter()
        out = tc.turret_controller(*args)
        dt = time.perf_counter() - t0
        stats['n'] += 1
        stats['sec'] += dt
        stats['max'] = max(stats['max'], dt)
        return out
    return ctrl, stats


def main():
    print('%6s %9s %8s %9s %9s %9s %10s %10s' % ('f [Hz]', 'success', 'hits', 'fh mean', 'fh max',
                                                  'rms med', 'call mean', 'call max'))
    for f in RATES:
        succ, hits, pulses, fh, rms = 0, 0, 0, [], []
        ctrl, st = timed_controller()
        dt_plant = min(1e-3, 1.0 / f)
        for s in SCENARIOS:
            tc.reset()
            m = run(s, f_ctrl=f, dt_plant=dt_plant, controller=ctrl)['metrics']
            succ += m['success']
            hits += m['n_hits']
            pulses += m['n_pulses']
            fh.append(m['t_first_hit'])
            rms.append(m['point_rms_m'])
        print('%6d %5d/%-3d %4d/%-4d %8.2fs %8.2fs %8.2fm %8.3fms %8.2fms' % (
            f, succ, len(SCENARIOS), hits, pulses, np.nanmean(fh), np.nanmax(fh), np.median(rms),
            1e3 * st['sec'] / st['n'], 1e3 * st['max']))
    tc.reset()
    m = run('B', f_ctrl=100, integrator='euler', dt_plant=1e-3)['metrics']
    print('\nforward-Euler plant at 1 ms, controller 100 Hz, scenario B: success=%s first hit %.2f s'
          % (m['success'], m['t_first_hit']))


if __name__ == '__main__':
    main()
