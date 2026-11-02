"""Robustness of the frozen design to how an unseen evaluation script might call it.

python scripts/interface_robustness.py --n 300
The instructor's script is not available, so each variant emulates a plausible convention:

  I1  reference: NaN + meas_valid=0 between radar ticks, 100 Hz
  I2  last measurement held (meas_valid=0) between ticks and during blackouts
  I3  last measurement repeated with meas_valid=1 on every call (stale duplicates)
  I4  irregular radar: a sample is taken at the call time once >= 50 ms have elapsed, 30 Hz calls
  I5  as I4 at 333 Hz (3 ms step, which does not divide 50 ms)
  I6  turret starts at a random angle U(-10, 85) deg
  I7  forward-Euler plant integration at 1 ms
  I8  reference convention at 1 kHz
Writes results/interface_robustness.json and results/interface_report.md.
"""
import argparse
import json
import math
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from guardian import params as P  # noqa: E402
from guardian.sim import run  # noqa: E402
from guardian.target_motion_generator import (blackout_intervals, target_motion_generator,  # noqa: E402
                                              target_truth)
from mc_statistics import boot_ci, wilson  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), '..', 'results')
VARIANTS = {
    'I1': ('Reference: NaN + valid=0 between ticks, 100 Hz', 100.0),
    'I2': ('Held last value, valid=0 between ticks', 100.0),
    'I3': ('Stale value repeated with valid=1 every call', 100.0),
    'I4': ('Irregular sampling at call time, 30 Hz calls', 30.0),
    'I5': ('Irregular sampling at call time, 333 Hz calls', 1000.0 / 3.0),
    'I6': ('Random initial turret angle U(-10, 85) deg', 100.0),
    'I7': ('Forward-Euler plant at 1 ms', 100.0),
    'I8': ('Reference convention at 1 kHz', 1000.0),
}


def _blocked(t, bos):
    return any(t_on - 1e-9 <= t < t_on + dur - 1e-9 for t_on, dur in bos)


def make_generator(variant, seed):
    """Return a measurement source emulating one evaluation-script convention."""
    if variant in ('I1', 'I6', 'I7', 'I8'):
        return None
    if variant in ('I2', 'I3'):
        last = {'m': (float('nan'), float('nan'))}

        def gen(t, scenario):
            r, th, valid = target_motion_generator(t, scenario)
            if valid:
                last['m'] = (r, th)
                return r, th, 1
            return last['m'][0], last['m'][1], (1 if variant == 'I3' and last['m'][0] == last['m'][0] else 0)
        return gen
    rng = np.random.default_rng(seed + 555)
    state = {'t_last': -1e9}

    def gen(t, scenario):
        if t - state['t_last'] < P.DT_RADAR - 1e-9:
            return float('nan'), float('nan'), 0
        if _blocked(t, blackout_intervals(scenario)):
            return float('nan'), float('nan'), 0
        state['t_last'] = t
        x, y = target_truth(t, scenario)[:2]
        return (math.hypot(x, y) + rng.normal() * P.SIGMA_R,
                math.atan2(y, x) + rng.normal() * P.SIGMA_TH, 1)
    return gen


def _one(job):
    variant, seed = job
    _, f_ctrl = VARIANTS[variant]
    theta0 = 0.0
    if variant == 'I6':
        theta0 = math.radians(np.random.default_rng(seed + 999).uniform(-10.0, 85.0))
    o = run(f'MC{seed}', f_ctrl=f_ctrl, theta0=theta0, generator=make_generator(variant, seed),
            integrator='euler' if variant == 'I7' else 'rk4')
    m = o['metrics']
    return dict(variant=variant, seed=seed, success=bool(m['success']),
                t_first_hit=m['t_first_hit'] if m['success'] else P.T_END,
                n_pulses=m['n_pulses'], n_hits=m['n_hits'], point_rms_m=m['point_rms_m'],
                theta0_deg=math.degrees(theta0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=300)
    a = ap.parse_args()
    jobs = [(v, s) for v in VARIANTS for s in range(a.n)]
    with Pool(max(1, (os.cpu_count() or 2) - 1)) as pool:
        rows = pool.map(_one, jobs, chunksize=4)
    by = {v: sorted([r for r in rows if r['variant'] == v], key=lambda r: r['seed']) for v in VARIANTS}
    ref = by['I1']
    out = {}
    for v, rs in by.items():
        succ = sum(r['success'] for r in rs)
        pul = sum(r['n_pulses'] for r in rs)
        hits = sum(r['n_hits'] for r in rs)
        fh = np.array([r['t_first_hit'] for r in rs])
        rms = np.array([r['point_rms_m'] for r in rs])
        d_fh = fh - np.array([r['t_first_hit'] for r in ref])
        d_rms = rms - np.array([r['point_rms_m'] for r in ref])

        def wil(d):
            return 1.0 if np.allclose(d, 0) else float(stats.wilcoxon(d, zero_method='pratt').pvalue)

        out[v] = dict(label=VARIANTS[v][0], n=len(rs), success=succ, success_ci=wilson(succ, len(rs)),
                      pulses=pul, misses=pul - hits, hit_rate=hits / pul, hit_ci=wilson(hits, pul),
                      fh_mean=float(fh.mean()), fh_ci=boot_ci(fh, np.mean, 4000), fh_max=float(fh.max()),
                      rms_median=float(np.median(rms)),
                      d_fh=float(d_fh.mean()), p_fh=wil(d_fh), d_rms=float(d_rms.mean()), p_rms=wil(d_rms))
    if 'I6' in by:
        th = np.array([r['theta0_deg'] for r in by['I6']])
        out['I6']['theta0_range'] = [float(th.min()), float(th.max())]
    with open(os.path.join(RESULTS, 'interface_robustness.json'), 'w') as f:
        json.dump(dict(summary=out, runs=rows), f, indent=1)
    L = ['# Evaluation-interface robustness', '',
         f'{a.n} random admissible engagements per variant (seeds 0-{a.n - 1}), frozen design.', '',
         '| ID | Convention | Intercepted | Pulse hit rate [95 % CI] | misses | First hit mean [CI] / max [s] | '
         'Pointing RMS median [m] | Δ first hit vs I1, p | Δ RMS vs I1, p |',
         '|---|---|---|---|---|---|---|---|---|']
    for v, s in out.items():
        L.append(f"| {v} | {s['label']} | {s['success']}/{s['n']} | {100 * s['hit_rate']:.2f} % "
                 f"[{100 * s['hit_ci'][0]:.2f}, {100 * s['hit_ci'][1]:.2f}] | {s['misses']} | "
                 f"{s['fh_mean']:.3f} [{s['fh_ci'][0]:.3f}, {s['fh_ci'][1]:.3f}] / {s['fh_max']:.2f} | "
                 f"{s['rms_median']:.3f} | {s['d_fh']:+.3f}, {s['p_fh']:.2g} | {s['d_rms']:+.3f}, {s['p_rms']:.2g} |")
    with open(os.path.join(RESULTS, 'interface_report.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')
    sys.stdout.reconfigure(encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
