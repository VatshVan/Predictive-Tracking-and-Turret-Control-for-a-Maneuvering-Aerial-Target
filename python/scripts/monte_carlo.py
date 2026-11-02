"""Monte Carlo evaluation over random admissible engagements ('MC<seed>' scenarios).

python scripts/monte_carlo.py --n 300 --f-ctrl 100 --tag baseline
python scripts/monte_carlo.py --n 300 --set use_reid=False --tag no_reid
"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from guardian import turret_controller as tc  # noqa: E402
from guardian.sim import run  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), '..', 'results')


def _parse_value(v):
    if v in ('True', 'False'):
        return v == 'True'
    try:
        return int(v)
    except ValueError:
        return float(v)


def _one(job):
    seed, f_ctrl, overrides, theta0 = job
    tc.CFG.update(overrides)
    o = run(f'MC{seed}', f_ctrl=f_ctrl, theta0=theta0)
    m = o['metrics']
    m['seed'] = seed
    m['pulses'] = [(p['t'], p['miss']) for p in o['pulses']]
    return m


def summarize(rows):
    succ = np.array([r['success'] for r in rows])
    fh = np.array([r['t_first_hit'] for r in rows])
    hits = sum(r['n_hits'] for r in rows)
    pulses = sum(r['n_pulses'] for r in rows)
    return dict(
        n=len(rows),
        success_rate=float(succ.mean()),
        t_first_hit_mean=float(np.nanmean(fh)),
        t_first_hit_p95=float(np.nanpercentile(fh, 95)),
        t_first_hit_max=float(np.nanmax(fh)),
        pulse_hit_rate=hits / max(pulses, 1),
        pulses_total=pulses,
        point_rms_median=float(np.median([r['point_rms_m'] for r in rows])),
        point_rms_p95=float(np.percentile([r['point_rms_m'] for r in rows], 95)),
        frac_in_Reng_mean=float(np.mean([r['frac_in_Reng'] for r in rows])),
        V_peak_max=float(np.max([r['V_peak'] for r in rows])),
        failed_seeds=[int(r['seed']) for r in rows if not r['success']],
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=200)
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--f-ctrl', type=float, default=100.0)
    ap.add_argument('--theta0-deg', type=float, default=0.0)
    ap.add_argument('--set', nargs='*', default=[], help='controller overrides key=value')
    ap.add_argument('--tag', default='baseline')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    overrides = {k: _parse_value(v) for k, v in (s.split('=', 1) for s in a.set)}
    jobs = [(s, a.f_ctrl, overrides, np.deg2rad(a.theta0_deg)) for s in range(a.start, a.start + a.n)]
    t0 = time.time()
    with Pool(a.workers) as pool:
        rows = pool.map(_one, jobs, chunksize=4)
    summ = summarize(rows)
    summ.update(tag=a.tag, f_ctrl=a.f_ctrl, overrides=overrides, wall_s=round(time.time() - t0, 1))
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, f'mc_{a.tag}.json'), 'w') as f:
        json.dump(dict(summary=summ, runs=rows), f, indent=1, default=float)
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()
