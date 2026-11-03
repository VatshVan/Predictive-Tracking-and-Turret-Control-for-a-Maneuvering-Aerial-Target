"""Monte Carlo evaluation over random admissible engagements ('MC<seed>' scenarios).

python scripts/monte_carlo.py --n 300 --f-ctrl 100 --tag baseline
python scripts/monte_carlo.py --n 300 --set use_reid=False --tag no_reid
python scripts/monte_carlo.py --n 1000 --prefix HG --set ref_mode=exact_taylor --out results/research --tag HG_exact3
--perturb-plant draws R, L, J_eq and b independently in +-20 % per seed (controller stays nominal);
--perturb-params R,L,J_EQ,B,K_T,K_E also perturbs the motor constants.
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
from guardian import turret_plant as tp  # noqa: E402
from guardian.sim import run  # noqa: E402
from guardian.target_motion_generator import _get  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), '..', 'results')


def _parse_value(v):
    if v in ('True', 'False'):
        return v == 'True'
    if v.startswith('(') and v.endswith(')'):
        return tuple(float(x) for x in v[1:-1].split(','))
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v


PLANT_NOM = dict(R=tp.R, L=tp.L, J_EQ=tp.J_EQ, B=tp.B, K_T=tp.K_T, K_E=tp.K_E)
PERTURB_DEFAULT = ('R', 'L', 'J_EQ', 'B')


def _perturb_plant(seed, on, names=PERTURB_DEFAULT):
    """Draw every factor in a fixed order (so a seed's R, L, J_EQ, B draws do not depend on
    which parameters are perturbed) and apply those in `names`."""
    rng = np.random.default_rng(770000 + seed)
    draws = {k: rng.uniform(0.8, 1.2) for k in PLANT_NOM}
    f = {k: (draws[k] if on and k in names else 1.0) for k in PLANT_NOM}
    for k, v in PLANT_NOM.items():
        setattr(tp, k, v * f[k])
    return f


def _one(job):
    seed, f_ctrl, overrides, theta0, prefix, perturb = job[:6]
    names = job[6] if len(job) > 6 else PERTURB_DEFAULT
    tc.CFG.update(overrides)
    scenario = f'{prefix}{seed}'
    factors = _perturb_plant(seed, perturb, names)
    o = run(scenario, f_ctrl=f_ctrl, theta0=theta0)
    m = o['metrics']
    m['seed'] = seed
    dbg = tc.get_debug() or {}
    m['n_proj'] = int(dbg.get('n_proj', 0))
    m['n_wait'] = int(dbg.get('n_wait', 0))
    if perturb:
        m['plant_factors'] = factors
    pulses = []
    for p in o['pulses']:
        k = int(round(p['t'] * f_ctrl))
        pulses.append((p['t'], p['miss'], float(o['p_hit'][k]), str(o['mode'][k]), p['offset']))
    m['pulses'] = pulses
    m['features'] = path_features(scenario, theta0)
    return m


def path_features(scenario, theta0=0.0):
    """Descriptors of how a random engagement was generated (for the experiment report)."""
    e = _get(scenario)
    tr = e['traj']
    k = tr['t'] <= 10.0 + 1e-9
    x, y, vx, vy, ax, ay = (tr[n][k] for n in ('x', 'y', 'vx', 'vy', 'ax', 'ay'))
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    los_rate = (x * vy - y * vx) / r ** 2
    jerk = np.hypot(np.diff(ax), np.diff(ay)) / (tr['t'][1] - tr['t'][0])
    bo = e['blackouts']
    return dict(
        r0=float(r[0]), th0_deg=float(np.degrees(th[0])), speed0=float(np.hypot(vx[0], vy[0])),
        r_min=float(r.min()), r_max=float(r.max()),
        v_max=float(np.hypot(vx, vy).max()), a_max=float(np.hypot(ax, ay).max()),
        a_rms=float(np.sqrt(np.mean(ax ** 2 + ay ** 2))), jerk_max=float(jerk.max()),
        los_rate_max=float(np.abs(los_rate).max()),
        th_min_deg=float(np.degrees(th.min())), th_max_deg=float(np.degrees(th.max())),
        n_blackouts=len(bo), blackouts=[(float(a), float(b)) for a, b in bo],
        blackout_total=float(sum(b for _, b in bo)),
        slew_deg=float(abs(np.degrees(th[0] - theta0))),
    )


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
    ap.add_argument('--prefix', default='MC', help="scenario family: 'MC' or 'HG'")
    ap.add_argument('--out', default=RESULTS)
    ap.add_argument('--perturb-plant', action='store_true')
    ap.add_argument('--perturb-params', default=','.join(PERTURB_DEFAULT),
                    help='plant parameters drawn in +-20 %% (from R,L,J_EQ,B,K_T,K_E)')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    overrides = {k: _parse_value(v) for k, v in (s.split('=', 1) for s in a.set)}
    names = tuple(a.perturb_params.split(','))
    jobs = [(s, a.f_ctrl, overrides, np.deg2rad(a.theta0_deg), a.prefix, a.perturb_plant, names)
            for s in range(a.start, a.start + a.n)]
    t0 = time.time()
    with Pool(a.workers) as pool:
        rows = pool.map(_one, jobs, chunksize=4)
    summ = summarize(rows)
    summ.update(tag=a.tag, f_ctrl=a.f_ctrl, overrides=overrides, prefix=a.prefix,
                perturb_plant=a.perturb_plant, perturb_params=list(names) if a.perturb_plant else [], wall_s=round(time.time() - t0, 1))
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, f'mc_{a.tag}.json'), 'w') as f:
        json.dump(dict(summary=summ, runs=rows), f, indent=1, default=float)
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()
