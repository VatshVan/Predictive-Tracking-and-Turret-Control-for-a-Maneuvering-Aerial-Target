"""Noise ensembles on the standardized and stress scenarios, paired across design variants.

python scripts/stress_ensemble.py --reps 50
Each scenario (A, B, C, D1-D6) is re-run with `reps` independent radar-noise realisations
('D6#k'), identical across variants, so differences are paired. Writes
results/stress_<tag>.json and results/stress_report.md.
"""
import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from guardian import turret_controller as tc  # noqa: E402
from guardian.sim import run  # noqa: E402
from statistics import boot_ci, holm, wilson  # noqa: E402

RESULTS = os.path.join(os.path.dirname(__file__), '..', 'results')
SCENARIOS = ['A', 'B', 'C', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6']
VARIANTS = {
    'final': ({}, 100.0),
    'order3': (dict(order=3, sigma_top=8.0, tau_top=1.5), 100.0),
    'psure0999': (dict(p_sure=0.999), 100.0),
    'pmin075': (dict(p_min=0.75), 100.0),
    'final_20Hz': ({}, 20.0),
    'order3_20Hz': (dict(order=3, sigma_top=8.0, tau_top=1.5), 20.0),
}
PAIRS = [('order3', 'final'), ('psure0999', 'final'), ('pmin075', 'final'), ('order3_20Hz', 'final_20Hz')]
T_END = 10.0
DEFAULT_CFG = dict(tc.CFG)


def _one(job):
    scenario, rep, overrides, f = job
    # pool workers are reused across variants: start every job from the default design
    tc.CFG.clear()
    tc.CFG.update(DEFAULT_CFG)
    tc.CFG.update(overrides)
    m = run(f'{scenario}#{rep}', f_ctrl=f)['metrics']
    return dict(scenario=scenario, rep=rep, success=bool(m['success']),
                t_first_hit=m['t_first_hit'] if m['success'] else T_END,
                n_pulses=m['n_pulses'], n_hits=m['n_hits'], point_rms_m=m['point_rms_m'],
                frac_in_Reng=m['frac_in_Reng'])


def run_variant(tag, reps, pool):
    path = os.path.join(RESULTS, f'stress_{tag}.json')
    if os.path.exists(path):
        with open(path) as f:
            rows = json.load(f)
        if len(rows) == len(SCENARIOS) * reps:
            return rows
    overrides, f_ctrl = VARIANTS[tag]
    jobs = [(s, r, overrides, f_ctrl) for s in SCENARIOS for r in range(reps)]
    rows = pool.map(_one, jobs, chunksize=4)
    with open(path, 'w') as fh:
        json.dump(rows, fh)
    return rows


def summary(rows):
    pul = sum(r['n_pulses'] for r in rows)
    hits = sum(r['n_hits'] for r in rows)
    succ = sum(r['success'] for r in rows)
    fh = np.array([r['t_first_hit'] for r in rows])
    rms = np.array([r['point_rms_m'] for r in rows])
    return dict(n=len(rows), success=succ, success_ci=wilson(succ, len(rows)), pulses=pul, misses=pul - hits,
                hit_rate=hits / pul, hit_ci=wilson(hits, pul), fh_mean=float(fh.mean()),
                fh_ci=boot_ci(fh, np.mean, 4000), rms_mean=float(rms.mean()), rms_ci=boot_ci(rms, np.mean, 4000))


def compare(var, base):
    key = {(r['scenario'], r['rep']): r for r in base}
    pr = [(key[(r['scenario'], r['rep'])], r) for r in var]
    d_fh = np.array([v['t_first_hit'] - b['t_first_hit'] for b, v in pr])
    d_rms = np.array([v['point_rms_m'] - b['point_rms_m'] for b, v in pr])
    d_miss = np.array([(v['n_pulses'] - v['n_hits']) - (b['n_pulses'] - b['n_hits']) for b, v in pr])

    def wil(d):
        return 1.0 if np.allclose(d, 0) else float(stats.wilcoxon(d, zero_method='pratt').pvalue)

    rng = np.random.default_rng(7)
    if np.any(d_miss):
        null = np.abs((rng.choice([-1, 1], size=(20000, len(d_miss))) * d_miss).sum(1))
        p_miss = float((np.sum(null >= abs(d_miss.sum())) + 1) / 20001)
    else:
        p_miss = 1.0
    per = {}
    for s in SCENARIOS:
        sel = [i for i, (b, _) in enumerate(pr) if b['scenario'] == s]
        per[s] = dict(d_rms=float(d_rms[sel].mean()), p_rms=wil(d_rms[sel]),
                      d_fh=float(d_fh[sel].mean()), p_fh=wil(d_fh[sel]),
                      miss_b=int(sum(pr[i][0]['n_pulses'] - pr[i][0]['n_hits'] for i in sel)),
                      miss_v=int(sum(pr[i][1]['n_pulses'] - pr[i][1]['n_hits'] for i in sel)))
    return dict(d_fh=float(d_fh.mean()), d_fh_ci=boot_ci(d_fh, np.mean, 4000), p_fh=wil(d_fh),
                d_rms=float(d_rms.mean()), d_rms_ci=boot_ci(d_rms, np.mean, 4000), p_rms=wil(d_rms),
                misses_base=int(sum(b['n_pulses'] - b['n_hits'] for b, _ in pr)),
                misses_var=int(sum(v['n_pulses'] - v['n_hits'] for _, v in pr)), p_miss=p_miss,
                per_scenario=per)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reps', type=int, default=50)
    a = ap.parse_args()
    with Pool(max(1, (os.cpu_count() or 2) - 1)) as pool:
        data = {t: run_variant(t, a.reps, pool) for t in VARIANTS}
    summ = {t: summary(r) for t, r in data.items()}
    comp = {f'{v} vs {b}': compare(data[v], data[b]) for v, b in PAIRS}
    for key in ('p_fh', 'p_rms', 'p_miss'):
        adj = holm({k: c[key] for k, c in comp.items()})
        for k in comp:
            comp[k][key + '_holm'] = adj[k]
    with open(os.path.join(RESULTS, 'stress_stats.json'), 'w') as f:
        json.dump(dict(summary=summ, compare=comp), f, indent=1)

    L = ['# Stress-scenario noise ensembles', '',
         f'Scenarios {", ".join(SCENARIOS)}; {a.reps} radar-noise realisations each, identical across variants.', '',
         '| Variant | runs | Intercepted | Pulse hit rate [95 % CI] | misses | First hit mean [s] [CI] | Pointing RMS mean [m] [CI] |',
         '|---|---|---|---|---|---|---|']
    for t, s in summ.items():
        L.append(f"| {t} | {s['n']} | {s['success']}/{s['n']} | {100 * s['hit_rate']:.2f} % "
                 f"[{100 * s['hit_ci'][0]:.2f}, {100 * s['hit_ci'][1]:.2f}] | {s['misses']} | {s['fh_mean']:.3f} "
                 f"[{s['fh_ci'][0]:.3f}, {s['fh_ci'][1]:.3f}] | {s['rms_mean']:.3f} [{s['rms_ci'][0]:.3f}, {s['rms_ci'][1]:.3f}] |")
    L += ['', '## Paired comparisons (variant minus reference, Holm-adjusted)', '',
          '| Comparison | Δ first hit [s] [CI] | p | Δ RMS [m] [CI] | p | misses ref → var | p |', '|---|---|---|---|---|---|---|']
    for k, c in comp.items():
        L.append(f"| {k} | {c['d_fh']:+.3f} [{c['d_fh_ci'][0]:+.3f}, {c['d_fh_ci'][1]:+.3f}] | {c['p_fh_holm']:.2g} | "
                 f"{c['d_rms']:+.3f} [{c['d_rms_ci'][0]:+.3f}, {c['d_rms_ci'][1]:+.3f}] | {c['p_rms_holm']:.2g} | "
                 f"{c['misses_base']} → {c['misses_var']} | {c['p_miss_holm']:.2g} |")
    L += ['', '## Per-scenario breakdown (unadjusted p)', '']
    for k, c in comp.items():
        L += [f'### {k}', '', '| Scenario | Δ RMS [m] | p | Δ first hit [s] | p | misses ref → var |', '|---|---|---|---|---|---|']
        for s, d in c['per_scenario'].items():
            L.append(f"| {s} | {d['d_rms']:+.3f} | {d['p_rms']:.2g} | {d['d_fh']:+.3f} | {d['p_fh']:.2g} | "
                     f"{d['miss_b']} → {d['miss_v']} |")
        L.append('')
    with open(os.path.join(RESULTS, 'stress_report.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')
    sys.stdout.reconfigure(encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
