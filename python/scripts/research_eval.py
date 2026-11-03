"""Statistics of the research campaign (scripts/run_research.sh -> results/research/mc_*.json).

python scripts/research_eval.py
Writes results/research/research_stats.json and results/research/research_report.md.

Families are kept apart: MC0-999 (random admissible engagements, the only source of the
failure bound) and HG0-999 (guided adversarial breaks, worst-case evidence, never pooled).
Each variant is paired with the final design of its own family on identical seeds, using the
tests of mc_statistics.py (Wilcoxon, sign-flip permutation, exact McNemar, Holm per family).
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from mc_statistics import describe, fmt_ci, holm, paired  # noqa: E402

RES = os.path.join(os.path.dirname(__file__), '..', 'results', 'research')
LABELS = {
    'final': 'Final design (deviation feedback, mixture P_hit)',
    'exact3': 'Exact Taylor n = 3 (T = 1/A_3)',
    'exact4': 'Exact Taylor n = 4 (T = 1/A_4, order-5 estimator)',
    'mm': 'Moment-matched (T = 1/Q, Q = A_3 + q4 x^4 + q5 x^5)',
    'blend0.05': 'Blend a = 0.05 (P = 0.95 A_3 + 0.05 D_dev)',
    'blend0.1': 'Blend a = 0.1',
    'blend0.2': 'Blend a = 0.2',
    'blend0.3': 'Blend a = 0.3',
    'blend0.5': 'Blend a = 0.5',
    'contact1': 'Contact order 1, P = (1 + s/6)^3 (decay bound)',
    'contact2': 'Contact order 2, P = ((s^2 + 8s + 32)/32)^2',
    'homotopy': 'Homotopy event-signature fire planner',
    'noreid': 'No class re-identification',
    'final_pert': 'Final design, plant +-20 % (R, L, J_eq, b)',
    'final_pert6': 'Final design, plant +-20 % (R, L, J_eq, b, K_t, K_e)',
    'mm_pert': 'Moment-matched, plant +-20 %',
}
BLENDS = ('blend0.05', 'blend0.1', 'blend0.2', 'blend0.3', 'blend0.5', 'contact1', 'contact2')
PAIRS = {'MC': [(v, 'final') for v in ('exact3', 'exact4', 'mm') + BLENDS + ('homotopy', 'noreid', 'final_pert', 'final_pert6')]
               + [('mm_pert', 'final_pert')],
         'HG': [(v, 'final') for v in ('exact3', 'exact4', 'mm') + BLENDS + ('homotopy', 'noreid')]}
HG_CLASSES = ('above', 'below', 'cross_down', 'cross_up', 'multiple')


def load(tag):
    path = os.path.join(RES, f'mc_{tag}.json')
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return sorted(json.load(f)['runs'], key=lambda r: r['seed'])


def counters(runs):
    proj = np.array([r.get('n_proj', 0) for r in runs])
    wait = np.array([r.get('n_wait', 0) for r in runs])
    return dict(proj_mean=float(proj.mean()), proj_runs_frac=float(np.mean(proj > 0)), proj_max=int(proj.max()),
                wait_mean=float(wait.mean()), wait_runs_frac=float(np.mean(wait > 0)))


def by_class(runs):
    """HG breakdown by intended class (seed % 5) and H-signature of each pulse offset."""
    out = {}
    for c, name in enumerate(HG_CLASSES):
        sel = [r for r in runs if r['seed'] % 5 == c]
        offs = [p[4] for r in sel for p in r['pulses'] if len(p) > 4]
        sig = np.sign(np.where(np.abs(offs) > 2.0, offs, 0.0)) if offs else np.zeros(0)
        pul = sum(r['n_pulses'] for r in sel)
        hits = sum(r['n_hits'] for r in sel)
        out[name] = dict(n=len(sel), success=sum(r['success'] for r in sel), misses=pul - hits,
                         hit_rate=hits / max(pul, 1), rms_median=float(np.median([r['point_rms_m'] for r in sel])),
                         above=int(np.sum(sig > 0)), below=int(np.sum(sig < 0)),
                         miss_mean=float(np.mean([p[1] for r in sel for p in r['pulses']])))
    return out


def calib(runs, edges=(0.0, 0.95, 0.99, 0.999, 1.0001)):
    """Predicted P_hit at each pulse against the observed hit frequency (hit: miss <= R_eng = 2 m)."""
    ph = np.array([p[2] for r in runs for p in r['pulses']])
    hit = np.array([p[1] <= 2.0 for r in runs for p in r['pulses']], float)
    bins = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (ph >= lo) & (ph < hi)
        if k.any():
            bins.append(dict(lo=lo, hi=min(hi, 1.0), n=int(k.sum()), pred=float(ph[k].mean()), obs=float(hit[k].mean())))
    return dict(n=len(ph), pred=float(ph.mean()), obs=float(hit.mean()), bins=bins)


def main():
    stats = {}
    L = ['# Research campaign: exact-prediction modes and beam-relative homotopy', '',
         'Controller at 100 Hz, plant RK4 at 1 ms. Each family has 1000 independent seeds (0-999); no reused '
         'trajectories or noise replicas are pooled. HG results are worst-case evidence and are not used in '
         'the failure bound.', '']
    for fam in ('MC', 'HG'):
        tags = [t for t in LABELS if load(f'{fam}_{t}') is not None]
        if 'final' not in tags:
            continue
        runs = {t: load(f'{fam}_{t}') for t in tags}
        desc = {t: describe(runs[t]) for t in tags}
        cnt = {t: counters(runs[t]) for t in tags}
        comp = {}
        for v, b in PAIRS[fam]:
            if v in runs and b in runs:
                comp[f'{v} vs {b}'] = paired(runs[b], runs[v])
        for key in ('p_fh', 'p_rms', 'p_miss', 'p_success'):
            adj = holm({k: c[key] for k, c in comp.items()})
            for k in comp:
                comp[k][key + '_holm'] = adj[k]
        stats[fam] = dict(describe=desc, paired=comp, counters=cnt, calib={t: calib(runs[t]) for t in tags})
        if fam == 'HG':
            stats[fam]['by_class'] = {t: by_class(runs[t]) for t in tags}

        title = 'MC0-999: random admissible engagements' if fam == 'MC' else 'HG0-999: guided adversarial breaks'
        L += [f'## {title}', '',
              '| Configuration | Intercepted | Success CI (Wilson) | Pulse hit rate | misses / pulses | First hit mean [s] | '
              'p95 [s] | Pointing RMS median [m] | p95 [m] | Time in R_eng |',
              '|---|---|---|---|---|---|---|---|---|---|']
        for t in tags:
            d = desc[t]
            L.append(f"| {LABELS[t]} | {d['success']}/{d['n']} | {fmt_ci(d['success_ci'], '%.4f')} | "
                     f"{100 * d['pulse_hit_rate']:.2f} % | {d['misses']} / {d['pulses']} | {d['fh_mean']:.3f} | "
                     f"{d['fh_p95']:.2f} | {d['rms_median']:.3f} | {d['rms_p95']:.3f} | {100 * d['inR_mean']:.2f} % |")
        if fam == 'MC':
            d = desc['final']
            L += ['', f"Final design: {d['success']}/{d['n']} intercepted; one-sided 97.5 % upper bound on the failure "
                  f"probability {100 * d['failure_upper95']:.2f} % (Wilson, MC0-999 only)."]
        L += ['', '### Paired comparisons (variant minus reference, Holm-adjusted within this family)', '',
              '| Comparison | Δ first hit [s] [CI] | p | Δ pointing RMS [m] [CI] | p | runs worse | misses ref → var | p | '
              'McNemar lost/gained | p |',
              '|---|---|---|---|---|---|---|---|---|---|']
        for k, c in comp.items():
            L.append(f"| {k} | {c['d_fh_mean']:+.3f} {fmt_ci(c['d_fh_ci'])} | {c['p_fh_holm']:.2g} | "
                     f"{c['d_rms_mean']:+.3f} {fmt_ci(c['d_rms_ci'])} | {c['p_rms_holm']:.2g} | "
                     f"{100 * c['rms_worse_frac']:.0f} % | {c['misses_base']} → {c['misses_var']} | {c['p_miss_holm']:.2g} | "
                     f"{c['mcnemar_lost']}/{c['mcnemar_gained']} | {c['p_success_holm']:.2g} |")
        L += ['', '### Counters', '', '| Configuration | envelope projections / run | runs with a projection | max | '
              'look-ahead waits / run | runs with a wait |', '|---|---|---|---|---|---|']
        for t in tags:
            c = cnt[t]
            L.append(f"| {t} | {c['proj_mean']:.2f} | {100 * c['proj_runs_frac']:.1f} % | {c['proj_max']} | "
                     f"{c['wait_mean']:.2f} | {100 * c['wait_runs_frac']:.1f} % |")
        c = stats[fam]['calib']['final']
        L += ['', '### Hit-probability calibration at the pulses (final design)', '',
              f"All pulses: predicted {100 * c['pred']:.2f} %, observed {100 * c['obs']:.2f} % (n = {c['n']}).", '',
              '| P_hit bin | pulses | mean predicted | observed |', '|---|---|---|---|']
        L += [f"| [{b['lo']:.3f}, {b['hi']:.3f}) | {b['n']} | {b['pred']:.4f} | {b['obs']:.4f} |" for b in c['bins']]
        if fam == 'HG':
            L += ['', '### HG by homotopy class (seed % 5)', '', '| Configuration | ' + ' | '.join(HG_CLASSES) + ' |',
                  '|---|' + '---|' * len(HG_CLASSES)]
            for t in tags:
                bc = stats[fam]['by_class'][t]
                L.append(f'| {t} | ' + ' | '.join(
                    f"{bc[c]['success']}/{bc[c]['n']}, {bc[c]['misses']} miss (+{bc[c]['above']}/-{bc[c]['below']}), "
                    f"mean miss {bc[c]['miss_mean']:.2f} m"
                    for c in HG_CLASSES) + ' |')
        L.append('')
    with open(os.path.join(RES, 'research_stats.json'), 'w') as f:
        json.dump(stats, f, indent=1)
    with open(os.path.join(RES, 'research_report.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')
    sys.stdout.reconfigure(encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
