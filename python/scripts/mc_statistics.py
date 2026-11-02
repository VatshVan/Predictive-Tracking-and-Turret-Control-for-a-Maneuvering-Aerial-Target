"""Statistical analysis of the Monte Carlo campaign (results/mc_S_*.json).

python scripts/mc_statistics.py
Writes results/stats.json and results/stats_report.md.

Per configuration
  interception success      k/n with Wilson 95 % CI
  pulse hit rate            pooled, with a cluster (run-level) bootstrap 95 % CI, since pulses
                            within one engagement are not independent
  first-hit time            mean / median / p95 / max, bootstrap 95 % CIs for mean and p95
  pointing RMS, time in R   bootstrap 95 % CIs
Paired comparison of each variant against the final baseline on identical seeds
  first-hit time, pointing RMS   Wilcoxon signed-rank (Pratt zero handling) + bootstrap CI of
                                 the mean paired difference
  misses per run                 two-sided sign-flip permutation test (20 000 permutations)
  interception success           exact McNemar (binomial test on discordant pairs)
  multiplicity                   Holm-Bonferroni within each metric family
A run without any hit is censored at t_end = 10 s for paired first-hit comparisons.
"""
import glob
import json
import math
import os
import sys

import numpy as np
from scipy import stats

RESULTS = os.path.join(os.path.dirname(__file__), '..', 'results')
BASE = 'S_final'
N_BOOT = 10000
N_PERM = 20000
RNG = np.random.default_rng(2044)
T_END = 10.0

LABELS = {
    'S_final': 'Final design',
    'S_order3': 'Order 3 (Singer acceleration)',
    'S_order2': 'Order 2 (constant velocity)',
    'S_noreid': 'No class re-identification',
    'S_ambush': 'Ambush hold enabled',
    'S_pdev10': 'Deviation pole p = 10',
    'S_pdev30': 'Deviation pole p = 30',
    'S_sigma40': 'Jerk noise sigma_j = 40',
    'S_pmin075': 'Coast threshold 0.75',
    'S_psure0999': 'Instant-fire threshold 0.999',
    'S_rate20': 'Controller at 20 Hz',
    'S_theta60': 'Turret starts at 60 deg',
}


def wilson(k, n, z=1.959964):
    if n == 0:
        return (float('nan'), float('nan'))
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))


def boot_ci(x, fn, n=N_BOOT):
    x = np.asarray(x)
    idx = RNG.integers(0, len(x), (n, len(x)))
    vals = np.array([fn(x[i]) for i in idx])
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def cluster_hit_rate_ci(hits, pulses, n=N_BOOT):
    hits, pulses = np.asarray(hits), np.asarray(pulses)
    idx = RNG.integers(0, len(hits), (n, len(hits)))
    rates = hits[idx].sum(1) / np.maximum(pulses[idx].sum(1), 1)
    return (float(np.percentile(rates, 2.5)), float(np.percentile(rates, 97.5)))


def holm(pvals):
    """Holm-Bonferroni adjusted p-values (dict name -> p)."""
    names = sorted(pvals, key=lambda k: pvals[k])
    m = len(names)
    adj, running = {}, 0.0
    for i, k in enumerate(names):
        running = max(running, min(1.0, (m - i) * pvals[k]))
        adj[k] = running
    return adj


def load(tag):
    with open(os.path.join(RESULTS, f'mc_{tag}.json')) as f:
        d = json.load(f)
    runs = sorted(d['runs'], key=lambda r: r['seed'])
    return d['summary'], runs


def describe(runs):
    n = len(runs)
    succ = np.array([r['success'] for r in runs])
    fh = np.array([r['t_first_hit'] for r in runs], dtype=float)
    fh_ok = fh[np.isfinite(fh)]
    hits = np.array([r['n_hits'] for r in runs])
    pul = np.array([r['n_pulses'] for r in runs])
    rms = np.array([r['point_rms_m'] for r in runs])
    inr = np.array([r['frac_in_Reng'] for r in runs])
    return dict(
        n=n,
        success=int(succ.sum()), success_rate=float(succ.mean()), success_ci=wilson(int(succ.sum()), n),
        failure_upper95=wilson(int((~succ).sum()), n)[1],
        pulses=int(pul.sum()), hits=int(hits.sum()), misses=int(pul.sum() - hits.sum()),
        pulse_hit_rate=float(hits.sum() / pul.sum()),
        pulse_hit_ci_wilson=wilson(int(hits.sum()), int(pul.sum())),
        pulse_hit_ci_cluster=cluster_hit_rate_ci(hits, pul),
        fh_mean=float(fh_ok.mean()), fh_mean_ci=boot_ci(fh_ok, np.mean),
        fh_median=float(np.median(fh_ok)), fh_p95=float(np.percentile(fh_ok, 95)),
        fh_p95_ci=boot_ci(fh_ok, lambda v: np.percentile(v, 95)), fh_max=float(fh_ok.max()),
        rms_median=float(np.median(rms)), rms_median_ci=boot_ci(rms, np.median),
        rms_p95=float(np.percentile(rms, 95)),
        inR_mean=float(inr.mean()), inR_mean_ci=boot_ci(inr, np.mean),
        V_peak=float(max(r['V_peak'] for r in runs)),
    )


def paired(base_runs, var_runs):
    b = {r['seed']: r for r in base_runs}
    pairs = [(b[r['seed']], r) for r in var_runs if r['seed'] in b]
    fb = np.array([p[0]['t_first_hit'] for p in pairs], dtype=float)
    fv = np.array([p[1]['t_first_hit'] for p in pairs], dtype=float)
    fb[~np.isfinite(fb)] = T_END
    fv[~np.isfinite(fv)] = T_END
    rb = np.array([p[0]['point_rms_m'] for p in pairs])
    rv = np.array([p[1]['point_rms_m'] for p in pairs])
    mb = np.array([p[0]['n_pulses'] - p[0]['n_hits'] for p in pairs])
    mv = np.array([p[1]['n_pulses'] - p[1]['n_hits'] for p in pairs])
    sb = np.array([p[0]['success'] for p in pairs])
    sv = np.array([p[1]['success'] for p in pairs])

    def wil(d):
        if np.allclose(d, 0):
            return 1.0
        return float(stats.wilcoxon(d, zero_method='pratt').pvalue)

    d_fh, d_rms, d_miss = fv - fb, rv - rb, mv - mb
    if np.any(d_miss != 0):
        signs = RNG.choice([-1, 1], size=(N_PERM, len(d_miss)))
        null = np.abs((signs * d_miss).sum(1))
        p_miss = float((np.sum(null >= abs(d_miss.sum())) + 1) / (N_PERM + 1))
    else:
        p_miss = 1.0
    b01 = int(np.sum(sb & ~sv))
    b10 = int(np.sum(~sb & sv))
    p_succ = float(stats.binomtest(b01, b01 + b10, 0.5).pvalue) if b01 + b10 > 0 else 1.0
    return dict(
        n_pairs=len(pairs),
        d_fh_mean=float(d_fh.mean()), d_fh_ci=boot_ci(d_fh, np.mean), p_fh=wil(d_fh),
        d_rms_mean=float(d_rms.mean()), d_rms_ci=boot_ci(d_rms, np.mean), p_rms=wil(d_rms),
        rms_worse_frac=float(np.mean(d_rms > 0)),
        misses_base=int(mb.sum()), misses_var=int(mv.sum()), p_miss=p_miss,
        mcnemar_lost=b01, mcnemar_gained=b10, p_success=p_succ,
    )


def fmt_ci(ci, f='%.3f'):
    return '[' + (f % ci[0]) + ', ' + (f % ci[1]) + ']'


def main():
    tags = sorted(os.path.basename(p)[3:-5] for p in glob.glob(os.path.join(RESULTS, 'mc_S_*.json')))
    tags = [BASE] + [t for t in tags if t != BASE]
    _, base_runs = load(BASE)
    desc, comp = {}, {}
    for t in tags:
        _, runs = load(t)
        desc[t] = describe(runs)
        if t != BASE:
            comp[t] = paired(base_runs, runs)
    for key in ('p_fh', 'p_rms', 'p_miss', 'p_success'):
        adj = holm({t: comp[t][key] for t in comp})
        for t in comp:
            comp[t][key + '_holm'] = adj[t]
    with open(os.path.join(RESULTS, 'stats.json'), 'w') as f:
        json.dump(dict(describe=desc, paired=comp, labels=LABELS), f, indent=1)

    L = ['# Monte Carlo statistics', '',
         'Random admissible engagements (`MC<seed>`), controller at 100 Hz unless stated, plant RK4 at 1 ms.',
         f'Baseline: {desc[BASE]["n"]} seeds; variants: seeds 0-499 paired with the baseline.', '',
         '## Per-configuration performance (95 % confidence intervals)', '',
         '| Configuration | n | Intercepted | Success CI (Wilson) | Pulse hit rate | Cluster-bootstrap CI | '
         'First hit mean [s] | CI | p95 [s] | max [s] | Pointing RMS median [m] | Time in R_eng |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for t in tags:
        d = desc[t]
        L.append(f"| {LABELS.get(t, t)} | {d['n']} | {d['success']}/{d['n']} | {fmt_ci(d['success_ci'], '%.4f')} | "
                 f"{100 * d['pulse_hit_rate']:.2f} % ({d['misses']} misses / {d['pulses']}) | "
                 f"{fmt_ci([100 * c for c in d['pulse_hit_ci_cluster']], '%.2f')} % | {d['fh_mean']:.3f} | "
                 f"{fmt_ci(d['fh_mean_ci'])} | {d['fh_p95']:.2f} | {d['fh_max']:.2f} | {d['rms_median']:.3f} | "
                 f"{100 * d['inR_mean']:.2f} % |")
    L += ['', f"With {desc[BASE]['success']}/{desc[BASE]['n']} interceptions, the one-sided 97.5 % upper bound on "
          f"the baseline failure probability is {100 * desc[BASE]['failure_upper95']:.2f} % (Wilson).", '',
          '## Paired comparison with the final design (variant minus baseline, Holm-adjusted p)', '',
          '| Variant | pairs | Δ first hit [s] | CI | p (Holm) | Δ pointing RMS [m] | CI | p (Holm) | '
          'misses base → variant | p (Holm) | McNemar lost/gained |',
          '|---|---|---|---|---|---|---|---|---|---|---|']
    for t in tags[1:]:
        c = comp[t]
        L.append(f"| {LABELS.get(t, t)} | {c['n_pairs']} | {c['d_fh_mean']:+.3f} | {fmt_ci(c['d_fh_ci'])} | "
                 f"{c['p_fh_holm']:.2g} | {c['d_rms_mean']:+.3f} | {fmt_ci(c['d_rms_ci'])} | {c['p_rms_holm']:.2g} | "
                 f"{c['misses_base']} → {c['misses_var']} | {c['p_miss_holm']:.2g} | "
                 f"{c['mcnemar_lost']}/{c['mcnemar_gained']} |")
    with open(os.path.join(RESULTS, 'stats_report.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')
    sys.stdout.reconfigure(encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
