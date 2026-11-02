"""Write paper/numbers.tex: every number quoted in the paper, taken from the result files.

python scripts/export_paper_numbers.py
python scripts/export_paper_numbers.py --results results/superseded --draft   # preliminary build
Also copies the figures used by the paper into paper/figures/.
"""
import argparse
import json
import math
import os
import re
import shutil

import numpy as np
from scipy import stats

HERE = os.path.dirname(__file__)
RES = os.path.join(HERE, '..', 'results')
PAPER = os.path.join(HERE, '..', '..', 'paper')
FIGS = ['scenario_C', 'mc_first_hit', 'mc_blackout_aligned', 'mc_calibration', 'mc_ablation_forest',
        'mc_error_bands']


def load(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)


def p_tex(p):
    """p-value in TeX: plain to 2 s.f. above 1e-3, else a power of ten."""
    if p >= 1e-3:
        return f'{p:.2g}'
    e = math.floor(math.log10(p))
    m = p / 10 ** e
    return f'{m:.1f}\\times10^{{{e}}}'


def main():
    global RES
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', default=RES)
    ap.add_argument('--draft', action='store_true', help='mark the PDF as preliminary')
    a = ap.parse_args()
    RES = a.results
    S = load('stats.json')
    d = S['describe']['S_final']
    P = S['paired']
    ST = load('stress_stats.json')
    E = load('ensemble_stats.json')
    if 'psure0999 vs final' not in ST['compare']:
        # older run: final used 0.999 and the variant 0.99, so flip the comparison
        c = ST['compare']['psure099 vs final']
        ST['compare']['psure0999 vs final'] = dict(
            d_fh=-c['d_fh'], p_fh_holm=c['p_fh_holm'],
            misses_base=c['misses_var'], misses_var=c['misses_base'])
    runs = load('mc_S_final.json')['runs']

    feats = [r['features'] for r in runs]
    fh = [r['t_first_hit'] for r in runs]
    rms = [r['point_rms_m'] for r in runs]
    miss = np.array([p[1] for r in runs for p in r['pulses']])
    ens = np.load(os.path.join(RES, 'ensemble.npz'))
    t_on = ens['blackouts'][:, 0, 0]
    n_bo = int(np.sum(np.isfinite(t_on) & (t_on >= 1.5) & (t_on + 1.5 <= 10.0)))

    c3, c3_20 = ST['compare']['order3 vs final'], ST['compare']['order3_20Hz vs final_20Hz']
    cps = ST['compare']['psure0999 vs final']
    sf = ST['summary']['final']
    r20 = S['describe']['S_rate20']

    N = {
        'mcN': d['n'],
        'mcFailUB': f"{100 * d['failure_upper95']:.2f}",
        'mcHit': f"{100 * d['pulse_hit_rate']:.2f}",
        'mcHitLo': f"{100 * d['pulse_hit_ci_cluster'][0]:.2f}",
        'mcHitHi': f"{100 * d['pulse_hit_ci_cluster'][1]:.2f}",
        'mcMiss': d['misses'], 'mcPulses': d['pulses'],
        'mcFh': f"{d['fh_mean']:.2f}", 'mcFhLo': f"{d['fh_mean_ci'][0]:.2f}", 'mcFhHi': f"{d['fh_mean_ci'][1]:.2f}",
        'mcFhNinetyFive': f"{d['fh_p95']:.2f}", 'mcFhMax': f"{d['fh_max']:.2f}",
        'mcRms': f"{d['rms_median']:.2f}", 'mcInR': f"{100 * d['inR_mean']:.1f}", 'mcVpk': f"{d['V_peak']:.0f}",
        'mcRhoSlew': f"{stats.spearmanr([f['slew_deg'] for f in feats], fh).statistic:.2f}",
        'mcRhoArms': f"{stats.spearmanr([f['a_rms'] for f in feats], rms).statistic:.2f}",
        'mcMissMed': f'{np.median(miss):.2f}', 'mcMissNinetyNine': f'{np.percentile(miss, 99):.2f}',
        'nBlackouts': n_bo,
        'ezTrack': f"{E['consistency']['track_Ez2']:.2f}", 'ezCoast': f"{E['consistency']['coast_Ez2']:.2f}",
        'brier': f"{E['calibration']['brier']:.4f}",
        'tzFh': f"{r20['fh_mean']:.2f}", 'tzHit': f"{100 * r20['pulse_hit_rate']:.2f}",
        'abOtwoRms': f"{P['S_order2']['d_rms_mean']:.2f}",
        'abOtwoMiss': f"{P['S_order2']['misses_base']}$\\to${P['S_order2']['misses_var']}",
        'abOtwoLost': P['S_order2']['mcnemar_lost'],
        'abPtenFh': f"{P['S_pdev10']['d_fh_mean']:.2f}",
        'abPthirtyRms': f"{abs(P['S_pdev30']['d_rms_mean']):.3f}", 'abPthirtyFh': f"{P['S_pdev30']['d_fh_mean']:.2f}",
        'abAmbRms': f"{P['S_ambush']['d_rms_mean']:.2f}",
        'abOthreeRms': f"{P['S_order3']['d_rms_mean']:+.3f}", 'abOthreePrms': p_tex(P['S_order3']['p_rms_holm']),
        'stHit': f"{100 * sf['hit_rate']:.2f}", 'stMiss': sf['misses'], 'stFh': f"{sf['fh_mean']:.2f}",
        'stOthreeRms': f"{c3['d_rms']:+.3f}", 'stOthreePrms': p_tex(c3['p_rms_holm']),
        'stDsixOthree': f"{c3['per_scenario']['D6']['d_rms']:.2f}",
        'stTwentyOthreeRms': f"{c3_20['d_rms']:+.3f}", 'stTwentyOthreeP': p_tex(c3_20['p_rms_holm']),
        'stTwentyMiss': f"{c3_20['misses_base']}$\\to${c3_20['misses_var']}",
        'stPsureFh': f"{cps['d_fh']:.3f}", 'stPsureP': p_tex(cps['p_fh_holm']),
        'stPsureMiss': f"{cps['misses_base']}$\\to${cps['misses_var']}",
    }

    # generalisation evidence: pooled nominal-setting engagements, slew model, interface study
    D = S['describe']
    pooled = [(D[t]['n'], D[t]['success']) for t in ('S_final', 'S_rate20', 'S_theta60')]
    pooled += [(ST['summary'][t]['n'], ST['summary'][t]['success']) for t in ('final', 'final_20Hz')]
    n_pool = sum(n for n, _ in pooled)
    k_pool = sum(k for _, k in pooled)
    z = 1.959964
    fail_ub = (z * z / (2 * n_pool) + z * math.sqrt(z * z / (4 * n_pool ** 2))) / (1 + z * z / n_pool) \
        if k_pool == n_pool else float('nan')
    slew = np.array([f['slew_deg'] for f in feats])
    fh_arr = np.array(fh, dtype=float)
    A = np.vstack([np.ones_like(slew), slew]).T
    coef = np.linalg.lstsq(A, fh_arr, rcond=None)[0]
    res = fh_arr - A @ coef
    N.update(pooledN=f'{n_pool:,}'.replace(',', '\\,'), pooledOK=f'{k_pool:,}'.replace(',', '\\,'),
             pooledFail=f'{100 * fail_ub:.2f}',
             slewA=f'{coef[0]:.2f}', slewB=f'{coef[1]:.4f}', slewR=f'{1 - res.var() / fh_arr.var():.2f}',
             slewPI=f'{1.96 * res.std(ddof=2):.2f}', slewWorst=f'{coef[0] + coef[1] * 95:.1f}')
    ipath = os.path.join(HERE, '..', 'results', 'interface_robustness.json')
    if os.path.exists(ipath):
        with open(ipath) as f:
            IS = json.load(f)['summary']
        N.update(ifN=f"{sum(s['n'] for s in IS.values()):,}".replace(',', '\\,'),
                 ifOK=f"{sum(s['success'] for s in IS.values()):,}".replace(',', '\\,'),
                 ifHitMin=f"{100 * min(s['hit_rate'] for s in IS.values()):.2f}",
                 ifVariants=len(IS))

    rows = []
    with open(os.path.join(HERE, '..', 'results', 'scenarios.txt')) as f:
        for line in f:
            m = re.match(r'(\w+): first hit ([\d.]+) s, (\d+)/(\d+) hits, pointing RMS ([\d.]+) m', line)
            if m:
                rows.append(f'{m[1]} & {float(m[2]):.2f} & {m[3]}/{m[4]} & {m[5]}\\\\')
    os.makedirs(PAPER, exist_ok=True)
    with open(os.path.join(PAPER, 'numbers.tex'), 'w', encoding='utf-8') as f:
        f.write('% Generated by python/scripts/export_paper_numbers.py -- do not edit by hand.\n')
        for k, v in N.items():
            f.write(f'\\newcommand{{\\{k}}}{{{v}}}\n')
        f.write('\\newcommand{\\scenarioRows}{%\n' + '\n'.join(rows) + '}\n')
        note = ('\\\\\\textcolor{red}{PRELIMINARY DRAFT: numbers from the previous run '
                '(instant-fire threshold 0.999); final numbers pending}') if a.draft else ''
        f.write('\\newcommand{\\draftnote}{' + note + '}\n')
    os.makedirs(os.path.join(PAPER, 'figures'), exist_ok=True)
    fig_dir = os.path.join(HERE, '..', 'results', 'figures')
    for name in FIGS:
        shutil.copy(os.path.join(fig_dir, f'{name}.pdf'), os.path.join(PAPER, 'figures', f'{name}.pdf'))
    print(f'wrote {len(N)} numbers, {len(rows)} scenario rows, {len(FIGS)} figures')


if __name__ == '__main__':
    main()
