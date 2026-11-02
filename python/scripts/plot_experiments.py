"""Experimental figures for the report and paper (Monte Carlo campaign).

Requires: results/mc_S_*.json (run_ablations.sh), results/ensemble.npz (mc_ensemble.py),
results/stats.json (statistics.py), results/rate_study.txt (rate_study.py).
python scripts/plot_experiments.py
"""
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from guardian import params as P  # noqa: E402
from guardian.target_motion_generator import _get  # noqa: E402
from plot_scenarios import BAND, C_EST, C_REF, C_TUR, INK, INK2, OUT  # noqa: E402

RES = os.path.join(os.path.dirname(__file__), '..', 'results')
C4 = '#eda100'         # slot 4 yellow (only where a 4th series is unavoidable)


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'{name}.{ext}'))
    plt.close(fig)
    print('wrote', name)


def load_runs(tag):
    with open(os.path.join(RES, f'mc_{tag}.json')) as f:
        return sorted(json.load(f)['runs'], key=lambda r: r['seed'])


# ---------------------------------------------------------------------------- path generation
def fig_paths(n_show=80):
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw=dict(width_ratios=[1.25, 1]))
    a = ax[0]
    th = np.linspace(np.radians(-10), np.radians(85), 200)
    for rr in (P.R_MIN, P.R_MAX):
        a.plot(rr * np.cos(th), rr * np.sin(th), color=INK2, lw=0.7, ls='--')
    for s in range(n_show):
        tr = _get(f'MC{s}')['traj']
        k = tr['t'] <= P.T_END
        a.plot(tr['x'][k], tr['y'][k], color=C_EST, lw=0.5, alpha=0.45)
        a.plot(tr['x'][0], tr['y'][0], 'o', ms=2, color=INK)
    a.plot(0, 0, '^', color=INK, ms=6)
    a.set_aspect('equal')
    a.set_xlabel('x [m]')
    a.set_ylabel('y [m]')
    a.set_title(f'(a) {n_show} random admissible trajectories (dots = start)', loc='left')
    a.text(150, 150, 'turret', fontsize=6.5, color=INK2)

    a = ax[1]
    for s in range(12):
        e = _get(f'MC{s}')
        tr = e['traj']
        k = tr['t'] <= P.T_END
        a.plot(tr['t'][k], np.hypot(tr['ax'][k], tr['ay'][k]), lw=0.8, color=C_REF, alpha=0.6)
    a.axhline(P.A_MAX_TGT, color=INK, lw=0.7, ls='--')
    a.text(0.2, P.A_MAX_TGT + 0.6, '$a_{max}$ = 20 m/s$^2$', fontsize=6.5, color=INK)
    a.set_xlabel('time [s]')
    a.set_ylabel('$|a|$ [m/s$^2$]')
    a.set_ylim(0, 23)
    a.set_title('(b) acceleration profiles, 12 runs', loc='left')
    save(fig, 'mc_paths')


def fig_features(runs):
    F = [r['features'] for r in runs]
    bo_on = [b[0] for f in F for b in f['blackouts']]
    bo_len = [b[1] for f in F for b in f['blackouts']]
    panels = [
        ([f['r0'] for f in F], 'initial range [m]'),
        ([f['th0_deg'] for f in F], 'initial elevation [deg]'),
        ([f['speed0'] for f in F], 'initial speed [m/s]'),
        ([f['a_max'] for f in F], 'peak |a| [m/s$^2$]'),
        ([f['jerk_max'] for f in F], 'peak |jerk| [m/s$^3$]'),
        ([f['los_rate_max'] for f in F], 'peak LOS rate [rad/s]'),
        (bo_on, 'blackout onset [s] (per interval)'),
        (bo_len, 'blackout duration [s] (per interval)'),
    ]
    fig, ax = plt.subplots(2, 4, figsize=(7.2, 3.3), gridspec_kw=dict(hspace=0.55, wspace=0.35))
    for a, (v, lab) in zip(ax.flat, panels):
        a.hist(v, bins=25, color=C_EST, edgecolor='white', linewidth=0.5)
        a.set_xlabel(lab)
        a.tick_params(labelsize=6)
    ax[0, 0].set_ylabel('count')
    ax[1, 0].set_ylabel('count')
    n2 = sum(f['n_blackouts'] == 2 for f in F)
    fig.suptitle(f'Monte Carlo design space: {len(F)} engagements, {len(bo_on)} blackout intervals '
                 f'({n2} runs with two)', x=0.02, ha='left', fontsize=8)
    save(fig, 'mc_features')


# ---------------------------------------------------------------------------- outcomes
def ecdf(v):
    v = np.sort(np.asarray(v))
    return v, np.arange(1, len(v) + 1) / len(v)


def fig_first_hit(base, r20, t60):
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.5))
    for runs, col, ls, lab in ((base, C_EST, '-', 'final design, 100 Hz'),
                               (r20, C_REF, (0, (4, 2)), '20 Hz controller'),
                               (t60, C_TUR, (0, (1, 1.2)), 'turret starts at 60°')):
        v = [r['t_first_hit'] for r in runs if np.isfinite(r['t_first_hit'])]
        x, y = ecdf(v)
        ax[0].step(x, y, where='post', color=col, ls=ls, lw=1.3, label=f'{lab} (n={len(runs)})')
    ax[0].set_xlabel('time of first successful pulse [s]')
    ax[0].set_ylabel('fraction of engagements')
    ax[0].set_xlim(0, 4)
    ax[0].legend(frameon=False, loc='center right', fontsize=6.5)
    ax[0].set_title('(a) interception-time ECDF', loc='left')

    slew = [r['features']['slew_deg'] for r in base]
    fh = [r['t_first_hit'] for r in base]
    ax[1].scatter(slew, fh, s=6, color=C_EST, alpha=0.5, lw=0)
    sl = np.linspace(0, 60, 50)
    ax[1].plot(sl, np.radians(sl) / 0.58, color=INK2, lw=0.8, ls='--')
    ax[1].text(44, 0.95, 'pure slew\nat 0.58 rad/s', fontsize=6.5, color=INK2, ha='left')
    rho = stats.spearmanr(slew, fh)
    ax[1].set_xlabel('initial slew |θ$_T$(0) − θ(0)| [deg]')
    ax[1].set_ylabel('first hit [s]')
    ax[1].set_title(f'(b) acquisition-limited (Spearman ρ = {rho.statistic:.2f})', loc='left')
    save(fig, 'mc_first_hit')


def fig_miss_hist(base):
    miss = np.array([p[1] for r in base for p in r['pulses']])
    fig, ax = plt.subplots(figsize=(3.5, 2.2))
    ax.hist(miss, bins=np.linspace(0, 3, 61), color=C_EST, edgecolor='white', linewidth=0.4)
    ax.axvline(P.R_ENG, color=INK, lw=0.8, ls='--')
    ax.set_yscale('log')
    ax.set_xlabel('beam-to-target miss distance at each pulse [m]')
    ax.set_ylabel('pulses (log)')
    n_out = int(np.sum(miss > P.R_ENG))
    ax.set_title(f'{len(miss)} pulses: median {np.median(miss):.2f} m, '
                 f'p99 {np.percentile(miss, 99):.2f} m, {n_out} > $R_{{eng}}$', loc='left')
    save(fig, 'mc_miss_hist')


def fig_difficulty(base):
    jerk = np.array([r['features']['jerk_max'] for r in base])
    arms = np.array([r['features']['a_rms'] for r in base])
    rms = np.array([r['point_rms_m'] for r in base])
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.4), sharey=True)
    for a, v, lab in ((ax[0], arms, 'RMS target acceleration [m/s$^2$]'),
                      (ax[1], jerk, 'peak target jerk [m/s$^3$]')):
        a.scatter(v, rms, s=6, color=C_EST, alpha=0.5, lw=0)
        rho = stats.spearmanr(v, rms)
        a.set_xlabel(lab)
        a.set_title(f'Spearman ρ = {rho.statistic:.2f}, p = {rho.pvalue:.1g}', loc='left')
    ax[0].set_ylabel('pointing RMS after 2 s [m]')
    save(fig, 'mc_difficulty')


# ---------------------------------------------------------------------------- ensemble time series
def fig_error_bands(E):
    t, err = E['t'], np.abs(E['err'])
    fig, ax = plt.subplots(figsize=(3.5, 2.4))
    q50, q90, q99 = (np.nanpercentile(err, q, axis=0) for q in (50, 90, 99))
    ax.fill_between(t, 0, q99, color=C_EST, alpha=0.15, lw=0, label='99th percentile')
    ax.fill_between(t, 0, q90, color=C_EST, alpha=0.30, lw=0, label='90th percentile')
    ax.plot(t, q50, color=C_EST, lw=1.2, label='median')
    ax.axhline(P.R_ENG, color=INK, lw=0.8, ls='--')
    ax.text(9.9, P.R_ENG + 0.1, '$R_{eng}$', ha='right', fontsize=6.5)
    ax.set_ylim(0, 6)
    ax.set_xlim(0, 10)
    ax.set_xlabel('time [s]')
    ax.set_ylabel('|cross-range pointing error| [m]')
    ax.set_title(f'Ensemble of {err.shape[0]} engagements', loc='left')
    ax.legend(frameon=False, loc='upper right')
    save(fig, 'mc_error_bands')


def fig_blackout_aligned(E):
    t, f = E['t'], float(E['f_ctrl'])
    rel = np.arange(-0.5, 1.5 + 1e-9, 1 / f)
    segs = {k: [] for k in ('err', 'est', 'sig')}
    for i in range(E['err'].shape[0]):
        t_on = E['blackouts'][i, 0, 0]
        if not np.isfinite(t_on) or t_on + 1.5 > 10.0 or t_on < 1.5:
            continue
        k0 = int(round((t_on - 0.5) * f))
        sl = slice(k0, k0 + len(rel))
        for k in segs:
            segs[k].append(np.abs(E[k][i, sl]) if k != 'sig' else E[k][i, sl])
    durs = E['blackouts'][:, 0, 1]
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    for k, col, lab in (('err', C_TUR, '|pointing error|'), ('est', C_EST, '|estimation error|'),
                        ('sig', C_REF, 'filter σ (cross-range)')):
        A = np.array(segs[k])
        ax.fill_between(rel, np.percentile(A, 10, 0), np.percentile(A, 90, 0), color=col, alpha=0.18, lw=0)
        ax.plot(rel, np.median(A, 0), color=col, lw=1.2, label=lab)
    ax.axvspan(0, np.nanmedian(durs), color=BAND, zorder=0)
    ax.axhline(P.R_ENG, color=INK, lw=0.7, ls='--')
    ax.set_xlabel('time since blackout onset [s]')
    ax.set_ylabel('[m]')
    ax.set_ylim(0, 4)
    ax.legend(frameon=False, loc='upper right')
    ax.set_title(f'{len(segs["err"])} blackouts aligned (median, 10–90 % band)', loc='left')
    save(fig, 'mc_blackout_aligned')


def fig_consistency(E):
    mode, est, sig = E['mode'], E['est'], E['sig']
    t = E['t']
    track = (mode == 1) & (t[None, :] >= 2.0) & np.isfinite(sig) & (sig > 0)
    coast = (mode == 2) & np.isfinite(sig) & (sig > 0)
    z_tr = (est[track] / sig[track]).astype(float)
    z_co = (est[coast] / sig[coast]).astype(float)
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.4))
    x = np.linspace(-5, 5, 400)
    for a, z, title in ((ax[0], z_tr, 'tracking'), (ax[1], z_co, 'coasting (measurement loss)')):
        a.hist(np.clip(z, -5, 5), bins=np.linspace(-5, 5, 81), density=True, color=C_EST,
               edgecolor='white', linewidth=0.3, label='normalised error')
        a.plot(x, stats.norm.pdf(x), color=INK, lw=1.0, ls='--', label='N(0,1)')
        a.set_xlabel('estimation error / filter σ')
        a.set_title(f'({"ab"[title != "tracking"]}) {title}: E[z²] = {np.mean(z ** 2):.2f}, '
                    f'|z|<2: {100 * np.mean(np.abs(z) < 2):.1f} %', loc='left')
    ax[0].set_ylabel('density')
    ax[0].legend(frameon=False, loc='upper left')
    save(fig, 'mc_consistency')
    return dict(track_Ez2=float(np.mean(z_tr ** 2)), coast_Ez2=float(np.mean(z_co ** 2)),
                track_in2=float(np.mean(np.abs(z_tr) < 2)), coast_in2=float(np.mean(np.abs(z_co) < 2)))


def fig_calibration(E):
    t = E['t']
    m = (E['mode'] >= 1) & (t[None, :] >= 1.0)
    p = E['phit'][m].astype(float)
    y = (np.abs(E['err'][m]) <= P.R_ENG).astype(float)
    edges = np.array([0, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999, 1.0001])
    idx = np.digitize(p, edges) - 1
    xs, ys, ns, lo, hi = [], [], [], [], []
    for b in range(len(edges) - 1):
        sel = idx == b
        n = int(sel.sum())
        if n < 30:
            continue
        k = int(y[sel].sum())
        xs.append(p[sel].mean())
        ys.append(k / n)
        ns.append(n)
        z = 1.96
        c = (k / n + z * z / (2 * n)) / (1 + z * z / n)
        h = z * np.sqrt((k / n) * (1 - k / n) / n + z * z / (4 * n * n)) / (1 + z * z / n)
        lo.append(c - h)
        hi.append(c + h)
    brier = float(np.mean((p - y) ** 2))
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    ax.plot([0, 1], [0, 1], color=INK2, lw=0.8, ls='--')
    yerr = [np.maximum(np.array(ys) - np.array(lo), 0), np.maximum(np.array(hi) - np.array(ys), 0)]
    ax.errorbar(xs, ys, yerr=yerr, fmt='o', ms=4,
                color=C_EST, mfc=C_EST, mec='white', capsize=2, lw=1)
    for x_, y_, n_ in zip(xs, ys, ns):
        if x_ < 0.95:
            ax.annotate(f'{n_}', (x_, y_), xytext=(4, -9), textcoords='offset points', fontsize=5.5, color=INK2)
    hi_bins = [(x_, n_) for x_, n_ in zip(xs, ns) if x_ >= 0.95]
    note = 'samples in top bins:\n' + '\n'.join(f'P≈{x_:.3f}: {n_}' for x_, n_ in hi_bins)
    ax.text(0.98, 0.04, note, transform=ax.transAxes, ha='right', va='bottom', fontsize=5.5, color=INK2)
    ax.set_xlabel('predicted hit probability $P_{hit}$')
    ax.set_ylabel('observed fraction within $R_{eng}$')
    ax.set_title(f'Calibration of $P_{{hit}}$ (Brier {brier:.4f}; n per bin)', loc='left')
    save(fig, 'mc_calibration')
    return dict(brier=brier, bins=list(zip(map(float, xs), map(float, ys), ns)))


# ---------------------------------------------------------------------------- ablations
def fig_forest(S):
    comp, labels = S['paired'], S['labels']
    order = sorted(comp, key=lambda t: comp[t]['d_rms_mean'])
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 3.0), sharey=True, gridspec_kw=dict(wspace=0.08))
    y = np.arange(len(order))
    for a, key, lab in ((ax[0], 'd_fh', 'Δ first-hit time [s]'), (ax[1], 'd_rms', 'Δ pointing RMS [m]')):
        for i, t in enumerate(order):
            c = comp[t]
            m, ci = c[key + '_mean'], c[key + '_ci']
            sig = c['p_' + key[2:] + '_holm'] < 0.05
            a.plot(ci, [i, i], color=C_EST if sig else INK2, lw=1.4)
            a.plot(m, i, 'o', ms=4.5, color=C_EST if sig else 'white', mec=C_EST if sig else INK2, mew=1)
        a.axvline(0, color=INK, lw=0.7)
        a.set_xlabel(lab)
    ax[1].set_xscale('symlog', linthresh=0.02, linscale=0.6)
    ticks = [-0.1, -0.02, 0, 0.02, 0.1, 0.5, 1.0]
    ax[1].set_xticks(ticks)
    ax[1].set_xticklabels(['−0.1', '−0.02', '0', '0.02', '0.1', '0.5', '1'], fontsize=6)
    ax[1].set_xlim(-0.15, 1.3)
    misses_b = [comp[t]['misses_base'] for t in order]
    misses_v = [comp[t]['misses_var'] for t in order]
    ax[2].barh(y + 0.18, misses_v, height=0.34, color=C_REF, label='variant')
    ax[2].barh(y - 0.18, misses_b, height=0.34, color=C_EST, label='final design')
    for i, t in enumerate(order):
        star = '*' if comp[t]['p_miss_holm'] < 0.05 else ''
        ax[2].text(misses_v[i] + 1.5, y[i] + 0.18, f'{misses_v[i]}{star}', va='center', fontsize=6, color=INK)
        ax[2].text(misses_b[i] + 1.5, y[i] - 0.18, f'{misses_b[i]}', va='center', fontsize=6, color=INK2)
    ax[2].set_xlim(0, max(misses_v + misses_b) * 1.25)
    ax[2].set_xlabel('missed pulses (* p < 0.05)')
    ax[2].legend(frameon=False, loc='lower right', fontsize=6)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels([labels.get(t, t) for t in order], fontsize=6.5)
    fig.suptitle('Paired one-factor ablations vs final design (500 seeds; filled = Holm-adjusted p < 0.05; '
                 'bars: 95 % bootstrap CI)', x=0.02, ha='left', fontsize=7.5)
    save(fig, 'mc_ablation_forest')


def fig_rate():
    path = os.path.join(RES, 'rate_study.txt')
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if parts and parts[0].isdigit():
                rows.append((int(parts[0]), float(parts[3][:-1]), float(parts[4][:-1]), float(parts[5][:-1]),
                             parts[2]))
    f_, fh_mean, fh_max, rms, hits = zip(*rows)
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.3))
    ax[0].plot(f_, fh_mean, 'o-', color=C_EST, ms=4, label='mean')
    ax[0].plot(f_, fh_max, 's--', color=C_REF, ms=4, label='max')
    ax[0].set_xscale('log')
    ax[0].set_xlabel('controller call rate [Hz]')
    ax[0].set_ylabel('first hit [s]')
    ax[0].legend(frameon=False)
    ax[0].set_title('(a) interception time (24 cases per rate)', loc='left')
    ax[1].plot(f_, rms, 'o-', color=C_TUR, ms=4)
    ax[1].set_xscale('log')
    ax[1].set_xlabel('controller call rate [Hz]')
    ax[1].set_ylabel('median pointing RMS [m]')
    ax[1].set_title('(b) pointing accuracy; all rates 24/24 intercepted', loc='left')
    for a in ax:
        a.axvspan(100, 1000, color=BAND, zorder=0)
    ax[1].text(300, max(rms) * 0.97, 'recommended', ha='center', fontsize=6.5, color=INK2)
    save(fig, 'rate_study')


def main():
    base = load_runs('S_final')
    fig_paths()
    fig_features(base)
    fig_first_hit(base, load_runs('S_rate20'), load_runs('S_theta60'))
    fig_miss_hist(base)
    fig_difficulty(base)
    E = dict(np.load(os.path.join(RES, 'ensemble.npz')))
    fig_error_bands(E)
    fig_blackout_aligned(E)
    cons = fig_consistency(E)
    cal = fig_calibration(E)
    with open(os.path.join(RES, 'stats.json')) as f:
        S = json.load(f)
    fig_forest(S)
    fig_rate()
    with open(os.path.join(RES, 'ensemble_stats.json'), 'w') as f:
        json.dump(dict(consistency=cons, calibration=cal), f, indent=1)
    print(json.dumps(dict(consistency=cons, calibration=cal), indent=1))


if __name__ == '__main__':
    main()
