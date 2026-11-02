"""Additional figures for the detailed experimental report.

python scripts/plot_report_extras.py
Needs results from run_all_experiments.sh (ensemble.npz, stress_*.json, mc_S_final.json).
"""
import json
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scipy.signal as sg  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from guardian import params as P  # noqa: E402
from guardian import turret_controller as tc  # noqa: E402
from guardian.target_motion_generator import _get  # noqa: E402
from plot_scenarios import BAND, C_EST, C_REF, C_TUR, INK, INK2, OUT  # noqa: E402

RES = os.path.join(os.path.dirname(__file__), '..', 'results')
SCEN = ['A', 'B', 'C', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6']
NAMES = {'A': 'A nominal', 'B': 'B maneuver', 'C': 'C loss', 'D1': 'D1 stealth break', 'D2': 'D2 weave',
         'D3': 'D3 close pass', 'D4': 'D4 repeated loss', 'D5': 'D5 combined', 'D6': 'D6 jerk test'}


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'{name}.{ext}'))
    plt.close(fig)
    print('wrote', name)


def fig_scenario_paths():
    fig, axs = plt.subplots(3, 3, figsize=(7.2, 6.0), sharex=True, sharey=True,
                            gridspec_kw=dict(hspace=0.3, wspace=0.08))
    th = np.linspace(np.radians(-10), np.radians(85), 200)
    for a, s in zip(axs.flat, SCEN):
        e = _get(s)
        tr = e['traj']
        k = tr['t'] <= P.T_END
        for rr in (P.R_MIN, P.R_MAX):
            a.plot(rr * np.cos(th), rr * np.sin(th), color=INK2, lw=0.5, ls='--')
        a.plot(tr['x'][k], tr['y'][k], color=C_EST, lw=1.4)
        for t_on, dur in e['blackouts']:
            m = (tr['t'] >= t_on) & (tr['t'] < t_on + dur)
            a.plot(tr['x'][m], tr['y'][m], color=C_REF, lw=3.0, solid_capstyle='butt')
        a.plot(tr['x'][0], tr['y'][0], 'o', color=INK, ms=3)
        a.plot(0, 0, '^', color=INK, ms=4)
        a.set_title(NAMES[s], loc='left', fontsize=7)
        a.set_aspect('equal')
        a.set_xlim(-100, 2100)
        a.set_ylim(-150, 1300)
        a.tick_params(labelsize=6)
    for a in axs[-1]:
        a.set_xlabel('x [m]')
    for a in axs[:, 0]:
        a.set_ylabel('y [m]')
    fig.suptitle('Standardized and stress trajectories (dot = start, orange = radar blackout, '
                 'dashed = 800/2000 m range limits)', x=0.02, ha='left', fontsize=7.5)
    save(fig, 'scenario_paths_grid')


def fig_scenario_accel():
    fig, axs = plt.subplots(3, 3, figsize=(7.2, 5.4), sharex=True, sharey=True,
                            gridspec_kw=dict(hspace=0.35, wspace=0.08))
    for a, s in zip(axs.flat, SCEN):
        e = _get(s)
        tr = e['traj']
        k = tr['t'] <= P.T_END
        for t_on, dur in e['blackouts']:
            a.axvspan(t_on, t_on + dur, color=BAND, lw=0, zorder=0)
        a.plot(tr['t'][k], tr['ax'][k], color=C_EST, lw=1.0, label='$a_x$')
        a.plot(tr['t'][k], tr['ay'][k], color=C_REF, lw=1.0, ls=(0, (4, 2)), label='$a_y$')
        a.plot(tr['t'][k], np.hypot(tr['ax'][k], tr['ay'][k]), color=INK, lw=0.8, label='$|a|$')
        a.axhline(P.A_MAX_TGT, color=INK2, lw=0.5, ls=':')
        a.axhline(-P.A_MAX_TGT, color=INK2, lw=0.5, ls=':')
        jerk = np.hypot(np.diff(tr['ax'][k]), np.diff(tr['ay'][k])).max() / 1e-3
        a.set_title(f'{NAMES[s]} (peak jerk {jerk:.0f} m/s$^3$)', loc='left', fontsize=6.5)
        a.set_ylim(-22, 22)
        a.tick_params(labelsize=6)
    axs[0, 0].legend(frameon=False, fontsize=6, loc='lower left', ncol=3)
    for a in axs[-1]:
        a.set_xlabel('time [s]')
    for a in axs[:, 0]:
        a.set_ylabel('[m/s$^2$]')
    fig.suptitle('Acceleration profiles of each scenario (shaded = blackout, dotted = $\\pm a_{max}$)',
                 x=0.02, ha='left', fontsize=7.5)
    save(fig, 'scenario_accel_grid')


def fig_taylor_roots():
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    markers = ['o', 's', 'D', '^', 'v', 'P']
    cols = [C_EST, C_TUR, C_REF, INK, '#e34948', '#6250d6']
    for n in range(1, 7):
        c = [P.TP ** k / math.factorial(k) for k in range(n, -1, -1)]
        r = np.roots(c)
        ax.plot(r.real, r.imag, markers[n - 1], color=cols[n - 1], ms=5, mfc='white' if n >= 5 else cols[n - 1],
                label=f'$A_{n}$' + (' (unstable)' if n >= 5 else ''))
    ax.axvline(0, color=INK, lw=0.8)
    ax.axvspan(0, 2.5, color='#fbe3dc', zorder=0)
    ax.set_xlim(-5.2, 2.2)
    ax.set_xlabel('Re [rad/s]')
    ax.set_ylabel('Im [rad/s]')
    ax.legend(frameon=False, fontsize=6, loc='upper left', ncol=2)
    ax.set_title('Roots of Taylor sections $A_n(s)$ of $e^{st_p}$, $t_p$ = 0.5 s', loc='left')
    save(fig, 'taylor_roots')


def fig_loop_responses():
    num = [P.N_GEAR * P.K_T]
    den = np.polymul([1.0, 0.0], np.polyadd(np.polymul([P.L_ARM, P.R_ARM], [P.J_EQ, P.B_TURRET]),
                                            [P.N_GEAR ** 2 * P.K_T * P.K_E]))
    kp, ti, tf = tc.CFG['kp'], tc.CFG['ti'], tc.CFG['tau_f']
    C_num = np.polymul([kp * ti, kp], [1.0])
    C_den = np.polymul([ti, 0.0], [tf, 1.0])
    L_num = np.polymul(C_num, num)
    L_den = np.polymul(C_den, den)
    T = sg.lti(L_num, np.polyadd(L_den, L_num))
    S = sg.lti(L_den, np.polyadd(L_den, L_num))
    vel = sg.lti(num, den[:-1])                       # theta_dot / V
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.3), gridspec_kw=dict(wspace=0.4))
    t = np.linspace(0, 0.6, 3000)
    _, y = sg.step(vel, T=t)
    ax[0].plot(t, 1e3 * y, color=C_EST)
    ax[0].set_xlabel('time [s]')
    ax[0].set_ylabel('$\\dot\\theta$ [mrad/s]')
    ax[0].set_title('(a) plant: 1 V step', loc='left')
    t = np.linspace(0, 3, 3000)
    _, y = sg.step(T, T=t)
    ax[1].plot(t, y, color=C_TUR)
    ax[1].axhline(1, color=INK2, lw=0.6, ls='--')
    ax[1].set_xlabel('time [s]')
    ax[1].set_ylabel('$\\theta / y_d$')
    ax[1].set_title('(b) feedback-only step', loc='left')
    w = np.logspace(-1, 3, 2000)
    _, mT = sg.freqresp(T, w)
    _, mS = sg.freqresp(S, w)
    ax[2].semilogx(w, 20 * np.log10(np.abs(mS)), color=C_REF, label='$|S|$')
    ax[2].semilogx(w, 20 * np.log10(np.abs(mT)), color=C_TUR, label='$|T|$')
    ax[2].axhline(0, color=INK2, lw=0.6)
    ax[2].set_xlabel('$\\omega$ [rad/s]')
    ax[2].set_ylabel('[dB]')
    ax[2].set_ylim(-60, 10)
    ax[2].legend(frameon=False, fontsize=6)
    ax[2].set_title('(c) sensitivity', loc='left')
    peak_S = 20 * np.log10(np.abs(mS).max())
    save(fig, 'loop_responses')
    return dict(max_S_dB=float(peak_S), Ms=float(np.abs(mS).max()), max_T_dB=float(20 * np.log10(np.abs(mT).max())))


def fig_ensemble_extras():
    E = dict(np.load(os.path.join(RES, 'ensemble.npz')))
    t = E['t']
    err, est, sig = np.abs(E['err']), np.abs(E['est']), E['sig']
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.4))
    frac = np.mean(err <= P.R_ENG, axis=0)
    ax[0].plot(t, 100 * frac, color=C_TUR)
    ax[0].set_ylim(0, 102)
    ax[0].set_xlabel('time [s]')
    ax[0].set_ylabel('engagements within $R_{eng}$ [%]')
    ax[0].set_title(f'(a) {err.shape[0]} engagements; mean after 2 s: {100 * frac[t >= 2].mean():.1f} %',
                    loc='left')
    m = t >= 0.5
    ax[1].fill_between(t[m], np.nanpercentile(est[:, m], 10, 0), np.nanpercentile(est[:, m], 90, 0),
                       color=C_EST, alpha=0.2, lw=0)
    ax[1].plot(t[m], np.nanmedian(est[:, m], 0), color=C_EST, label='|estimation error| (median, 10-90 %)')
    ax[1].plot(t[m], np.nanmedian(sig[:, m], 0), color=C_REF, ls=(0, (4, 2)), label='filter σ (median)')
    ax[1].set_ylim(0, 1.6)
    ax[1].set_xlabel('time [s]')
    ax[1].set_ylabel('cross-range [m]')
    ax[1].legend(frameon=False, fontsize=6)
    ax[1].set_title('(b) estimator convergence and accuracy', loc='left')
    save(fig, 'mc_ensemble_extras')
    return dict(frac_after2=float(frac[t >= 2].mean()),
                t_med_in=float(t[np.argmax(np.median(err, 0) <= P.R_ENG)]),
                t_90_in=float(t[np.argmax(np.percentile(err, 90, 0) <= P.R_ENG)]),
                est_med_after2=float(np.nanmedian(est[:, t >= 2])),
                sig_med_after2=float(np.nanmedian(sig[:, t >= 2])))


def fig_stress_per_scenario():
    def load(tag):
        with open(os.path.join(RES, f'stress_{tag}.json')) as f:
            return json.load(f)
    variants = [('final', 'order 4 (final), 100 Hz', C_EST), ('order3', 'order 3, 100 Hz', C_REF),
                ('final_20Hz', 'order 4, 20 Hz', C_TUR), ('order3_20Hz', 'order 3, 20 Hz', INK2)]
    fig, ax = plt.subplots(2, 1, figsize=(7.2, 4.0), sharex=True, gridspec_kw=dict(hspace=0.15))
    x = np.arange(len(SCEN))
    wbar = 0.2
    out = {}
    for j, (tag, lab, col) in enumerate(variants):
        rows = load(tag)
        for i, (key, a) in enumerate((('point_rms_m', ax[0]), ('t_first_hit', ax[1]))):
            means, lo, hi = [], [], []
            for s in SCEN:
                v = np.array([r[key] for r in rows if r['scenario'] == s])
                mu = v.mean()
                se = v.std(ddof=1) / math.sqrt(len(v))
                means.append(mu)
                lo.append(1.96 * se)
                hi.append(1.96 * se)
            a.bar(x + (j - 1.5) * wbar, means, wbar, yerr=[lo, hi], color=col, label=lab,
                  error_kw=dict(lw=0.6, capsize=1.5), edgecolor='white', linewidth=0.4)
            out.setdefault(tag, {})[key] = dict(zip(SCEN, map(float, means)))
        out[tag]['misses'] = {s: int(sum(r['n_pulses'] - r['n_hits'] for r in rows if r['scenario'] == s))
                              for s in SCEN}
    ax[0].set_ylabel('pointing RMS [m]')
    ax[0].set_ylim(0, 1.3)
    ax[1].set_ylabel('first hit [s]')
    ax[1].set_xticks(x)
    ax[1].set_xticklabels([NAMES[s] for s in SCEN], fontsize=6.5, rotation=20)
    ax[0].legend(frameon=False, fontsize=6, ncol=4, loc='upper left')
    ax[0].set_title('Stress ensembles: 50 noise realisations per scenario (mean ± 95 % CI)', loc='left')
    save(fig, 'stress_per_scenario')
    return out


def fig_pulses():
    with open(os.path.join(RES, 'mc_S_final.json')) as f:
        runs = json.load(f)['runs']
    npul = [r['n_pulses'] for r in runs]
    ph = np.array([p[2] for r in runs for p in r['pulses']])
    modes = [p[3] for r in runs for p in r['pulses']]
    gaps = []
    for r in runs:
        ts = [p[0] for p in r['pulses']]
        gaps += list(np.diff(ts))
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.2), gridspec_kw=dict(wspace=0.4))
    ax[0].hist(npul, bins=np.arange(min(npul) - 0.5, max(npul) + 1.5), color=C_EST, edgecolor='white')
    ax[0].set_xlabel('pulses per engagement')
    ax[0].set_ylabel('engagements')
    ax[0].set_title('(a) firing cadence', loc='left')
    ax[1].hist(np.clip(1 - ph, 1e-7, 1), bins=np.logspace(-7, 0, 36), color=C_REF, edgecolor='white')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('$1-P_{hit}$ at the fire instant')
    ax[1].set_title('(b) confidence when firing', loc='left')
    ax[2].hist(gaps, bins=np.linspace(1.0, 2.0, 41), color=C_TUR, edgecolor='white')
    ax[2].set_xlabel('time between pulses [s]')
    ax[2].set_title('(c) wait beyond 1 s recharge', loc='left')
    save(fig, 'mc_pulses')
    g = np.array(gaps)
    return dict(pulses_mean=float(np.mean(npul)), pulses_min=int(min(npul)), pulses_max=int(max(npul)),
                coast_shots=modes.count('coast'), track_shots=modes.count('track'),
                gap_median=float(np.median(g)), gap_frac_le_1p01=float(np.mean(g <= 1.01)),
                gap_p95=float(np.percentile(g, 95)), frac_phit_ge_0999=float(np.mean(ph >= 0.999)))


def fig_interface():
    path = os.path.join(RES, 'interface_robustness.json')
    if not os.path.exists(path):
        return None
    with open(path) as f:
        S = json.load(f)['summary']
    ids = list(S)
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=True, gridspec_kw=dict(wspace=0.06))
    y = np.arange(len(ids))[::-1]
    for i, v in enumerate(ids):
        s = S[v]
        ax[0].errorbar(s['fh_mean'], y[i], xerr=[[s['fh_mean'] - s['fh_ci'][0]], [s['fh_ci'][1] - s['fh_mean']]],
                       fmt='o', color=C_EST, ms=4, capsize=2)
        ax[1].barh(y[i], 100 * s['hit_rate'], color=C_TUR, height=0.6)
        ax[1].text(100 * s['hit_rate'] + 0.02, y[i], f"{s['success']}/{s['n']} int., {s['misses']} miss",
                   va='center', fontsize=6)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels([f"{v}: {S[v]['label']}" for v in ids], fontsize=6)
    ax[0].set_xlabel('first hit mean [s] (95 % CI)')
    ax[1].set_xlim(99.0, 100.6)
    ax[1].set_xlabel('pulse hit rate [%]')
    fig.suptitle('Robustness to evaluation-script conventions (300 engagements each)', x=0.02, ha='left',
                 fontsize=7.5)
    save(fig, 'interface_robustness')
    return True


def main():
    fig_scenario_paths()
    fig_scenario_accel()
    fig_taylor_roots()
    out = dict(loop=fig_loop_responses(), ensemble=fig_ensemble_extras(),
               stress=fig_stress_per_scenario(), pulses=fig_pulses())
    fig_interface()
    with open(os.path.join(RES, 'report_extras.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != 'stress'}, indent=1))


if __name__ == '__main__':
    main()
