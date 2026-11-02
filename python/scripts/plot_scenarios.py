"""Time-domain figures for the standardized (A, B, C) and stress (D1-D6) scenarios.

python scripts/plot_scenarios.py               # all scenarios at 1 kHz
python scripts/plot_scenarios.py A B C --f-ctrl 100
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from guardian import params as P  # noqa: E402
from guardian.sim import run  # noqa: E402
from guardian.target_motion_generator import blackout_intervals  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'figures')

INK = '#0b0b0b'
INK2 = '#52514e'
GRID = '#e4e3df'
BAND = '#efeee9'
C_EST = '#2a78d6'      # slot 1 blue
C_REF = '#eb6834'      # slot 2 orange
C_TUR = '#1baf7a'      # slot 3 aqua

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8,
    'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.2,
    'savefig.dpi': 300, 'savefig.bbox': 'tight',
})


def _shade_blackouts(ax, scenario):
    for (t_on, dur) in blackout_intervals(scenario):
        ax.axvspan(t_on, t_on + dur, color=BAND, zorder=0, lw=0)


def plot_scenario(scenario, f_ctrl):
    o = run(scenario, f_ctrl=f_ctrl)
    t = o['t']
    m = o['metrics']
    deg = np.degrees
    fig, ax = plt.subplots(4, 1, figsize=(3.5, 5.6), sharex=True,
                           gridspec_kw=dict(height_ratios=[1.4, 1.2, 0.8, 0.8], hspace=0.12))

    a = ax[0]
    _shade_blackouts(a, scenario)
    a.plot(t, deg(o['th_true']), color=INK, lw=1.0, label='true $\\theta_T$')
    a.plot(t, deg(o['th_hat']), color=C_EST, lw=1.0, ls=(0, (4, 2)), label='estimate $\\hat\\theta$')
    a.plot(t, deg(o['r_ref']), color=C_REF, lw=1.0, ls=(0, (1, 1.5)), label='reference $\\hat\\theta(t+t_p)$')
    a.plot(t, deg(o['th']), color=C_TUR, lw=1.4, alpha=0.9, label='turret $\\theta$')
    a.set_ylabel('elevation [deg]')
    a.legend(loc='lower right', frameon=False, ncol=2, handlelength=2.2)
    a.set_title(f'Scenario {scenario}: first hit {m["t_first_hit"]:.2f} s, '
                f'{m["n_hits"]}/{m["n_pulses"]} pulses hit', loc='left', color=INK)

    a = ax[1]
    _shade_blackouts(a, scenario)
    err = o['r_true'] * np.sin(o['th'] - o['th_true'])
    a.axhspan(-P.R_ENG, P.R_ENG, color=BAND, zorder=0, lw=0)
    a.axhline(P.R_ENG, color=INK2, lw=0.6, ls='--')
    a.axhline(-P.R_ENG, color=INK2, lw=0.6, ls='--')
    a.plot(t, err, color=C_TUR, lw=1.0, label='pointing error')
    for p in o['pulses']:
        k = int(round(p['t'] * f_ctrl))
        if p['hit']:
            a.plot(p['t'], err[k], 'o', ms=4.5, mfc=INK, mec='white', mew=0.8, zorder=5)
        else:
            a.plot(p['t'], np.clip(err[k], -6, 6), 'x', ms=5, color=INK, mew=1.2, zorder=5)
    a.set_ylim(-6, 6)
    a.set_ylabel('cross-range [m]')
    a.text(0.99, 0.95, '$\\pm R_{eng}$ band; dots = pulses', transform=a.transAxes,
           ha='right', va='top', color=INK2, fontsize=6.5)

    a = ax[2]
    _shade_blackouts(a, scenario)
    a.plot(t, o['V'], color=C_EST, lw=0.9)
    a.set_ylabel('$V_m$ [V]')
    a.set_ylim(-650, 650)

    a = ax[3]
    _shade_blackouts(a, scenario)
    a.plot(t, o['p_hit'], color=C_REF, lw=0.9)
    a.set_ylabel('$P_{hit}$')
    a.set_ylim(-0.05, 1.05)
    a.set_xlabel('time [s]')
    a.set_xlim(0, P.T_END)

    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'scenario_{scenario}.{ext}'))
    plt.close(fig)
    return m


def plot_trajectories(scenarios):
    from guardian.target_motion_generator import _get
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    th = np.linspace(np.radians(-10), np.radians(85), 100)
    for rr in (P.R_MIN, P.R_MAX):
        ax.plot(rr * np.cos(th), rr * np.sin(th), color=INK2, lw=0.6, ls='--')
    for s in scenarios:
        tr = _get(s)['traj']
        k = tr['t'] <= P.T_END
        ax.plot(tr['x'][k], tr['y'][k], color=C_EST if s in 'ABC' else INK2,
                lw=1.2 if s in 'ABC' else 0.7)
        ax.annotate(s, (tr['x'][0], tr['y'][0]), fontsize=6.5, color=INK,
                    xytext=(2, 2), textcoords='offset points')
    ax.plot(0, 0, marker='^', color=INK, ms=6)
    ax.set_aspect('equal')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title('Target trajectories (start labelled); turret at origin', loc='left')
    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'trajectories.{ext}'))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('scenarios', nargs='*', default=['A', 'B', 'C', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6'])
    ap.add_argument('--f-ctrl', type=float, default=1000.0)
    a = ap.parse_args()
    plot_trajectories(a.scenarios)
    for s in a.scenarios:
        m = plot_scenario(s, a.f_ctrl)
        print(f"{s}: first hit {m['t_first_hit']:.3f} s, {m['n_hits']}/{m['n_pulses']} hits, "
              f"pointing RMS {m['point_rms_m']:.2f} m, in-R {m['frac_in_Reng']:.3f}, Vpk {m['V_peak']:.0f} V")


if __name__ == '__main__':
    main()
