"""Analytical design checks for the paper: plant, feedback margins, reference generator.

python scripts/design_analysis.py
Prints the numbers quoted in the paper and writes figures to results/figures/.
"""
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scipy.signal as sg  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from guardian import params as P  # noqa: E402
from guardian import turret_controller as tc  # noqa: E402
from plot_scenarios import C_EST, C_REF, C_TUR, INK, INK2, OUT  # noqa: E402

NUM = [P.N_GEAR * P.K_T]
DEN = np.polymul([1.0, 0.0], np.polyadd(np.polymul([P.L_ARM, P.R_ARM], [P.J_EQ, P.B_TURRET]),
                                        [P.N_GEAR ** 2 * P.K_T * P.K_E]))


def plant_report():
    poles = np.roots(DEN)
    pc = poles[np.argmax(poles.imag)]
    wn = abs(pc)
    print('Plant theta/V = %g / (%s)' % (NUM[0], np.array2string(DEN, precision=4)))
    print('  poles', np.round(poles, 3), ' wn = %.2f rad/s (%.2f Hz), zeta = %.3f' % (wn, wn / 2 / np.pi, -pc.real / wn))
    h_max = min(-2 * p.real / abs(p) ** 2 for p in poles if abs(p) > 0)
    print('  forward-Euler stability limit dt < %.2f ms' % (1e3 * h_max))
    print('  plant inverse V = %.3f y\'\'\' + %.3f y\'\' + %.2f y\'' % (DEN[0] / NUM[0], DEN[1] / NUM[0], DEN[2] / NUM[0]))


def margins(L, w):
    mag = np.abs(L)
    ph = np.unwrap(np.angle(L))
    ic = np.where(np.diff(np.sign(mag - 1.0)))[0]
    wc = w[ic[0]]
    pm = 180.0 + np.degrees(ph[ic[0]])
    ip = np.where(np.diff(np.sign(ph + np.pi)))[0]
    gm = [(-20 * np.log10(mag[i]), w[i]) for i in ip]
    return wc, pm, gm


def loop_continuous(w):
    s = 1j * w
    kp, ti, tf = tc.CFG['kp'], tc.CFG['ti'], tc.CFG['tau_f']
    C = kp * (1 + 1 / (ti * s)) / (tf * s + 1)
    return C * np.polyval(NUM, s) / np.polyval(DEN, s)


def loop_discrete(w, fs):
    """Implemented C(z) (backward-Euler integral, exact first-order filter) with the ZOH plant."""
    T = 1.0 / fs
    numd, dend, _ = sg.cont2discrete((NUM, DEN), T, method='zoh')
    z = np.exp(1j * w * T)
    Pz = np.polyval(np.squeeze(numd), z) / np.polyval(dend, z)
    kp, ti, tf = tc.CFG['kp'], tc.CFG['ti'], tc.CFG['tau_f']
    a = math.exp(-T / tf)
    F = (1 - a) * z / (z - a)
    I = T * z / (z - 1)
    return F * kp * (1 + I / ti) * Pz


def feedback_report():
    print('\nFeedback C(s) = Kp (1 + 1/(Ti s)) / (tau_f s + 1): Kp=%g V/rad, Ti=%g s, tau_f=%g s'
          % (tc.CFG['kp'], tc.CFG['ti'], tc.CFG['tau_f']))
    rows = []
    for label, fs in (('continuous', None), ('1000 Hz', 1000.0), ('100 Hz', 100.0), ('50 Hz', 50.0), ('20 Hz', 20.0)):
        w = np.logspace(-2, np.log10(np.pi * fs * 0.999 if fs else 3000.0), 200000)
        L = loop_continuous(w) if fs is None else loop_discrete(w, fs)
        wc, pm, gm = margins(L, w)
        gmin = min(gm, key=lambda g: g[0]) if gm else (float('inf'), float('nan'))
        rows.append((label, wc, pm, gmin[0], gmin[1]))
        print('  %-10s  wc = %.2f rad/s  PM = %5.1f deg  GM = %5.1f dB at %.1f rad/s' % rows[-1])
    return rows


def bode_figure():
    w = np.logspace(-1, np.log10(3000), 4000)
    s = 1j * w
    Pw = np.polyval(NUM, s) / np.polyval(DEN, s)
    L = loop_continuous(w)
    L20 = loop_discrete(w[w < np.pi * 20], 20.0)
    fig, ax = plt.subplots(2, 1, figsize=(3.5, 3.6), sharex=True, gridspec_kw=dict(hspace=0.1))
    ax[0].semilogx(w, 20 * np.log10(np.abs(Pw)), color=INK2, lw=1.0, ls='--', label='plant $G_p$')
    ax[0].semilogx(w, 20 * np.log10(np.abs(L)), color=C_EST, lw=1.2, label='loop $G_cG_p$ (continuous)')
    ax[0].semilogx(w[w < np.pi * 20], 20 * np.log10(np.abs(L20)), color=C_REF, lw=1.0, ls=(0, (1, 1.5)),
                   label='loop, 20 Hz ZOH')
    ax[0].axhline(0, color=INK, lw=0.6)
    ax[0].set_ylabel('magnitude [dB]')
    ax[0].set_ylim(-160, 60)
    ax[0].legend(frameon=False, loc='lower left')
    ax[1].semilogx(w, np.degrees(np.unwrap(np.angle(L))), color=C_EST, lw=1.2)
    ax[1].semilogx(w[w < np.pi * 20], np.degrees(np.unwrap(np.angle(L20))), color=C_REF, lw=1.0, ls=(0, (1, 1.5)))
    ax[1].axhline(-180, color=INK, lw=0.6)
    ax[1].set_ylabel('phase [deg]')
    ax[1].set_xlabel('$\\omega$ [rad/s]')
    ax[1].set_ylim(-400, -60)
    wn = abs(np.roots(DEN)).max()
    for a in ax:
        a.axvline(wn, color=INK2, lw=0.6, ls=':')
    ax[0].text(wn * 1.1, 30, 'resonance\n%.1f rad/s' % wn, fontsize=6.5, color=INK2)
    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'bode.{ext}'))
    plt.close(fig)


def reference_report():
    tp = P.TP
    print('\nTaylor sections A_n(s) of exp(s tp), tp = %.2f s' % tp)
    for n in range(1, 7):
        c = [tp ** k / math.factorial(k) for k in range(n, -1, -1)]
        r = np.roots(c)
        z = min(-x.real / abs(x) for x in r)
        print('  n=%d  %-8s min damping %.3f  roots %s' % (n, 'Hurwitz' if np.all(r.real < 0) else 'UNSTABLE', z,
                                                         np.array2string(np.round(r, 2))))
    p = tc.CFG['p_dev']
    d3 = tp ** 3 / 6
    K = [d3 * p ** 3 - 1, 3 * d3 * p ** 2 - tp, 3 * d3 * p - tp ** 2 / 2]
    print('Deviation feedback K(s) = %.3f + %.3f s + %.4f s^2  (triple pole at -%g)' % (K[0], K[1], K[2], p))

    # deviation transient: unit initial deviation, A3 alone vs A3 + K
    t = np.linspace(0, 3, 600)
    a3 = [d3, tp ** 2 / 2, tp, 1.0]
    dev = [d3, d3 * 3 * p, d3 * 3 * p ** 2, d3 * p ** 3]
    fig, ax = plt.subplots(figsize=(3.5, 1.9))
    for den, col, lab in ((a3, C_REF, '$1/A_3(s)$ (prediction-matched lag only)'),
                          (dev, C_EST, '$1/(A_3+K)(s)$, triple pole $-%g$' % p)):
        sys_ = sg.lti([den[-1]], den)          # deviation decays like 1 - step
        _, ystep = sg.step(sys_, T=t)
        ax.plot(t, 1 - ystep, color=col, lw=1.2, label=lab)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xlabel('time [s]')
    ax.set_ylabel('deviation $\\delta(t)/\\delta(0)$')
    ax.legend(frameon=False, loc='upper right')
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'reference_deviation.{ext}'))
    plt.close(fig)


if __name__ == '__main__':
    plant_report()
    feedback_report()
    reference_report()
    bode_figure()
    print('\nfigures written to', os.path.abspath(OUT))
