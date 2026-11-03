"""Figures for the research chapter of the experimental report (docs/research_lemma_homotopy.md).

python scripts/plot_research.py
Needs results/research/ from research_contact_order.py, run_research.sh, run_research2.sh and
run_deconfound.sh (plus research_deconfound.py).
"""
import json
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
from plot_scenarios import C_EST, C_REF, C_TUR, INK, INK2, OUT  # noqa: E402

RES = os.path.join(os.path.dirname(__file__), '..', 'results', 'research')
TP = 0.5
C_4 = '#8a5cd0'
C_5 = '#c2477a'
DESIGNS = [  # tag, label, colour
    ('final', 'final', C_EST), ('blend0.2', r'blend $\alpha$=0.2', INK2), ('contact1', 'contact 1', C_TUR),
    ('contact2', 'contact 2', C_4), ('exact3', 'exact $n$=3', C_REF), ('mm', 'moment-matched', C_5),
    ('exact4', 'exact $n$=4', INK),
]


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'{name}.{ext}'))
    plt.close(fig)
    print('wrote', name)


def runs(tag, folder=RES):
    with open(os.path.join(folder, f'mc_{tag}.json')) as f:
        return json.load(f)['runs']


def totals(rs):
    pul = sum(r['n_pulses'] for r in rs)
    hits = sum(r['n_hits'] for r in rs)
    return dict(pulses=pul, hits=hits, misses=pul - hits, rate=(pul - hits) / pul,
                rms=float(np.median([r['point_rms_m'] for r in rs])),
                inR=float(np.mean([r['frac_in_Reng'] for r in rs])))


def fig_contact_order():
    with open(os.path.join(RES, 'contact_order.json')) as f:
        rows = json.load(f)['rows']
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw=dict(wspace=0.3))
    ms = sorted({r['m'] for r in rows})
    for fx, mk, lab in ((0, 'o', 'Taylor section ($f$=0)'), (None, 's', 'best completion, $f\\leq6$')):
        z = []
        for m in ms:
            sel = [r['zeta_min'] for r in rows if r['m'] == m and (fx is None or r['f'] == fx)]
            z.append(max(sel))
        a.plot(ms, z, marker=mk, color=C_EST if fx == 0 else C_REF, lw=1.2, label=lab)
    a.axhline(0, color=INK2, lw=0.6)
    a.axvspan(2.5, 4.5, color='#efeee9', zorder=0)
    a.set_xlabel('contact order $m$  ($Q - e^{x} = O(x^{m+1})$)')
    a.set_ylabel(r'min. root damping $\zeta_{\min}$')
    a.set_title('(a) Stable contact order is capped at 4--5')
    a.legend(frameon=False, loc='lower left')
    polys = {
        'exact $n$=3': [TP ** k / math.factorial(k) for k in range(4)],
        'exact $n$=4': [TP ** k / math.factorial(k) for k in range(5)],
        'moment-matched': [TP ** k / math.factorial(k) for k in range(4)] + [0.7414 / 24 * TP ** 4,
                                                                            0.4697 / 120 * TP ** 5],
        'contact 1': [1.0, 0.5, 1 / 12, 1 / 216],
        'contact 2': [1.0, 0.5, 0.125, 1 / 64, 1 / 1024],
        'final ($-20$, triple)': [1.0, 0.15, 0.0075, 0.000125],
    }
    cols = [C_REF, INK, C_5, C_TUR, C_4, C_EST]
    for (lab, p), c in zip(polys.items(), cols):
        r = np.roots(p[::-1])
        b.scatter(r.real, r.imag, s=18, color=c, label=lab, zorder=3)
    b.axvline(-6, color=INK2, lw=0.6, ls='--')
    b.text(-6.3, 6.6, r'$\sigma \leq N/t_p$', ha='right', fontsize=7, color=INK2)
    b.axvline(0, color=INK2, lw=0.6)
    b.set_xlim(-22, 1)
    b.set_xlabel('Re $s$ [rad/s]')
    b.set_ylabel('Im $s$ [rad/s]')
    b.set_title('(b) Deviation poles of each reference law')
    b.legend(frameon=False, fontsize=6, loc='upper left')
    save(fig, 'research_contact_order')


def fig_tradeoff():
    mc = {t: totals(runs(f'MC_{t}')) for t, _, _ in DESIGNS}
    hg = {t: totals(runs(f'HG_{t}')) for t, _, _ in DESIGNS}
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw=dict(wspace=0.42))
    offs = {'contact1': (-4, 6), 'contact2': (-6, -11), 'exact3': (4, -10), 'mm': (4, 4), 'exact4': (-40, 6)}
    for t, lab, c in DESIGNS:
        a.scatter(mc[t]['rms'], 100 * hg[t]['rate'], color=c, s=26, zorder=3)
        a.annotate(lab, (mc[t]['rms'], 100 * hg[t]['rate']), xytext=offs.get(t, (4, 3)),
                   textcoords='offset points', fontsize=6.5, color=c)
    a.set_xlabel('random engagements: pointing RMS median [m]')
    a.set_ylabel('adversarial: misses per pulse [%]')
    a.set_title('(a) Slow poles: fewer misses per shot, worse pointing')
    x = np.arange(len(DESIGNS))
    w = 0.38
    a_mc = [mc[t]['hits'] for t, _, _ in DESIGNS]
    a_hg = [hg[t]['hits'] for t, _, _ in DESIGNS]
    b.bar(x - w / 2, a_mc, w, color=C_EST, label='MC0--999')
    b.bar(x + w / 2, a_hg, w, color=C_REF, label='HG0--999')
    b.set_xticks(x, [lab for _, lab, _ in DESIGNS], rotation=35, ha='right', fontsize=6.5)
    b.set_ylim(7000, 9400)
    b.set_ylabel('hits delivered (1000 engagements)')
    b.set_title('(b) ...but the final design delivers the most hits')
    b.legend(frameon=False, loc='upper right')
    save(fig, 'research_tradeoff')


def fig_deconfound():
    with open(os.path.join(RES, 'deconfound', 'deconfound.json')) as f:
        M = json.load(f)['matrix']
    names = ['final', 'contact1', 'exact3']
    labels = ['final', 'contact 1', 'exact $n$=3']
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw=dict(wspace=0.45))
    for a, key, title, fmt, cmap in (
            (axs[0], 'miss_rate', '(a) misses per pulse [%]', lambda v: f'{100 * v:.2f}', 'Oranges'),
            (axs[1], 'hits', '(b) hits delivered (200 engagements)', lambda v: f'{v:.0f}', 'Blues')):
        Z = np.zeros((3, 3))
        for i, att in enumerate(names):
            for j, dfn in enumerate(names):
                c = M[f'{att}|{dfn}']
                Z[i, j] = c['miss_rate'] if key == 'miss_rate' else c['pulses'] - c['misses']
        a.imshow(Z, cmap=cmap, aspect='auto')
        a.grid(False)
        mid = 0.5 * (Z.min() + Z.max())
        for i in range(3):
            for j in range(3):
                a.text(j, i, fmt(Z[i, j]), ha='center', va='center', fontsize=7.5,
                       color='white' if Z[i, j] > mid else INK,
                       fontweight='bold' if i == j else 'normal')
        a.set_xticks(range(3), labels)
        a.set_yticks(range(3), labels)
        a.set_xlabel('design under attack')
        a.set_ylabel('adversary tuned against')
        a.set_title(title)
    save(fig, 'research_deconfound')


def main():
    fig_contact_order()
    fig_tradeoff()
    fig_deconfound()


if __name__ == '__main__':
    main()
