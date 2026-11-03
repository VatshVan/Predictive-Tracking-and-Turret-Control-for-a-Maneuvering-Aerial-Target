"""Attacker x defender matrix for the closed-loop HG adversary (scripts/run_deconfound.sh).

python scripts/research_deconfound.py
Writes results/research/deconfound/deconfound_report.md and deconfound.json.

Rows: the design the adversary was tuned against (final = HG0-199 from the main campaign).
Columns: the design under attack. Same 200 base engagements and radar noise in every cell, so
cells in one row are paired (same trajectories), and the diagonal is each design under an
adversary tuned against itself.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from mc_statistics import holm, paired  # noqa: E402

RES = os.path.join(os.path.dirname(__file__), '..', 'results', 'research')
OUT = os.path.join(RES, 'deconfound')
N = 200
DESIGNS = ('final', 'contact1', 'exact3')


def load(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        runs = [r for r in json.load(f)['runs'] if r['seed'] < N]
    return sorted(runs, key=lambda r: r['seed'])


def cell(attacker, design):
    if attacker == 'final':
        return load(os.path.join(RES, f'mc_HG_{design}.json'))
    return load(os.path.join(OUT, f'mc_HX{attacker}_{design}.json'))


def summ(runs):
    pul = sum(r['n_pulses'] for r in runs)
    mis = pul - sum(r['n_hits'] for r in runs)
    return dict(misses=mis, pulses=pul, miss_rate=mis / pul, succ=sum(r['success'] for r in runs),
                miss_mean=float(np.mean([p[1] for r in runs for p in r['pulses']])),
                rms_median=float(np.median([r['point_rms_m'] for r in runs])))


def main():
    M, C = {}, {}
    for a in DESIGNS:
        for d in DESIGNS:
            runs = cell(a, d)
            if runs is not None:
                M[(a, d)] = summ(runs)
        base = cell(a, 'final')
        for d in DESIGNS[1:]:
            v = cell(a, d)
            if base is not None and v is not None:
                C[(a, d)] = paired(base, v)
    adj = holm({k: c['p_miss'] for k, c in C.items()})
    L = ['# Attacker x defender matrix (HG adversary, seeds 0-199)', '',
         'Rows: design the adversary was tuned against. Columns: design under attack. Cells: missed pulses / '
         'pulses (miss rate), mean miss distance, interceptions. Diagonal = each design against its own adversary.', '',
         '| adversary tuned on ↓ / design → | ' + ' | '.join(DESIGNS) + ' |', '|---|' + '---|' * len(DESIGNS)]
    for a in DESIGNS:
        row = []
        for d in DESIGNS:
            m = M.get((a, d))
            row.append('—' if m is None else
                       f"{m['misses']}/{m['pulses']} ({100 * m['miss_rate']:.2f} %), {m['miss_mean']:.2f} m, "
                       f"{m['succ']}/{N}")
        L.append(f'| {a} | ' + ' | '.join(row) + ' |')
    L += ['', '## Paired with the final design under the same adversary (Holm over all cells)', '',
          '| adversary | design | misses final → design | p (Holm) | Δ first hit [s] | Δ RMS [m] |', '|---|---|---|---|---|---|']
    for (a, d), c in C.items():
        L.append(f"| {a} | {d} | {c['misses_base']} → {c['misses_var']} | {adj[(a, d)]:.2g} | "
                 f"{c['d_fh_mean']:+.3f} | {c['d_rms_mean']:+.3f} |")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'deconfound.json'), 'w') as f:
        json.dump(dict(matrix={f'{a}|{d}': v for (a, d), v in M.items()},
                       paired={f'{a}|{d}': v for (a, d), v in C.items()}), f, indent=1)
    with open(os.path.join(OUT, 'deconfound_report.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')
    sys.stdout.reconfigure(encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
