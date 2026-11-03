"""Pole-shift study: how far left can the reference-model poles go?

python scripts/research_pole_shift.py
Writes results/research/pole_shift.json and results/research/pole_shift.md.

In the linear band every reference mode reduces to deviation dynamics P(D)(y_d - theta_hat) = 0,
excited only by the estimator's correction steps, so the poles of P set how fast corrections
are absorbed.

Decay bound. If P(s) = sum p_k s^k with p_0 = 1 and p_1 = tp (first moment of e^{tp s}), then
sum_i 1/s_i = -p_1/p_0 = -tp. A pole with decay rate sigma = -Re s_i satisfies |Re(1/s_i)| <= 1/sigma,
so tp <= N / sigma_min, i.e. sigma_min <= N / tp. A degree-N loop that keeps even first-order
prediction fidelity cannot have every pole faster than N / tp (6 rad/s for N = 3, tp = 0.5 s).

Part 1: best decay with contact order 3 kept and f free coefficients (degree 3 + f).
Part 2: weighted average P_a = (1 - a) A_3 + a d3 (s + p)^3, the 'blend' reference mode.
"""
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import differential_evolution

RESULTS = os.path.join(os.path.dirname(__file__), '..', 'results', 'research')
TP, P_DEV = 0.5, 20.0
A3 = [TP ** k / math.factorial(k) for k in range(4)]


def poles_info(asc):
    r = np.roots(np.array(asc)[::-1])
    return dict(decay=float(-r.real.max()), zeta=float((-r.real / np.abs(r)).min()),
                poles=[[float(p.real), float(p.imag)] for p in r])


def best_decay(f, zeta_floor=0.3, restarts=3):
    """Maximise the slowest decay rate with q_0..q_3 fixed (contact order 3) and zeta >= floor."""
    sc = np.array([TP ** (3 + j) / math.factorial(3 + j) for j in range(1, f + 1)])

    def cost(u):
        i = poles_info(A3 + list(sc * np.exp(u)))
        return -i['decay'] + 50.0 * max(0.0, zeta_floor - i['zeta'])

    best = min((differential_evolution(cost, [(-15.0, 8.0)] * f, seed=s, maxiter=400, popsize=30, tol=1e-10)
                for s in range(restarts)), key=lambda r: r.fun)
    asc = A3 + list(sc * np.exp(best.x))
    return dict(f=f, degree=3 + f, bound=(3 + f) / TP, **poles_info(asc))


def main():
    os.makedirs(RESULTS, exist_ok=True)
    contact = [dict(f=0, degree=3, bound=3 / TP, **poles_info(A3))] + [best_decay(f) for f in range(2, 7)]
    d3 = A3[3]
    D_dev = [d3 * c for c in (P_DEV ** 3, 3 * P_DEV ** 2, 3 * P_DEV, 1.0)]
    blends = []
    for a in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0):
        P = [(1 - a) * q + a * d for q, d in zip(A3, D_dev)]
        blends.append(dict(a=a, p0=P[0], p1_over_p0=P[1] / P[0], **poles_info(P)))
    with open(os.path.join(RESULTS, 'pole_shift.json'), 'w') as f:
        json.dump(dict(contact=contact, blends=blends), f, indent=1)
    L = ['# Pole-shift study', '', f'Decay bound with first-moment fidelity: sigma_min <= N / tp = 2N rad/s.', '',
         '## Contact order 3 kept, f free coefficients (zeta >= 0.3)', '',
         '| degree N | bound N/tp [rad/s] | best decay [rad/s] | zeta_min |', '|---|---|---|---|']
    for c in contact:
        L.append(f"| {c['degree']} | {c['bound']:.0f} | {c['decay']:.2f} | {c['zeta']:.3f} |")
    L += ['', '## Weighted average P_a = (1 - a) A_3 + a d3 (s + 20)^3', '',
          '| a | slowest decay [rad/s] | zeta_min | P(0) | P1/P0 (tp = 0.5 is exact) |', '|---|---|---|---|---|']
    for b in blends:
        L.append(f"| {b['a']} | {b['decay']:.2f} | {b['zeta']:.3f} | {b['p0']:.1f} | {b['p1_over_p0']:.3f} |")
    with open(os.path.join(RESULTS, 'pole_shift.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')
    sys.stdout.reconfigure(encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
