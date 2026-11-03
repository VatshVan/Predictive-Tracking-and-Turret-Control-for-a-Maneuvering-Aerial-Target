"""Contact-order study: how well can a Hurwitz polynomial agree with e^x at x = 0?

python scripts/research_contact_order.py
Writes results/research/contact_order.json and results/research/contact_order.md.

Q(x) = sum_k q_k x^k with q_k = 1/k! for k <= m (contact order m with e^x: Q - e^x = O(x^(m+1)))
and f = 0..6 free higher coefficients q_{m+1..m+f}. For each (m, f) the free coefficients are
chosen to maximise the minimum root damping  zeta_min = min_i (-Re x_i / |x_i|).

Necessary condition, independent of the free part: Q is Hurwitz iff its reciprocal
x^N Q(1/x) is, whose leading coefficients are q_0, q_1, ...; the Hurwitz minor Delta_j of the
reciprocal only involves q_0..q_(2j-1). So Delta_j with 2j-1 <= m are fixed by the contact order:
Delta_1 = 1, Delta_2 = 1/3, Delta_3 = 1/45, ... If one of them is <= 0, no completion of any degree
is Hurwitz.

Also: exact-prediction reference models 1/A_n(tp s) for n = 3, 4, 5 (poles, damping, step
deviation response), showing n = 5 is unstable.
"""
import json
import math
import os
import sys
from fractions import Fraction

import numpy as np
from scipy.optimize import differential_evolution

RESULTS = os.path.join(os.path.dirname(__file__), '..', 'results', 'research')
TP = 0.5


def hurwitz_minors(desc):
    """Leading principal minors of the Hurwitz matrix of a0 s^N + a1 s^(N-1) + ... (exact)."""
    a = list(desc)
    N = len(a) - 1
    H = [[Fraction(0)] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            k = 2 * (j + 1) - (i + 1)
            if 0 <= k <= N:
                H[i][j] = Fraction(a[k])
    out = []
    for n in range(1, N + 1):
        out.append(_det([row[:n] for row in H[:n]]))
    return out


def _det(M):
    M = [row[:] for row in M]
    n = len(M)
    d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n):
                M[r][k] -= f * M[c][k]
    return d


def fixed_minors(m):
    """Hurwitz minors of the reciprocal polynomial fixed by contact order m."""
    q = [Fraction(1, math.factorial(k)) for k in range(m + 1)]
    out = []
    for j in range(1, m + 2):
        if 2 * j - 1 > m:
            break
        # minors of a degree-(2j-1) truncation equal those of any longer polynomial for this j
        out.append(hurwitz_minors(q[:2 * j])[j - 1])
    return out


def zeta_min(asc):
    r = np.roots(asc[::-1])
    return float(np.min(-r.real / np.abs(r)))


def best_completion(m, f, restarts=3):
    """Differential evolution over log-scaled free coefficients, best of several restarts."""
    base = [1.0 / math.factorial(k) for k in range(m + 1)]
    if f == 0:
        return base, zeta_min(base)
    scale = np.array([1.0 / math.factorial(m + j) for j in range(1, f + 1)])

    def cost(u):
        asc = base + list(scale * np.exp(u))
        return -zeta_min(asc)

    best = None
    for seed in range(restarts):
        res = differential_evolution(cost, [(-20.0, 8.0)] * f, seed=seed, tol=1e-10, maxiter=500,
                                     popsize=30, polish=True)
        if best is None or res.fun < best.fun:
            best = res
    asc = base + list(scale * np.exp(best.x))
    return asc, zeta_min(asc)


def ref_model(n):
    asc = [TP ** k / math.factorial(k) for k in range(n + 1)]
    poles = np.roots(asc[::-1])
    # unit initial deviation of the reference-model output: free response of 1/A_n
    A = np.zeros((n, n))
    A[:-1, 1:] = np.eye(n - 1)
    A[-1, :] = -np.array(asc[:n]) / asc[n]
    t = np.linspace(0, 3.0, 301)
    from scipy.linalg import expm
    x0 = np.zeros(n)
    x0[0] = 1.0
    y = np.array([(expm(A * tt) @ x0)[0] for tt in t])
    return dict(n=n, poles=[[float(p.real), float(p.imag)] for p in poles],
                zeta_min=float(np.min(-poles.real / np.abs(poles))),
                slowest_decay_s=float(-1.0 / np.max(poles.real)) if np.max(poles.real) < 0 else None,
                hurwitz=bool(np.all(poles.real < 0)),
                dev_abs_at_1s=float(abs(y[100])), dev_abs_at_3s=float(abs(y[-1])),
                dev_peak_after_0p5s=float(np.max(np.abs(y[50:]))))


def main():
    os.makedirs(RESULTS, exist_ok=True)
    rows = []
    for m in range(3, 9):
        fm = fixed_minors(m)
        for f in range(7):
            asc, z = best_completion(m, f)
            rows.append(dict(m=m, f=f, degree=m + f, zeta_min=z, hurwitz=z > 0,
                             free_coeffs=asc[m + 1:],
                             free_over_taylor=[c * math.factorial(m + 1 + j) for j, c in enumerate(asc[m + 1:])]))
        rows[-1]['fixed_minors'] = [str(x) for x in fm]
    taylor = {n: bool(np.all(np.roots([1.0 / math.factorial(k) for k in range(n, -1, -1)]).real < 0))
              for n in range(1, 11)}
    refs = [ref_model(n) for n in (3, 4, 5)]
    out = dict(rows=rows, taylor_hurwitz=taylor, fixed_minors={m: [str(x) for x in fixed_minors(m)] for m in range(1, 13)},
               ref_models=refs)
    with open(os.path.join(RESULTS, 'contact_order.json'), 'w') as fh:
        json.dump(out, fh, indent=1)

    L = ['# Contact-order study', '',
         'Pure Taylor sections T_n(x) Hurwitz: ' + ', '.join(f'n={n}: {"yes" if v else "no"}' for n, v in taylor.items()), '',
         '## Hurwitz minors fixed by the contact order (reciprocal polynomial)', '',
         '| m | fixed minors Delta_1..Delta_j (2j-1 <= m) |', '|---|---|']
    for m in range(1, 13):
        L.append(f"| {m} | {', '.join(out['fixed_minors'][m])} |")
    L += ['', '## Best minimum root damping with f free higher coefficients', '',
          '| contact order m | f=0 (Taylor) | f=1 | f=2 | f=3 | f=4 | f=5 | f=6 |', '|---|---|---|---|---|---|---|---|']
    for m in range(3, 9):
        zs = [r['zeta_min'] for r in rows if r['m'] == m]
        L.append(f'| {m} | ' + ' | '.join(f'{z:+.4f}' for z in zs) + ' |')
    L += ['', 'zeta_min > 0 means a Hurwitz completion exists. Free coefficients (as multiples of the Taylor value '
          '1/k!) for the best completions:', '']
    for r in rows:
        if r['f'] > 0:
            L.append(f"- m={r['m']}, f={r['f']}: zeta_min={r['zeta_min']:+.4f}, q/q_taylor = "
                     + ', '.join(f'{c:.4g}' for c in r['free_over_taylor']))
    L += ['', '## Exact-prediction reference models 1/A_n(tp s), tp = 0.5 s', '',
          '| n | poles [rad/s] | zeta_min | slowest time constant [s] | |deviation| at 1 s | at 3 s |',
          '|---|---|---|---|---|---|']
    for r in refs:
        poles = ', '.join(f'{p[0]:.2f}{p[1]:+.2f}j' for p in r['poles'])
        tc_ = f"{r['slowest_decay_s']:.3f}" if r['slowest_decay_s'] else 'unstable'
        L.append(f"| {r['n']} | {poles} | {r['zeta_min']:+.3f} | {tc_} | {r['dev_abs_at_1s']:.3g} | {r['dev_abs_at_3s']:.3g} |")
    with open(os.path.join(RESULTS, 'contact_order.md'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L) + '\n')
    sys.stdout.reconfigure(encoding='utf-8')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
