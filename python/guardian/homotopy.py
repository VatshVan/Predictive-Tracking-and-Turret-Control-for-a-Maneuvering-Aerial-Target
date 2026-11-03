"""Beam-relative homotopy classes and guided adversarial engagements ('HG<seed>' scenarios).

Homotopy classes (T-MPC mapping, de Groot et al. 2024, Def. 1 and App. A). Work in the
state-time plane (t, e), e = signed cross-range offset of the target from the beam. The
obstacles are the shooter's firing events O_k = {t_k} x [-R_ENG, R_ENG]. A target path is a graph
e(t), so its winding number around each O_k is fixed by the side on which it passes the event,
and two paths are homotopic iff they pass every event on the same side. The class of a path is
its H-signature h_k in {+1 above, 0 through (= hit), -1 below}.

Guided adversary. Within the admissible envelope, after a blackout the beam follows the radar,
so a target can leave the tube at a firing event only through estimator lag. HG<seed> realizes
the class seed % 5 with no 'through' events:
  above      (+, +, ..., +)          below    (-, -, ..., -)
  cross_down (+, ..., +, -, ..., -)  cross_up (-, ..., -, +, ..., +)
  multiple   (+, -, +, -, ...)
The final design's firing events are found by running it on the base trajectory. Before each
event k, the normal-to-LOS acceleration reverses (quintic, finite jerk) toward the side h_k at
the largest admissible magnitude, a lead delta ahead of the event. delta is chosen from a grid to
maximise the mean miss of the final design, with one fixed-point re-run on the new events.
Trajectories violating the operating envelope are rejected and resampled. Cases are cached in
results/hg_cache/ so every design variant faces identical trajectories.
"""
import os

import numpy as np

from . import params as P

CLASSES = ('above', 'below', 'cross_down', 'cross_up', 'multiple')
LEADS = (0.10, 0.20, 0.35, 0.50, 0.70)
CACHE = os.path.join(os.path.dirname(__file__), '..', 'results', 'hg_cache')
# designs the adversary can be tuned against: 'HG<seed>' targets the final design,
# 'HX<target>_<seed>' the named design (same base engagement and radar noise as HG<seed>)
TARGETS = {
    'final': {},
    'contact1': dict(ref_mode='contact'),
    'exact3': dict(ref_mode='exact_taylor', ref_n=3),
}


def signature(offsets, r_eng=P.R_ENG):
    """H-signature of offsets at the firing events: +1 above, 0 through, -1 below."""
    e = np.asarray(offsets, dtype=float)
    return np.where(e > r_eng, 1, np.where(e < -r_eng, -1, 0)).astype(int)


def class_signs(cls, n_ev):
    """Intended side at each of n_ev events for class index cls."""
    if cls == 0:
        return [1] * n_ev
    if cls == 1:
        return [-1] * n_ev
    h = n_ev // 2
    if cls == 2:
        return [1] * h + [-1] * (n_ev - h)
    if cls == 3:
        return [-1] * h + [1] * (n_ev - h)
    return [1 if k % 2 == 0 else -1 for k in range(n_ev)]


def classify(sig):
    """Class index of a realised signature (events passed 'through' are skipped); -1 if all hits."""
    s = [v for v in sig if v != 0]
    if not s:
        return -1
    n_cross = sum(1 for a, b in zip(s, s[1:]) if a != b)
    if n_cross == 0:
        return 0 if s[0] > 0 else 1
    if n_cross == 1:
        return 2 if s[0] > 0 else 3
    return 4


def _cap(a0, u):
    """Largest lam >= 0 with |a0 + lam u| <= 0.97 A_MAX (u unit vector)."""
    amax = 0.97 * P.A_MAX_TGT
    au = a0 @ u
    return max(-au + np.sqrt(max(au * au - a0 @ a0 + amax ** 2, 0.0)), 0.0)


def _schedule(t, events, signs, lead, ramp=0.10, hold=0.05):
    """Sign pattern s(t): -h_k until event k - lead, then +h_k until event k + hold."""
    from .target_motion_generator import _smooth
    s = np.zeros_like(t)
    prev_end = events[0] - 1.0
    for tk, hk in zip(events, signs):
        t_rev = tk - lead
        pre = (t >= prev_end) & (t < t_rev)
        s[pre] = -hk
        post = (t >= t_rev) & (t < tk + hold)
        s[post] = hk
        prev_end = tk + hold
    # quintic smoothing of the steps (finite jerk)
    k = max(1, int(round(ramp / (t[1] - t[0]))))
    w = _smooth(np.linspace(0.0, 1.0, k + 2))
    w = np.diff(w)
    return np.convolve(s, w / w.sum(), mode='same')


def _base(rng, t):
    """Random base engagement: mild maneuver, one or two blackouts."""
    from .target_motion_generator import _smooth
    r0 = rng.uniform(1500.0, 2000.0)
    th0 = np.deg2rad(rng.uniform(8.0, 55.0))
    speed = rng.uniform(70.0, 130.0)
    heading = np.pi + th0 + np.deg2rad(rng.uniform(-40.0, 40.0))
    p0 = (r0 * np.cos(th0), r0 * np.sin(th0))
    v0 = (speed * np.cos(heading), speed * np.sin(heading))
    f = rng.uniform(0.1, 0.6, (2, 2))
    phi = rng.uniform(0, 2 * np.pi, (2, 2))
    amp = rng.uniform(0.0, 0.15 * P.A_MAX_TGT, (2, 2))
    g = _smooth((t - 0.5) / 0.5)
    ax = g * np.sum(amp[0][:, None] * np.sin(2 * np.pi * f[0][:, None] * t + phi[0][:, None]), 0)
    ay = g * np.sum(amp[1][:, None] * np.sin(2 * np.pi * f[1][:, None] * t + phi[1][:, None]), 0)
    blackouts = [(rng.uniform(P.T_BOO_MIN, 8.5), rng.uniform(0.1, P.T_BO_MAX))]
    if rng.random() < 0.3:
        t2 = rng.uniform(blackouts[0][0] + blackouts[0][1] + 0.5, 9.5)
        if t2 < 9.5:
            blackouts.append((t2, rng.uniform(0.1, P.T_BO_MAX)))
    return dict(p0=p0, v0=v0, blackouts=blackouts, accel_grid=(ax, ay))


def _simulate(d, target='final'):
    """Target design on definition d; returns the pulse list (fire events with signed offsets)."""
    from . import target_motion_generator as tmg
    from . import turret_controller as tc
    from .sim import run
    key = '__hg_probe__'
    saved = dict(tc.CFG)
    tc.CFG.clear()
    tc.CFG.update(tc.DEFAULT_CFG)
    tc.CFG.update(TARGETS[target])
    tmg._CACHE.pop(key, None)
    tmg._PROBE = d
    try:
        return run(key, f_ctrl=100.0, log=False)['pulses']
    finally:
        tmg._PROBE = None
        tmg._CACHE.pop(key, None)
        tc.CFG.clear()
        tc.CFG.update(saved)
        tc.reset()


def _attack(base, t, events, cls, lead):
    from .target_motion_generator import _integrate
    traj0 = _integrate(base, t)
    th = np.arctan2(traj0['y'], traj0['x'])
    nx, ny = -np.sin(th), np.cos(th)                      # LOS normal along the base path
    s = _schedule(t, events, class_signs(cls, len(events)), lead)
    ax, ay = base['accel_grid']
    lam = np.array([_cap(np.array([ax[i], ay[i]]), np.array([nx[i], ny[i]]) * (1 if s[i] >= 0 else -1))
                    for i in range(0, len(t), 50)])
    lam = np.interp(t, t[::50], lam)
    d = dict(base)
    d['accel_grid'] = (ax + lam * s * nx, ay + lam * s * ny)
    return d


def generate(seed, target='final'):
    """Build HG<seed> by closed-loop adversarial search against the target design."""
    from .target_motion_generator import DT_GRID, T_GRID_END, _envelope_ok, _integrate
    rng = np.random.default_rng(900001 + seed)
    cls = seed % 5
    t = np.arange(0.0, T_GRID_END + DT_GRID / 2, DT_GRID)
    for _ in range(200):
        base = _base(rng, t)
        traj = _integrate(base, t)
        if not _envelope_ok(traj, elev_margin=True):
            continue
        events = [p['t'] for p in _simulate(base, target)]
        if len(events) < 2:
            continue
        best = None
        for lead in LEADS:
            d = _attack(base, t, events, cls, lead)
            if not _envelope_ok(_integrate(d, t), elev_margin=True):
                continue
            pulses = _simulate(d, target)
            score = float(np.mean([p['miss'] for p in pulses])) if pulses else 0.0
            if best is None or score > best[0]:
                best = (score, lead, d)
        if best is None:
            continue
        # fixed point: re-aim at the events the attacked engagement actually produced
        ev1 = [p['t'] for p in _simulate(best[2], target)]
        d = _attack(base, t, ev1, cls, best[1]) if len(ev1) >= 2 else best[2]
        traj = _integrate(d, t)
        if not _envelope_ok(traj, elev_margin=True):
            d, traj = best[2], _integrate(best[2], t)
        pulses = _simulate(d, target)
        sig = signature([p['offset'] for p in pulses])
        return dict(ax=d['accel_grid'][0], ay=d['accel_grid'][1], p0=np.array(d['p0']), v0=np.array(d['v0']),
                    blackouts=np.array(d['blackouts']), cls=cls, lead=best[1],
                    events=np.array(ev1), signature=sig, realised=classify(sig),
                    miss_mean=float(np.mean([p['miss'] for p in pulses])) if pulses else 0.0)
    raise RuntimeError(f"HG seed {seed}: no admissible adversarial trajectory found")


def hg_def(seed, target='final'):
    """Scenario definition for HG<seed> / HX<target>_<seed> (generated once, then cached)."""
    folder = CACHE if target == 'final' else os.path.join(CACHE, target)
    path = os.path.join(folder, f'HG{seed}.npz')
    if os.path.exists(path):
        z = np.load(path)
        g = {k: z[k] for k in z.files}
    else:
        g = generate(seed, target)
        os.makedirs(folder, exist_ok=True)
        tmp = path + f'.{os.getpid()}.tmp.npz'
        np.savez(tmp, **g)
        os.replace(tmp, path)
    d = dict(p0=tuple(g['p0']), v0=tuple(g['v0']), blackouts=[tuple(b) for b in np.atleast_2d(g['blackouts'])],
             accel_grid=(g['ax'], g['ay']))
    d['hg'] = dict(cls=int(g['cls']), name=CLASSES[int(g['cls'])], lead=float(g['lead']),
                   signature=[int(v) for v in g['signature']], realised=int(g['realised']),
                   miss_mean=float(g['miss_mean']))
    return d

