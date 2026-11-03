"""Target trajectories and radar measurements for the standardized and stress scenarios.

Scenario identifiers
  'A'   nominal non-maneuvering target                         (spec Scenario A)
  'B'   smooth bounded maneuver in x and y from t = 1 s          (spec Scenario B)
  'C'   nominal target, radar loss on [1.0, 1.5) s                (spec Scenario C)
  'D1'  stealth break: maneuver begins inside a late blackout
  'D2'  weave: sinusoidal near-a_max jinking in x and y
  'D3'  close-in pass: r -> 800 m, elevation approaching +85 deg
  'D4'  repeated blackouts on the Scenario B maneuver
  'D5'  combined: weave + late blackout + close range
  'D6'  jerk test: near bang-bang acceleration reversals
  'MC<seed>'  random admissible trajectory and blackout pattern (e.g. 'MC17')
  'HG<seed>'  closed-loop adversary in homotopy class seed % 5 against the final design
  'HX<target>_<seed>'  the same adversary tuned against another design (homotopy.TARGETS)

Every trajectory is generated once on a 1 ms grid, checked against the operating
envelope (800 <= r <= 2000 m, |v| <= 150 m/s, |a| <= 20 m/s^2) and cached.
Measurements are produced at radar ticks t = k * 0.05 s; measurement noise is indexed
by the tick number so results do not depend on the caller's step size.
"""
import numpy as np

from . import params as P

DT_GRID = 1e-3
T_GRID_END = P.T_END + 1.0

_CACHE = {}
_PROBE = None


def _smooth(tau):
    """Quintic smoothstep, C2-continuous, 0 for tau <= 0 and 1 for tau >= 1."""
    s = np.clip(tau, 0.0, 1.0)
    return s**3 * (10.0 - 15.0 * s + 6.0 * s**2)


def _window(t, t_on, t_off, ramp):
    """Smooth pulse: rises over [t_on, t_on+ramp], falls over [t_off, t_off+ramp]."""
    return _smooth((t - t_on) / ramp) - _smooth((t - t_off) / ramp)


def _scenario_def(scenario):
    """Initial state, acceleration profile and blackout intervals for a named scenario."""
    if scenario in ('A', 'C'):
        d = dict(p0=(1700.0, 700.0), v0=(-120.0, 0.0),
                 accel=lambda t: (0.0 * t, 0.0 * t))
        d['blackouts'] = [(1.0, 0.5)] if scenario == 'C' else []
        return d
    if scenario in ('B', 'D4'):
        def accel(t):
            g = _smooth((t - 1.0) / 0.5)
            return (g * 8.0 * np.sin(2 * np.pi * 0.15 * (t - 1.0)),
                    g * 14.0 * np.sin(2 * np.pi * 0.30 * (t - 1.0)))
        d = dict(p0=(1700.0, 700.0), v0=(-120.0, 0.0), accel=accel)
        d['blackouts'] = ([(1.0, 0.5), (3.2, 0.4), (5.5, 0.5), (7.5, 0.4)]
                          if scenario == 'D4' else [])
        return d
    if scenario == 'D1':
        def accel(t):
            w = _window(t, 3.1, 4.4, 0.3)
            return (w * 9.0, w * 17.0)
        return dict(p0=(1900.0, 500.0), v0=(-125.0, 0.0), accel=accel,
                    blackouts=[(3.0, 0.5)])
    if scenario in ('D2', 'D5'):
        def accel(t):
            g = _smooth((t - 0.5) / 0.5)
            return (g * 13.0 * np.sin(2 * np.pi * t / 1.6),
                    g * 14.0 * np.cos(2 * np.pi * t / 1.2))
        if scenario == 'D2':
            return dict(p0=(1800.0, 800.0), v0=(-110.0, 0.0), accel=accel, blackouts=[])
        return dict(p0=(1300.0, 700.0), v0=(-80.0, 0.0), accel=accel,
                    blackouts=[(4.0, 0.5)])
    if scenario == 'D3':
        def accel(t):
            g = _smooth((t - 1.0) / 0.5)
            return (0.0 * t, g * 4.0 * np.sin(2 * np.pi * 0.4 * (t - 1.0)))
        return dict(p0=(1450.0, 835.0), v0=(-135.0, 0.0), accel=accel, blackouts=[])
    if scenario == 'D6':
        def accel(t):
            # +-18 m/s^2 reversals every 1.2 s with 0.15 s transitions (peak jerk ~450 m/s^3)
            s = np.zeros_like(t)
            edges = 1.0 + 0.6 + 1.2 * np.arange(8)
            s += _smooth((t - 1.0) / 0.15)
            for k, te in enumerate(edges):
                s += (-2.0 if k % 2 == 0 else 2.0) * _smooth((t - te) / 0.15)
            return (0.0 * t, 18.0 * s)
        return dict(p0=(1750.0, 850.0), v0=(-125.0, 0.0), accel=accel, blackouts=[])
    if isinstance(scenario, str) and scenario.startswith('MC'):
        return _monte_carlo_def(int(scenario[2:]))
    if scenario == '__hg_probe__':
        return dict(_PROBE)                     # candidate trajectory of the HG adversary search
    if isinstance(scenario, str) and scenario.startswith('HX'):
        from .homotopy import hg_def
        target, _, seed = scenario[2:].rpartition('_')
        return hg_def(int(seed), target)
    if isinstance(scenario, str) and scenario.startswith('HG'):
        from .homotopy import hg_def
        return hg_def(int(scenario[2:]))
    raise ValueError(f"unknown scenario {scenario!r}")


def _monte_carlo_def(seed):
    """Random admissible engagement: smooth bounded maneuvers and random radar loss."""
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, T_GRID_END + DT_GRID / 2, DT_GRID)
    for _ in range(500):
        r0 = rng.uniform(1500.0, 2000.0)
        th0 = np.deg2rad(rng.uniform(8.0, 55.0))
        speed = rng.uniform(70.0, 150.0)
        heading = np.pi + th0 + np.deg2rad(rng.uniform(-40.0, 40.0))
        p0 = (r0 * np.cos(th0), r0 * np.sin(th0))
        v0 = (speed * np.cos(heading), speed * np.sin(heading))

        n_sin = 3
        f = rng.uniform(0.1, 1.0, (2, n_sin))
        phi = rng.uniform(0, 2 * np.pi, (2, n_sin))
        amp = rng.uniform(0.2, 1.0, (2, n_sin))
        t_on = rng.uniform(0.5, 3.0)
        has_break = rng.random() < 0.5
        tb = rng.uniform(1.0, 7.0)
        tb_len = rng.uniform(0.5, 2.0)
        tb_ramp = rng.uniform(0.15, 0.6)
        bdir = rng.uniform(0, 2 * np.pi)
        a_peak = rng.uniform(0.4, 1.0) * P.A_MAX_TGT * 0.97

        g = _smooth((t - t_on) / 0.5)
        ax = g * np.sum(amp[0][:, None] * np.sin(2 * np.pi * f[0][:, None] * t + phi[0][:, None]), 0)
        ay = g * np.sum(amp[1][:, None] * np.sin(2 * np.pi * f[1][:, None] * t + phi[1][:, None]), 0)
        if has_break:
            w = _window(t, tb, tb + tb_len, tb_ramp)
            ax = ax + 1.5 * w * np.cos(bdir)
            ay = ay + 1.5 * w * np.sin(bdir)
        scale = a_peak / max(np.max(np.hypot(ax, ay)), 1e-9)
        ax_s, ay_s = ax * scale, ay * scale

        blackouts = []
        t1 = rng.uniform(P.T_BOO_MIN, 8.5)
        blackouts.append((t1, rng.uniform(0.1, P.T_BO_MAX)))
        if rng.random() < 0.3:
            t2 = rng.uniform(blackouts[0][0] + blackouts[0][1] + 0.5, 9.5)
            if t2 < 9.5:
                blackouts.append((t2, rng.uniform(0.1, P.T_BO_MAX)))

        d = dict(p0=p0, v0=v0, blackouts=blackouts,
                 accel_grid=(ax_s, ay_s))
        traj = _integrate(d, t)
        if _envelope_ok(traj, elev_margin=True):
            d['traj'] = traj
            return d
    raise RuntimeError(f"MC seed {seed}: no admissible trajectory found")


def _integrate(d, t):
    if 'accel_grid' in d:
        ax, ay = d['accel_grid']
    else:
        ax, ay = d['accel'](t)
        ax = np.broadcast_to(ax, t.shape).astype(float)
        ay = np.broadcast_to(ay, t.shape).astype(float)
    dt = t[1] - t[0]

    def cumtrapz(f):
        out = np.empty_like(f)
        out[0] = 0.0
        out[1:] = np.cumsum(0.5 * (f[1:] + f[:-1]) * dt)
        return out

    vx = d['v0'][0] + cumtrapz(ax)
    vy = d['v0'][1] + cumtrapz(ay)
    x = d['p0'][0] + cumtrapz(vx)
    y = d['p0'][1] + cumtrapz(vy)
    return dict(t=t, x=x, y=y, vx=vx, vy=vy, ax=ax, ay=ay)


def envelope_report(traj, t_end=P.T_END):
    m = traj['t'] <= t_end + 1e-9
    r = np.hypot(traj['x'][m], traj['y'][m])
    th = np.arctan2(traj['y'][m], traj['x'][m])
    return dict(r_min=r.min(), r_max=r.max(),
                v_max=np.hypot(traj['vx'][m], traj['vy'][m]).max(),
                a_max=np.hypot(traj['ax'][m], traj['ay'][m]).max(),
                th_min_deg=np.degrees(th.min()), th_max_deg=np.degrees(th.max()))


def _envelope_ok(traj, elev_margin=False):
    e = envelope_report(traj)
    ok = (e['r_min'] >= P.R_MIN and e['r_max'] <= P.R_MAX
          and e['v_max'] <= P.V_MAX_TGT and e['a_max'] <= P.A_MAX_TGT)
    if elev_margin:
        ok = ok and e['th_min_deg'] >= -8.0 and e['th_max_deg'] <= 83.0
    return ok


def _get(scenario):
    """Cached trajectory + noise. 'NAME#k' reuses NAME's trajectory with noise realisation k."""
    if scenario in _CACHE:
        return _CACHE[scenario]
    base, _, rep = scenario.partition('#')
    d = _scenario_def(base)
    t = np.arange(0.0, T_GRID_END + DT_GRID / 2, DT_GRID)
    traj = d['traj'] if 'traj' in d else _integrate(d, t)
    if not _envelope_ok(traj):
        raise ValueError(f"scenario {scenario} violates the operating envelope: {envelope_report(traj)}")
    n_ticks = int(round(T_GRID_END / P.DT_RADAR)) + 1
    if base.startswith('MC'):
        seed = int(base[2:]) + 10007
    elif base.startswith('HG'):
        seed = int(base[2:]) + 50021
    elif base.startswith('HX'):
        seed = int(base.rpartition('_')[2]) + 50021      # same noise as HG<seed>
    else:
        seed = sum(ord(c) * 131**k for k, c in enumerate(base)) % (2**32)
    if rep:
        seed = (seed + 7919 * (int(rep) + 1)) % (2**32)
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal((n_ticks, 2)) * np.array([P.SIGMA_R, P.SIGMA_TH])
    entry = dict(traj=traj, blackouts=d['blackouts'], noise=noise)
    _CACHE[scenario] = entry
    return entry


def target_truth(t, scenario):
    """True target state (x, y, vx, vy, ax, ay) at time t, linearly interpolated on the 1 ms grid."""
    e = _get(scenario)['traj']
    k = min(int(t / DT_GRID), len(e['t']) - 2)
    w = (t - e['t'][k]) / DT_GRID
    return tuple((1 - w) * e[n][k] + w * e[n][k + 1] for n in ('x', 'y', 'vx', 'vy', 'ax', 'ay'))


def blackout_intervals(scenario):
    return list(_get(scenario)['blackouts'])


def target_motion_generator(t, scenario):
    """Return (r_meas, theta_meas, meas_valid). A new measurement exists only at radar ticks
    outside blackout intervals; otherwise r_meas = theta_meas = NaN and meas_valid = 0."""
    e = _get(scenario)
    k = int(round(t / P.DT_RADAR))
    if abs(t - k * P.DT_RADAR) > 1e-6:
        return float('nan'), float('nan'), 0
    tk = k * P.DT_RADAR
    for (t_on, dur) in e['blackouts']:
        if t_on - 1e-9 <= tk < t_on + dur - 1e-9:
            return float('nan'), float('nan'), 0
    x, y = target_truth(tk, scenario)[:2]
    r = np.hypot(x, y) + e['noise'][k, 0]
    th = np.arctan2(y, x) + e['noise'][k, 1]
    return float(r), float(th), 1
