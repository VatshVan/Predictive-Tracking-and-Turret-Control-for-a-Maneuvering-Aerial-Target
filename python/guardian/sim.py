"""Closed-loop engagement simulation mirroring the evaluation call pattern.

Each controller step k at t_k = k / f_ctrl:
  [r, th, valid] = target_motion_generator(t_k, scenario)
  [V, fire]      = turret_controller(t_k, r, th, valid, theta, theta_dot)
  pulse judged at t_k if fire rises 0 -> 1 and >= 1 s since the last pulse
  plant integrated from t_k to t_{k+1} with V held (RK4 or forward Euler, step dt_plant),
  mechanical elevation stops applied to the state.
"""
import math

import numpy as np

from . import params as P
from . import turret_controller as tc
from .geometry import beam_miss_distance
from .target_motion_generator import target_motion_generator, target_truth
from .turret_plant import turret_plant


def _rk4(t, s, V, h):
    def f(tt, ss):
        return np.array(turret_plant(tt, V, ss[0], ss[1], ss[2]))
    k1 = f(t, s)
    k2 = f(t + h / 2, s + h / 2 * k1)
    k3 = f(t + h / 2, s + h / 2 * k2)
    k4 = f(t + h, s + h * k3)
    return s + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def _euler(t, s, V, h):
    return s + h * np.array(turret_plant(t, V, s[0], s[1], s[2]))


def _stops(s):
    th, thd, i = s
    if th <= P.TH_MIN and thd < 0:
        return np.array([P.TH_MIN, 0.0, i])
    if th >= P.TH_MAX and thd > 0:
        return np.array([P.TH_MAX, 0.0, i])
    return s


def run(scenario, f_ctrl=1000.0, dt_plant=1e-3, integrator='rk4', theta0=0.0,
        t_end=P.T_END, controller=None, log=True):
    """Simulate one engagement. Returns a dict with time histories and metrics."""
    if controller is None:
        tc.reset()
        controller = tc.turret_controller
    dt_c = 1.0 / f_ctrl
    n_steps = int(round(t_end / dt_c))
    n_sub = max(1, int(round(dt_c / dt_plant)))
    h = dt_c / n_sub
    step = _rk4 if integrator == 'rk4' else _euler

    s = np.array([theta0, 0.0, 0.0])
    fire_prev, t_pulse = 0, -1e9
    pulses = []
    L = {k: [] for k in ('t', 'th', 'thd', 'i', 'V', 'th_true', 'r_true', 'x', 'y',
                         'th_hat', 'r_ref', 'y_d', 'u_ff', 'u_fb', 'p_hit', 'cons', 'sig',
                         'mode', 'valid')}
    for k in range(n_steps + 1):
        t = k * dt_c
        r_m, th_m, valid = target_motion_generator(t, scenario)
        V, fire = controller(t, r_m, th_m, valid, float(s[0]), float(s[1]))
        x, y = target_truth(t, scenario)[:2]
        if fire == 1 and fire_prev == 0 and t - t_pulse >= P.T_RECHARGE - 1e-9:
            miss = beam_miss_distance(s[0], x, y)
            pulses.append(dict(t=t, miss=miss, hit=miss <= P.R_ENG))
            t_pulse = t
        fire_prev = fire
        if log:
            dbg = tc.get_debug() or {}
            L['t'].append(t); L['th'].append(s[0]); L['thd'].append(s[1]); L['i'].append(s[2])
            L['V'].append(V); L['x'].append(x); L['y'].append(y)
            L['th_true'].append(math.atan2(y, x)); L['r_true'].append(math.hypot(x, y))
            L['valid'].append(valid)
            for key in ('th_hat', 'r_ref', 'y_d', 'u_ff', 'u_fb', 'p_hit', 'cons', 'sig', 'mode'):
                L[key].append(dbg.get(key, float('nan')))
        if k == n_steps:
            break
        for j in range(n_sub):
            s = _stops(step(t + j * h, s, V, h))
            if not np.all(np.isfinite(s)) or abs(s[2]) > 1e6:
                raise FloatingPointError(f"plant diverged at t={t + j * h:.4f}")
    out = {k: np.array(v) for k, v in L.items()} if log else {}
    out['pulses'] = pulses
    out['metrics'] = metrics(out, pulses) if log else _pulse_metrics(pulses)
    return out


def _pulse_metrics(pulses):
    hits = [p for p in pulses if p['hit']]
    return dict(n_pulses=len(pulses), n_hits=len(hits),
                success=len(hits) > 0,
                t_first_hit=hits[0]['t'] if hits else float('nan'),
                hit_rate=len(hits) / len(pulses) if pulses else float('nan'),
                miss_median=float(np.median([p['miss'] for p in pulses])) if pulses else float('nan'))


def metrics(L, pulses):
    m = _pulse_metrics(pulses)
    err = L['r_true'] * np.sin(L['th'] - L['th_true'])          # cross-range pointing error [m]
    after = L['t'] >= 2.0
    m['point_rms_m'] = float(np.sqrt(np.mean(err[after] ** 2)))
    m['point_p95_m'] = float(np.percentile(np.abs(err[after]), 95))
    m['frac_in_Reng'] = float(np.mean(np.abs(err[after]) <= P.R_ENG))
    est = L['r_true'] * np.sin(L['th_hat'] - L['th_true'])
    ok = after & np.isfinite(est)
    m['est_rms_m'] = float(np.sqrt(np.mean(est[ok] ** 2)))
    m['V_peak'] = float(np.max(np.abs(L['V'])))
    return m
