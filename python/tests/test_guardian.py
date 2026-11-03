import math

import numpy as np
import pytest

from guardian import params as P
from guardian import turret_controller as tc
from guardian.geometry import beam_miss_distance
from guardian.sim import run
from guardian.target_motion_generator import (_get, envelope_report, target_motion_generator,
                                              target_truth)
from guardian.turret_plant import turret_plant

SCENARIOS = ['A', 'B', 'C', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6']


def test_plant_matches_derived_transfer_function():
    # Linearise turret_plant numerically and compare its poles with
    # theta/V = 1000 / (s (320 s^2 + 3350 s + 1.0015e6)).
    eps = 1e-6
    A = np.zeros((3, 3))
    base = np.array(turret_plant(0.0, 0.0, 0.0, 0.0, 0.0))
    for k in range(3):
        s = np.zeros(3)
        s[k] = eps
        A[:, k] = (np.array(turret_plant(0.0, 0.0, s[0], s[1], s[2])) - base) / eps
    poles = np.sort_complex(np.linalg.eigvals(A))
    expected = np.sort_complex(np.concatenate([[0.0], np.roots([320, 3350, 1.0015e6])]))
    assert np.allclose(poles, expected, atol=1e-3)


def test_plant_voltage_saturates():
    assert turret_plant(0, 1e4, 0, 0, 0)[2] == pytest.approx(600.0 / 0.1)


@pytest.mark.parametrize('scenario', SCENARIOS + [f'MC{i}' for i in range(20)])
def test_scenarios_respect_envelope(scenario):
    e = envelope_report(_get(scenario)['traj'])
    assert P.R_MIN <= e['r_min'] and e['r_max'] <= P.R_MAX
    assert e['v_max'] <= P.V_MAX_TGT and e['a_max'] <= P.A_MAX_TGT


def test_blackout_and_radar_ticks():
    assert target_motion_generator(1.0, 'C')[2] == 0
    assert target_motion_generator(1.45, 'C')[2] == 0
    assert target_motion_generator(1.5, 'C')[2] == 1
    assert target_motion_generator(0.95, 'C')[2] == 1
    assert target_motion_generator(0.951, 'A')[2] == 0      # between ticks


def test_beam_miss_distance():
    assert beam_miss_distance(0.0, 100.0, 1.5) == pytest.approx(1.5)
    assert beam_miss_distance(math.pi / 4, 100.0, 100.0) == pytest.approx(0.0, abs=1e-9)
    assert beam_miss_distance(0.0, -10.0, 0.0) == pytest.approx(10.0)  # behind the turret


def test_los_derivatives_match_finite_differences():
    x, y = 1200.0, 600.0
    vx, vy, ax, ay, jx, jy = -110.0, 15.0, 8.0, -12.0, 30.0, 20.0

    def theta(s):
        px = x + vx * s + ax * s * s / 2 + jx * s ** 3 / 6
        py = y + vy * s + ay * s * s / 2 + jy * s ** 3 / 6
        return math.atan2(py, px)

    st = np.array([x, vx, ax, jx, y, vy, ay, jy])
    _, th0, th1, th2, th3 = tc._los(st, 4)
    h = 1e-3
    fd1 = (theta(h) - theta(-h)) / (2 * h)
    fd2 = (theta(h) - 2 * theta(0) + theta(-h)) / h ** 2
    fd3 = (theta(2 * h) - 2 * theta(h) + 2 * theta(-h) - theta(-2 * h)) / (2 * h ** 3)
    assert th0 == pytest.approx(theta(0))
    assert th1 == pytest.approx(fd1, rel=1e-6)
    assert th2 == pytest.approx(fd2, rel=1e-4)
    assert th3 == pytest.approx(fd3, rel=1e-3)


def test_taylor_sections_hurwitz_up_to_degree_4():
    for n in range(1, 7):
        c = [P.TP ** k / math.factorial(k) for k in range(n, -1, -1)]
        stable = bool(np.all(np.roots(c).real < 0))
        assert stable == (n <= 4)


def test_controller_resets_between_runs():
    tc.reset()
    a = run('A', f_ctrl=50.0)['metrics']
    b = run('A', f_ctrl=50.0)['metrics']
    assert a == b


@pytest.mark.parametrize('f_ctrl', [20.0, 50.0, 1000.0])
def test_standard_scenarios_intercept(f_ctrl):
    for s in ['A', 'B', 'C']:
        m = run(s, f_ctrl=f_ctrl)['metrics']
        assert m['success'], (s, f_ctrl, m)
        assert m['t_first_hit'] < 2.5
        assert m['V_peak'] <= 600.0 + 1e-9


def test_euler_plant_integration_at_1ms():
    m = run('B', f_ctrl=100.0, integrator='euler')['metrics']
    assert m['success']


def test_los_fourth_derivative_matches_finite_differences():
    x, y = 1200.0, 600.0
    d = dict(vx=-110.0, vy=15.0, ax=8.0, ay=-12.0, jx=30.0, jy=20.0, sx=-40.0, sy=25.0)

    def theta(s):
        px = x + d['vx'] * s + d['ax'] * s ** 2 / 2 + d['jx'] * s ** 3 / 6 + d['sx'] * s ** 4 / 24
        py = y + d['vy'] * s + d['ay'] * s ** 2 / 2 + d['jy'] * s ** 3 / 6 + d['sy'] * s ** 4 / 24
        return math.atan2(py, px)

    st = np.array([x, d['vx'], d['ax'], d['jx'], d['sx'], y, d['vy'], d['ay'], d['jy'], d['sy']])
    th4 = tc._los(st, 5, nd=4)[5]
    h = 1e-2
    fd4 = (theta(2 * h) - 4 * theta(h) + 6 * theta(0) - 4 * theta(-h) + theta(-2 * h)) / h ** 4
    assert th4 == pytest.approx(fd4, rel=1e-3)


@pytest.mark.parametrize('kw, hurwitz', [
    (dict(ref_mode='exact_taylor', ref_n=3), True),
    (dict(ref_mode='exact_taylor', ref_n=4, order=5), True),
    (dict(ref_mode='exact_taylor', ref_n=5, order=6), False),
    (dict(ref_mode='moment_matched'), True),
    (dict(ref_mode='contact'), True),
    (dict(ref_mode='contact', ref_poly=(1.0, 0.5, 0.125, 1 / 64, 1 / 1024)), True),
])
def test_exact_reference_models(kw, hurwitz):
    saved = dict(tc.CFG)
    try:
        tc.CFG.update(kw)
        S = tc._init(0.0, 0.0)
        n = len(S['A_m']) - 1
        if kw['ref_mode'] == 'contact':
            # contact order k = number of leading coefficients shared with A_3, minus one
            k = len(kw.get('ref_poly', tc.DEFAULT_CFG['ref_poly'])) - 3
            n = k
        # T A_n = 1 + O(s^(contact order + 1)): Q agrees with A_n up to the contact order
        assert S['D'][:n + 1] == pytest.approx(S['A_m'][:n + 1])
        assert bool(np.all(np.roots(S['D'][::-1]).real < 0)) == hurwitz
    finally:
        tc.CFG.clear()
        tc.CFG.update(saved)


def test_homotopy_signature_and_classes():
    from guardian.homotopy import class_signs, classify, signature
    assert list(signature([3.0, 1.0, -0.5, -2.5])) == [1, 0, 0, -1]
    for cls in range(5):
        assert classify(class_signs(cls, 6)) == cls      # each intended class is distinct
    assert classify([0, 0, 0]) == -1                      # every event passed through: all hits
    assert classify([1, 0, -1]) == 2                      # 'through' events do not change the class
    assert list(tc._h_signature([2.5, -2.5, 0.0])) == [tc.H_ABOVE, tc.H_BELOW, tc.H_THROUGH]


def test_homotopy_planner_probabilities():
    saved = dict(tc.CFG)
    try:
        tc.CFG.update(phit_mode='homotopy')
        tc.reset()
        run('A', f_ctrl=100.0, t_end=3.0)
        S = tc._S
        taus, P_thr, sgn, w, sig, rng = tc._h_plan(S, S['t_prev'], math.atan2(S['x'][S['n']], S['x'][0]), 0.5)
        assert taus[0] == 0.0 and len(taus) == 11
        assert np.all((P_thr >= 0.0) & (P_thr <= 1.0))
        assert w.sum() == pytest.approx(1.0) and sgn.shape == (len(w), len(taus))
    finally:
        tc.CFG.clear()
        tc.CFG.update(saved)


def test_blend_interpolates_reference_polynomials():
    saved = dict(tc.CFG)
    try:
        polys = {}
        for a in (0.0, 0.5, 1.0):
            tc.CFG.update(ref_mode='blend', blend_a=a)
            polys[a] = np.array(tc._init(0.0, 0.0)['D'])
        A3 = np.array([P.TP ** k / math.factorial(k) for k in range(4)])
        assert polys[0.0] == pytest.approx(A3)
        assert polys[0.5] == pytest.approx(0.5 * (polys[0.0] + polys[1.0]))
        assert np.roots(polys[1.0][::-1]).real == pytest.approx([-20.0] * 3, abs=1e-3)
    finally:
        tc.CFG.clear()
        tc.CFG.update(saved)


@pytest.mark.parametrize('kw', [dict(ref_mode='exact_taylor', ref_n=3), dict(ref_mode='moment_matched'),
                                dict(ref_mode='blend', blend_a=0.3), dict(phit_mode='homotopy')])
def test_research_modes_intercept(kw):
    saved = dict(tc.CFG)
    try:
        tc.CFG.update(kw)
        for s in ['A', 'C']:
            m = run(s, f_ctrl=100.0)['metrics']
            assert m['success'], (kw, s, m)
            assert m['V_peak'] <= 600.0 + 1e-9
    finally:
        tc.CFG.clear()
        tc.CFG.update(saved)
