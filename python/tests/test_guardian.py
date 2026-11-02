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
