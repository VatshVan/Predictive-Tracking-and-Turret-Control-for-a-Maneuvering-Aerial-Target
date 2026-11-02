"""Elevation-drive plant, derived from first principles (assignment Section 2).

Armature:   L di/dt   = V - R i - N K_e theta_dot
Turret:     J_eq th'' = N K_t i - b theta_dot,     J_eq = J_turret + N^2 J_m

theta(s)/V(s) = N K_t / ( s [ (L s + R)(J_eq s + b) + N^2 K_t K_e ] )
              = 1000 / ( s (320 s^2 + 3350 s + 1.0015e6) )

The supply voltage is clamped to +-600 V (docs/power_stage.md). The mechanical
elevation stops are applied by the integrator in sim.py because they act on the
state, not on its derivative.
"""
J_EQ = 1200.0 + 100.0**2 * 0.20
B = 1500.0
N = 100.0
K_T = 10.0
K_E = 10.0
R = 1.0
L = 0.10
V_MAX = 600.0


def turret_plant(t, V_m, theta_turret, theta_dot_turret, i):
    """Return [theta_dot, theta_ddot, i_dot] for the given motor voltage and state."""
    V = min(max(V_m, -V_MAX), V_MAX)
    i_dot = (V - R * i - N * K_E * theta_dot_turret) / L
    theta_ddot = (N * K_T * i - B * theta_dot_turret) / J_EQ
    return theta_dot_turret, theta_ddot, i_dot
