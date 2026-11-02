function [theta_dot_turret, theta_ddot_turret, i_dot] = turret_plant(t, V_m, theta_turret, theta_dot_turret, i) %#ok<INUSL>
%TURRET_PLANT Elevation-drive dynamics derived from first principles.
%
%   Armature:  L di/dt     = V - R i - N K_e theta_dot
%   Turret:    J_eq th''   = N K_t i - b theta_dot,   J_eq = J_turret + N^2 J_m = 3200 kg m^2
%
%   theta(s)/V(s) = N K_t / ( s [ (L s + R)(J_eq s + b) + N^2 K_t K_e ] )
%                 = 1000 / ( s (320 s^2 + 3350 s + 1.0015e6) )
%
%   Poles: 0, -5.23 +- 55.70j (8.90 Hz, zeta = 0.094). Integrate with a step <= 1 ms
%   (forward Euler is unstable above 3.34 ms). The supply is limited to +-600 V
%   (docs/power_stage.md). Mechanical stops (-10..+85 deg) act on the state and are
%   applied by the integrator (see run_engagement.m).

J_eq = 1200 + 100^2 * 0.20;
b = 1500; N = 100; K_t = 10; K_e = 10; R = 1.0; L = 0.10;
V_MAX = 600;

V = min(max(V_m, -V_MAX), V_MAX);
i_dot = (V - R * i - N * K_e * theta_dot_turret) / L;
theta_ddot_turret = (N * K_t * i - b * theta_dot_turret) / J_eq;
end
