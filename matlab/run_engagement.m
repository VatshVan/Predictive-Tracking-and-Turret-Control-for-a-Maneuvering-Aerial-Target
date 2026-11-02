function out = run_engagement(scenario, f_ctrl, dt_plant, theta0)
%RUN_ENGAGEMENT Closed-loop simulation of one engagement.
%
%   out = run_engagement(scenario)                 % 100 Hz controller, 1 ms RK4 plant
%   out = run_engagement(scenario, f_ctrl, dt_plant, theta0)
%
%   Each controller step k at t_k = k / f_ctrl:
%     [r, th, valid] = target_motion_generator(t_k, scenario)
%     [V, fire]      = turret_controller(t_k, r, th, valid, theta, theta_dot)
%     a pulse is judged at t_k when fire rises 0 -> 1 and >= 1 s has passed since the
%     last pulse: hit if the target lies within 2 m of the beam ray at the turret angle
%     the plant is integrated to t_{k+1} with V held (RK4, step dt_plant <= 1 ms),
%     with mechanical stops at -10 and +85 deg.

if nargin < 2 || isempty(f_ctrl), f_ctrl = 100; end
if nargin < 3 || isempty(dt_plant), dt_plant = 1e-3; end
if nargin < 4 || isempty(theta0), theta0 = 0; end

T_END = 10; R_ENG = 2; T_RECHARGE = 1;
TH_MIN = deg2rad(-10); TH_MAX = deg2rad(85);
dt_c = 1 / f_ctrl;
n_steps = round(T_END / dt_c);
n_sub = max(1, round(dt_c / dt_plant));
h = dt_c / n_sub;

clear turret_controller
s = [theta0; 0; 0];
fire_prev = 0; t_pulse = -1e9;
pulses = zeros(0, 3);                       % [t, miss, hit]
L = nan(n_steps + 1, 12);
for k = 0:n_steps
    t = k * dt_c;
    [r_m, th_m, valid, truth] = target_motion_generator(t, scenario);
    [V, fire, dbg] = turret_controller(t, r_m, th_m, valid, s(1), s(2));
    x = truth(1); y = truth(2);
    if fire == 1 && fire_prev == 0 && t - t_pulse >= T_RECHARGE - 1e-9
        miss = beam_miss(s(1), x, y);
        pulses(end+1, :) = [t, miss, miss <= R_ENG]; %#ok<AGROW>
        t_pulse = t;
    end
    fire_prev = fire;
    L(k+1, :) = [t, s(1), s(2), s(3), V, atan2(y, x), hypot(x, y), dbg.th_hat, dbg.r_ref, ...
                 dbg.y_d, dbg.p_hit, valid];
    if k == n_steps, break, end
    for j = 0:n_sub-1
        s = rk4(t + j*h, s, V, h);
        if s(1) <= TH_MIN && s(2) < 0, s(1) = TH_MIN; s(2) = 0; end
        if s(1) >= TH_MAX && s(2) > 0, s(1) = TH_MAX; s(2) = 0; end
    end
end
out.t = L(:,1); out.th = L(:,2); out.thd = L(:,3); out.i = L(:,4); out.V = L(:,5);
out.th_true = L(:,6); out.r_true = L(:,7); out.th_hat = L(:,8); out.r_ref = L(:,9);
out.y_d = L(:,10); out.p_hit = L(:,11); out.valid = L(:,12);
out.pulses = pulses;
hits = pulses(pulses(:,3) == 1, :);
err = out.r_true .* sin(out.th - out.th_true);
after = out.t >= 2;
out.metrics = struct( ...
    'n_pulses', size(pulses, 1), 'n_hits', size(hits, 1), 'success', ~isempty(hits), ...
    't_first_hit', first_or_nan(hits), ...
    'point_rms_m', sqrt(mean(err(after).^2)), ...
    'frac_in_Reng', mean(abs(err(after)) <= R_ENG), ...
    'V_peak', max(abs(out.V)));
end

function v = first_or_nan(hits)
if isempty(hits), v = NaN; else, v = hits(1, 1); end
end

function d = beam_miss(theta_beam, x, y)
u = [cos(theta_beam), sin(theta_beam)];
if x*u(1) + y*u(2) < 0
    d = hypot(x, y);
else
    d = abs(x*u(2) - y*u(1));
end
end

function s = rk4(t, s, V, h)
k1 = f(t, s, V);
k2 = f(t + h/2, s + h/2*k1, V);
k3 = f(t + h/2, s + h/2*k2, V);
k4 = f(t + h, s + h*k3, V);
s = s + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

function ds = f(t, s, V)
[a, b, c] = turret_plant(t, V, s(1), s(2), s(3));
ds = [a; b; c];
end
