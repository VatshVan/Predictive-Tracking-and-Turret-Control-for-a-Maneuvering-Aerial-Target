function [r_meas, theta_meas, meas_valid, truth, info] = target_motion_generator(t, scenario)
%TARGET_MOTION_GENERATOR Target trajectory and 20 Hz radar measurements.
%
%   [r_meas, theta_meas, meas_valid] = target_motion_generator(t, scenario)
%   [~, ~, ~, truth, info] = ...   also returns the true state at t
%                                  (truth = [x y vx vy ax ay]) and info.blackouts.
%
%   scenario: 'A' (or 1)  nominal non-maneuvering target
%             'B' (or 2)  smooth bounded maneuver in x and y from t = 1 s
%             'C' (or 3)  nominal target, radar loss on [1.0, 1.5) s
%             'D1' stealth break: maneuver begins inside a late blackout
%             'D2' weave: near-a_max sinusoidal jinking in x and y
%             'D3' close-in pass: r -> 850 m, elevation approaching +85 deg
%             'D4' repeated blackouts on the Scenario B maneuver
%             'D5' combined: weave + late blackout + close range
%             'D6' jerk test: near bang-bang acceleration reversals
%
%   Trajectories are integrated once on a 1 ms grid, checked against the operating
%   envelope (800 <= r <= 2000 m, |v| <= 150 m/s, |a| <= 20 m/s^2) and cached.
%   Measurements exist only at radar ticks t = k*0.05 s outside blackouts; otherwise
%   r_meas = theta_meas = NaN and meas_valid = 0. Noise is indexed by tick number and
%   drawn from a private RandStream, so the global RNG is untouched.
%   Python reference: python/guardian/target_motion_generator.py

persistent cache
if isempty(cache)
    cache = containers.Map();
end
if isnumeric(scenario)
    names = {'A', 'B', 'C'};
    scenario = names{scenario};
end
if ~isKey(cache, scenario)
    cache(scenario) = build(scenario);
end
e = cache(scenario);

DT_RADAR = 0.05;
k = round(t / DT_RADAR);
r_meas = NaN; theta_meas = NaN; meas_valid = 0;
if abs(t - k * DT_RADAR) <= 1e-6
    tk = k * DT_RADAR;
    blocked = false;
    for b = 1:size(e.blackouts, 1)
        if tk >= e.blackouts(b,1) - 1e-9 && tk < e.blackouts(b,1) + e.blackouts(b,2) - 1e-9
            blocked = true;
        end
    end
    if ~blocked
        s = interp_state(e, tk);
        r_meas = hypot(s(1), s(2)) + e.noise(k+1, 1);
        theta_meas = atan2(s(2), s(1)) + e.noise(k+1, 2);
        meas_valid = 1;
    end
end
if nargout > 3
    truth = interp_state(e, t);
    info.blackouts = e.blackouts;
end
end

% ------------------------------------------------------------------------------
function s = interp_state(e, t)
dt = 1e-3;
k = min(floor(t / dt), numel(e.t) - 2) + 1;
w = (t - e.t(k)) / dt;
s = (1 - w) * e.X(k, :) + w * e.X(k+1, :);
end

function e = build(scenario)
dt = 1e-3;
t = (0:dt:11.0)';
z = zeros(size(t));
switch scenario
    case {'A', 'C'}
        p0 = [1700, 700]; v0 = [-120, 0]; ax = z; ay = z;
        if strcmp(scenario, 'C'), bo = [1.0, 0.5]; else, bo = zeros(0, 2); end
    case {'B', 'D4'}
        g = smooth5((t - 1.0) / 0.5);
        ax = g .* 8 .* sin(2*pi*0.15*(t - 1.0));
        ay = g .* 14 .* sin(2*pi*0.30*(t - 1.0));
        p0 = [1700, 700]; v0 = [-120, 0];
        if strcmp(scenario, 'D4')
            bo = [1.0 0.5; 3.2 0.4; 5.5 0.5; 7.5 0.4];
        else
            bo = zeros(0, 2);
        end
    case 'D1'
        w = smooth5((t - 3.1) / 0.3) - smooth5((t - 4.4) / 0.3);
        ax = 9 * w; ay = 17 * w;
        p0 = [1900, 500]; v0 = [-125, 0]; bo = [3.0, 0.5];
    case {'D2', 'D5'}
        g = smooth5((t - 0.5) / 0.5);
        ax = g .* 13 .* sin(2*pi*t/1.6);
        ay = g .* 14 .* cos(2*pi*t/1.2);
        if strcmp(scenario, 'D2')
            p0 = [1800, 800]; v0 = [-110, 0]; bo = zeros(0, 2);
        else
            p0 = [1300, 700]; v0 = [-80, 0]; bo = [4.0, 0.5];
        end
    case 'D3'
        g = smooth5((t - 1.0) / 0.5);
        ax = z; ay = g .* 4 .* sin(2*pi*0.4*(t - 1.0));
        p0 = [1450, 835]; v0 = [-135, 0]; bo = zeros(0, 2);
    case 'D6'
        sgn = smooth5((t - 1.0) / 0.15);
        edges = 1.0 + 0.6 + 1.2 * (0:7);
        for kk = 1:numel(edges)
            if mod(kk - 1, 2) == 0, c = -2; else, c = 2; end
            sgn = sgn + c * smooth5((t - edges(kk)) / 0.15);
        end
        ax = z; ay = 18 * sgn;
        p0 = [1750, 850]; v0 = [-125, 0]; bo = zeros(0, 2);
    otherwise
        error('unknown scenario %s', scenario);
end
vx = v0(1) + cumtrapz(t, ax); vy = v0(2) + cumtrapz(t, ay);
x = p0(1) + cumtrapz(t, vx);  y = p0(2) + cumtrapz(t, vy);
e.t = t;
e.X = [x, y, vx, vy, ax, ay];
e.blackouts = bo;
m = t <= 10 + 1e-9;
r = hypot(x(m), y(m));
assert(min(r) >= 800 && max(r) <= 2000, 'scenario %s violates the range envelope', scenario);
assert(max(hypot(vx(m), vy(m))) <= 150, 'scenario %s violates the speed limit', scenario);
assert(max(hypot(ax(m), ay(m))) <= 20, 'scenario %s violates the acceleration limit', scenario);
seed = mod(sum(double(scenario) .* 131 .^ (0:numel(scenario)-1)), 2^32);
rs = RandStream('mt19937ar', 'Seed', seed);
n_ticks = round(11.0 / 0.05) + 1;
e.noise = randn(rs, n_ticks, 2) .* [1.0, deg2rad(0.02)];
end

function s = smooth5(tau)
% Quintic smoothstep (C2): 0 for tau <= 0, 1 for tau >= 1.
s = min(max(tau, 0), 1);
s = s.^3 .* (10 - 15*s + 6*s.^2);
end
