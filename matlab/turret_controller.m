function [V_m, fire_cmd, dbg] = turret_controller(t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret)
%TURRET_CONTROLLER Guardian predictive-feedback controller for the turret elevation axis.
%
%   [V_m, fire_cmd] = turret_controller(t, r_meas, theta_meas, meas_valid, ...
%                                       theta_turret, theta_dot_turret)
%
%   Pipeline:
%     1. Target-state estimation: converted-measurement Kalman filter (unbiased
%        polar->Cartesian conversion), Markov-jerk motion model (position..jerk per
%        axis), envelope projection, track-jump detection.
%     2. Target prediction: exact LOS-angle derivatives (complex log-derivatives) and
%        the Taylor advance A3(s) = 1 + tp s + (tp s)^2/2 + (tp s)^3/6, tp = 0.5 s.
%     3. Predicted elevation reference r(t) = A3(d/dt) theta_hat(t).
%     4. 2-DOF LTI control: reference generator (A3 + K)(D) y_d = r + K(D) theta_hat
%        (the A3 lag cancels the tp lead; K places the transient poles), exact
%        plant-inverse feed-forward, PI + roll-off feedback C(s).
%     5. Fire control: Gaussian-mixture hit probability over envelope-bounded
%        maneuver classes, class consensus, fire-now-or-wait rule, 1 s recharge.
%     6. Measurement loss: coasting and maneuver-class re-identification on recovery.
%
%   State is kept in a persistent struct and re-initialised automatically when t
%   decreases or t == 0 (or after `clear turret_controller`).
%   The optional third output returns internal signals for plotting.
%
%   Python reference implementation: python/guardian/turret_controller.py

persistent S
if isempty(S) || t < S.t_prev - 1e-12 || t == 0
    S = init_state(theta_turret, theta_dot_turret);
end
C = S.cfg;
if S.first
    dt = 0;
else
    dt = t - S.t_prev;
end
S.first = false;
n = S.n;

% 1. advance the reference generator over the last interval (sub-steps <= 1 ms)
if dt > 0 && ~isempty(S.derivs_prev)
    n_sub = max(1, ceil(dt / 1e-3 - 1e-9));
    h = dt / n_sub;
    m = S.m_ref;
    for k = 0:n_sub-1
        y3 = ref_jerk(S, m, taylor_adv(S.derivs_prev, k * h));
        m(3) = m(3) + h * y3;
        m(2) = m(2) + h * m(3);
        m(1) = m(1) + h * m(2);
    end
    S.m_ref = m;
end

% 2. estimation
S = estimator_step(S, t, r_meas, theta_meas, meas_valid);

% 3. nowcast + tp-prediction
mode = 'acquire';
y = S.m_ref(1); y1 = S.m_ref(2); y2 = S.m_ref(3);
derivs = [y, 0, 0, 0];
p_hit = 0; cons = 0; sig = NaN; th_hat = NaN; rng_hat = NaN;
if S.trk
    tau = max(t - S.t_upd, 0);
    [F, Q] = disc_target(S, tau);
    x_now = F * S.x;
    P_now = F * S.P * F' + Q;
    [rng_hat, th_hat, th1, th2, th3] = los(x_now, n);
    derivs = [th_hat, th1, th2, th3];
    mode = 'track';
    if tau > 1.5 * S.DT_RADAR
        mode = 'coast';
        if C.use_ambush && tau > 0.25
            if isempty(S.ambush)
                S.ambush = ambush_angle(S, S.x, tau, th1);
            end
            if ~isempty(S.ambush)
                derivs = [S.ambush, 0, 0, 0];
                mode = 'ambush';
            end
        end
    else
        S.ambush = [];
    end
    [p_hit, cons, sig] = hit_probability(S, x_now, P_now, tau, theta_turret);
end
lo = S.TH_MIN + S.TH_MARGIN;
hi = S.TH_MAX - S.TH_MARGIN;
if ~(derivs(1) >= lo && derivs(1) <= hi)
    derivs = [min(max(derivs(1), lo), hi), 0, 0, 0];
end
% predicted elevation reference: Taylor advance of the nowcast by tp
r_pred = sum(S.A_m .* derivs(1:numel(S.A_m)));
S.derivs_prev = derivs;

% 4-5. reference generator + 2-DOF control
y3 = ref_jerk(S, [y, y1, y2], derivs);
u_ff = S.FF(1) * y3 + S.FF(2) * y2 + S.FF(3) * y1;
e = y - theta_turret;
if dt > 0
    if ~(S.sat ~= 0 && S.sat * e > 0)
        S.xi = S.xi + e * dt;                 % conditional integration (anti-windup)
    end
    a = exp(-dt / C.tau_f);
    S.xf = a * S.xf + (1 - a) * C.kp * (e + S.xi / C.ti);
end
u = u_ff + S.xf;
if u > S.V_LIM
    V_m = S.V_LIM; S.sat = 1;
elseif u < -S.V_LIM
    V_m = -S.V_LIM; S.sat = -1;
else
    V_m = u; S.sat = 0;
end

% 6. fire control: fire-now-or-wait over the next recharge window
fire_cmd = 0;
ready = (t - S.t_fire >= S.T_RECHARGE + 1e-6) && S.trk && S.n_upd >= C.n_upd_min;
if ready && cons >= C.q_cons
    if strcmp(mode, 'coast') || strcmp(mode, 'ambush')
        if p_hit >= C.p_min && p_hit < S.p_prev
            fire_cmd = 1;                       % certainty can only decay: shoot now
        end
    elseif p_hit >= C.p_sure
        fire_cmd = 1;
    elseif p_hit >= C.p_fire
        if isempty(S.t_armed)
            S.t_armed = t; S.p_armed_max = p_hit;
        end
        S.p_armed_max = max(S.p_armed_max, p_hit);
        if p_hit < S.p_armed_max - 0.02 || t - S.t_armed >= C.t_wait
            fire_cmd = 1;                       % peak passed, or waited long enough
        end
    else
        S.t_armed = [];
    end
end
if fire_cmd
    S.t_fire = t;
    S.t_armed = [];
end
S.p_prev = p_hit;
S.t_prev = t;

if nargout > 2
    dbg = struct('mode', mode, 'th_hat', th_hat, 'r_hat', rng_hat, 'r_ref', r_pred, ...
                 'y_d', y, 'u_ff', u_ff, 'u_fb', S.xf, 'p_hit', p_hit, 'cons', cons, 'sig', sig);
end
end

% =============================================================================== setup
function cfg = default_cfg()
cfg = struct( ...
    'order', 4, ...          % states per axis: 2 = CV, 3 = Singer accel, 4 = Markov jerk
    'sigma_top', 20.0, ...   % std of the top-derivative Markov process [m/s^3]
    'tau_top', 0.4, ...      % its correlation time [s]
    'q_cv', 200.0, ...       % white-acceleration intensity for order 2 [m^2/s^3]
    'tau_fast', 0.02, ...    % extra reference pole for orders < 4 [s]
    'kp', 4000.0, ...        % feedback PI gain [V/rad]
    'ti', 1.0, ...           % integral time [s]
    'tau_f', 0.03, ...       % feedback roll-off [s]
    'p_dev', 20.0, ...       % deviation-dynamics triple pole [rad/s]
    'omega_ref_max', 0.58, ...  % acquisition rate limit [rad/s]
    'alpha_ref_max', 2.5, ...   % acquisition acceleration limit [rad/s^2]
    'use_reid', true, ...    % maneuver-class re-identification after measurement loss
    'use_ambush', false, ... % ambush hold during long measurement loss
    'p_sure', 0.99, ...      % tracking: fire immediately above this hit probability
    'p_fire', 0.90, ...      % tracking: arm; fire at the peak or after t_wait
    't_wait', 0.30, ...      % tracking: longest wait once armed [s]
    'p_min', 0.95, ...       % coasting: fire only while certainty is near tracking quality
    'q_cons', 0.80, ...      % minimum fraction of maneuver classes inside R_ENG
    'n_upd_min', 8, ...      % radar updates before the track may fire
    'chi2_jump', 40.0);      % innovation gate for track-jump detection
end

function S = init_state(theta, theta_dot)
cfg = default_cfg();
S.cfg = cfg;
S.TP = 0.5; S.R_ENG = 2.0; S.T_RECHARGE = 1.0; S.DT_RADAR = 0.05;
S.SIG_R = 1.0; S.SIG_TH = deg2rad(0.02);
S.V_TGT_MAX = 150.0; S.A_TGT_MAX = 20.0;
S.TH_MIN = deg2rad(-10); S.TH_MAX = deg2rad(85); S.TH_MARGIN = deg2rad(0.3);
S.FF = [0.32, 3.35, 1001.5];             % plant inverse: V = 0.32 y''' + 3.35 y'' + 1001.5 y'
S.V_LIM = 600.0;

n = cfg.order;
S.n = n;
A = zeros(n);
for k = 1:n-1
    A(k, k+1) = 1;
end
if n == 2
    q = cfg.q_cv;
else
    alpha = 1 / cfg.tau_top;
    A(n, n) = -alpha;
    q = 2 * alpha * cfg.sigma_top^2;
end
S.A = A; S.q = q;
m = n - 1;
S.A_m = S.TP .^ (0:m) ./ factorial(0:m);
D = S.A_m;
for k = 1:(3 - m)
    D = conv(D, [1, cfg.tau_fast]);       % ascending coefficients
end
S.D = D;
H = zeros(2, 2*n);
H(1, 1) = 1; H(2, n+1) = 1;
S.H = H;
S.classes = [0, 0];
for mag = [0.5, 1.0] * S.A_TGT_MAX
    for k = 0:7
        S.classes(end+1, :) = mag * [cos(k*pi/4), sin(k*pi/4)];
    end
end
S.t_prev = 0; S.first = true;
S.trk = false; S.x = []; S.P = []; S.t_upd = -1; S.n_upd = 0;
S.last_meas = [NaN, NaN]; S.n_jump = 0;
S.m_ref = [theta, theta_dot, 0]; S.derivs_prev = [];
S.xi = 0; S.xf = 0; S.sat = 0;
S.t_fire = -1e9; S.p_prev = 0; S.t_armed = []; S.p_armed_max = 0; S.ambush = [];
end

% =============================================================================== models
function [F, Q] = disc_target(S, dt)
% Van Loan discretisation of the 2-axis target model.
n = S.n;
Qc = zeros(n); Qc(n, n) = S.q;
M = [-S.A, Qc; zeros(n), S.A'];
E = expm(M * dt);
F1 = E(n+1:end, n+1:end)';
Q1 = F1 * E(1:n, n+1:end);
Q1 = 0.5 * (Q1 + Q1');
F = kron(eye(2), F1);
Q = kron(eye(2), Q1);
end

function [z, R] = converted(S, r, th)
% Unbiased converted measurement and covariance (Lerro & Bar-Shalom / Mo et al.).
lam = exp(-S.SIG_TH^2 / 2);
lam2 = exp(-2 * S.SIG_TH^2);
c = cos(th); s = sin(th);
z = [r * c; r * s] / lam;
k = (lam^-2 - 2) * r * r;
h = 0.5 * (r * r + S.SIG_R^2);
R11 = k * c * c + h * (1 + lam2 * cos(2 * th));
R22 = k * s * s + h * (1 - lam2 * cos(2 * th));
R12 = k * c * s + h * lam2 * sin(2 * th);
R = [R11, R12; R12, R22];
end

function [rng, th0, th1, th2, th3] = los(x, n)
% Range, LOS angle and its first three derivatives from the Cartesian state.
zs = zeros(1, 4);
for k = 0:3
    if k < n
        zs(k+1) = complex(x(k+1), x(n+k+1));
    end
end
z0 = zs(1);
w1 = zs(2) / z0; w2 = zs(3) / z0; w3 = zs(4) / z0;
th0 = atan2(imag(z0), real(z0));
th1 = imag(w1);
th2 = imag(w2 - w1 * w1);
th3 = imag(w3 - 3 * w1 * w2 + 2 * w1^3);
rng = abs(z0);
end

function x = project_envelope(S, x)
lims = [1.05 * S.V_TGT_MAX, 1.15 * S.A_TGT_MAX, 400];
for k = 1:3
    if k >= S.n
        break
    end
    mag = hypot(x(k+1), x(S.n+k+1));
    if mag > lims(k)
        x(k+1) = x(k+1) * lims(k) / mag;
        x(S.n+k+1) = x(S.n+k+1) * lims(k) / mag;
    end
end
end

% =============================================================================== estimator
function S = kf_update(S, z, R)
H = S.H; P = S.P;
nu = z - H * S.x;
Sm = H * P * H' + R;
K = P * H' / Sm;
S.x = S.x + K * nu;
I_KH = eye(size(P, 1)) - K * H;
S.P = I_KH * P * I_KH' + K * R * K';
S.x = project_envelope(S, S.x);
end

function S = track_init(S, z, R, t)
n = S.n;
x = zeros(2*n, 1);
x(1) = z(1); x(n+1) = z(2);
P = zeros(2*n);
P([1, n+1], [1, n+1]) = R;
sig0 = [NaN, S.V_TGT_MAX, S.A_TGT_MAX, S.cfg.sigma_top];
for k = 1:n-1
    P(k+1, k+1) = sig0(k+1)^2;
    P(n+k+1, n+k+1) = sig0(k+1)^2;
end
S.trk = true; S.x = x; S.P = P; S.t_upd = t; S.n_upd = 1; S.n_jump = 0;
end

function S = reidentify(S, z, R, F, Q)
% Maneuver-class re-identification after measurement loss (T-MPC-style consistent
% selection: the nominal hypothesis carries the larger prior weight).
n = S.n;
x_b = S.x;
Pp = F * S.P * F' + Q;
H = S.H;
Sinv = inv(H * Pp * H' + R);
n_cls = size(S.classes, 1);
best_cost = Inf; best_x = [];
for i = 1:n_cls
    xb = x_b;
    if i == 1
        prior = 0.5;
    else
        prior = 0.5 / (n_cls - 1);
        if n >= 3
            xb(3) = S.classes(i, 1); xb(n+3) = S.classes(i, 2);
            if n >= 4
                xb(4) = 0; xb(n+4) = 0;
            end
        else
            continue
        end
    end
    xi = F * xb;
    nu = z - H * xi;
    cost = 0.5 * (nu' * Sinv * nu) - log(prior);
    if cost < best_cost
        best_cost = cost; best_x = xi;
    end
end
S.x = best_x;
S.P = Pp;
end

function S = estimator_step(S, t, r_meas, th_meas, valid)
if ~valid || isnan(r_meas) || isnan(th_meas)
    return
end
if isequal([r_meas, th_meas], S.last_meas)
    return                                      % held duplicate, not a new measurement
end
S.last_meas = [r_meas, th_meas];
[z, R] = converted(S, r_meas, th_meas);
if ~S.trk
    S = track_init(S, z, R, t);
    return
end
gap = t - S.t_upd;
if gap <= 0
    return
end
[F, Q] = disc_target(S, gap);
x_pred = F * S.x;
P_pred = F * S.P * F' + Q;
nu = z - S.H * x_pred;
Sm = S.H * P_pred * S.H' + R;
d2 = nu' * (Sm \ nu);
if d2 > S.cfg.chi2_jump
    S.n_jump = S.n_jump + 1;
    if S.n_jump >= 2                            % persistent jump: new target, restart track
        S = track_init(S, z, R, t);
    end
    return
end
S.n_jump = 0;
if gap > 1.5 * S.DT_RADAR && S.cfg.use_reid
    S = reidentify(S, z, R, F, Q);
else
    S.x = x_pred; S.P = P_pred;
end
S = kf_update(S, z, R);
S.t_upd = t;
S.n_upd = S.n_upd + 1;
end

% =============================================================================== reference
function y3 = ref_jerk(S, m, derivs)
% Linear law: y''' = phi - 3p d2 - 3p^2 d1 - p^3 d0, d_j = y^(j) - theta_hat^(j),
% realised as a cascade (gains p/3, p, 3p) with saturated rate/acceleration demands.
D = S.D;
p = S.cfg.p_dev;
r_pred = sum(S.A_m .* derivs(1:numel(S.A_m)));
phi = (r_pred - D(1) * derivs(1) - D(2) * derivs(2) - D(3) * derivs(3)) / D(4);
om = S.cfg.omega_ref_max; al = S.cfg.alpha_ref_max;
y1_des = min(max(derivs(2) - (p / 3) * (m(1) - derivs(1)), -om), om);
y2_des = min(max(derivs(3) - p * (m(2) - y1_des), -al), al);
y3 = phi - 3 * p * (m(3) - y2_des);
end

function d = taylor_adv(derivs, s)
d = [derivs(1) + derivs(2) * s + derivs(3) * s * s / 2 + derivs(4) * s^3 / 6, ...
     derivs(2) + derivs(3) * s + derivs(4) * s * s / 2, ...
     derivs(3) + derivs(4) * s, ...
     derivs(4)];
end

% =============================================================================== fire control
function [p_hit, cons, sig] = hit_probability(S, x_now, P_now, tau, theta_turret)
n = S.n;
if n >= 3
    a_hat = [x_now(3), x_now(n+3)];
else
    a_hat = [0, 0];
end
px = x_now(1); py = x_now(n+1);
th = atan2(py, px);
nvec = [-sin(th); cos(th)];
Ppos = P_now([1, n+1], [1, n+1]);
sig = sqrt(max(nvec' * Ppos * nvec, 1e-6));
n_cls = size(S.classes, 1);
p_hit = 0; inside = 0;
for i = 1:n_cls
    if i == 1
        w = 0.5; qx = px; qy = py;
    else
        w = 0.5 / (n_cls - 1);
        qx = px + 0.5 * (S.classes(i, 1) - a_hat(1)) * tau * tau;
        qy = py + 0.5 * (S.classes(i, 2) - a_hat(2)) * tau * tau;
    end
    d = hypot(qx, qy) * sin(atan2(qy, qx) - theta_turret);
    p_hit = p_hit + w * (norm_cdf((S.R_ENG - d) / sig) - norm_cdf((-S.R_ENG - d) / sig));
    if abs(d) <= S.R_ENG
        inside = inside + 1;
    end
end
cons = inside / n_cls;
end

function p = norm_cdf(x)
p = 0.5 * (1 + erf(x / sqrt(2)));
end

function th_star = ambush_angle(S, x_upd, gap, th_dot)
% Angle ahead of the target that every maneuver class crosses within the next second.
n = S.n;
if th_dot >= 0, sgn = 1; else, sgn = -1; end
t_settle = 0.3; window = 1.0;
lo = -1e9; hi = 1e9;
for i = 1:size(S.classes, 1)
    for which = 1:2
        if which == 1, tau = gap + t_settle; else, tau = gap + window; end
        [F, ~] = disc_target(S, tau);
        xb = x_upd;
        if i > 1 && n >= 3
            xb(3) = S.classes(i, 1); xb(n+3) = S.classes(i, 2);
            if n >= 4
                xb(4) = 0; xb(n+4) = 0;
            end
        end
        xp = F * xb;
        th = sgn * atan2(xp(n+1), xp(1));
        if which == 1
            lo = max(lo, th);
        else
            hi = min(hi, th);
        end
    end
end
if hi > lo
    th_star = sgn * 0.5 * (lo + hi);
else
    th_star = [];
end
end
