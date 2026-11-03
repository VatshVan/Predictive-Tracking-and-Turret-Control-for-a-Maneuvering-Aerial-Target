"""Guardian predictive-feedback turret controller.

[V_m, fire_cmd] = turret_controller(t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret)

Pipeline (assignment Fig. 2):
  1. Target-state estimation: converted-measurement Kalman filter (unbiased polar->Cartesian
     conversion) with a Markov-jerk motion model (order 4: position..jerk per axis), envelope
     projection (|v| <= 150, |a| <= 20), track-jump detection.
  2. Target prediction: LOS angle derivatives (exact, via complex log-derivatives) and the
     Taylor advance A_3(s) = 1 + tp s + (tp s)^2/2 + (tp s)^3/6 at tp = 0.5 s.
  3. Predicted elevation reference r(t) = A_3(d/dt) theta_hat(t), clamped to the elevation limits.
  4. 2-DOF LTI controller: reference model T(s) = 1/A_3(s) with exact plant-inverse
     feed-forward, plus PI + roll-off feedback C(s). Because T(s) A_3(s) = 1, the turret tracks
     the nowcast theta_hat(t) while the controller only ever sees the tp-predicted reference.
  5. Fire control: Gaussian-mixture hit probability over maneuver classes (envelope-bounded
     acceleration hypotheses), class consensus, recharge timing, fire-at-peak rule.
  6. Blackout handling: coast, maneuver-class re-identification at recovery (T-MPC style),
     optional ambush hold.

All state lives in the module-level dict _S (maps to a MATLAB `persistent` struct). The state
is re-initialised automatically when t goes backwards or t == 0, so consecutive evaluation runs
need no explicit reset.
"""
import math

import numpy as np
from scipy.linalg import expm
from scipy.special import erf

# ----------------------------------------------------------------------------- parameters
TP = 0.50                 # prediction horizon [s] (prescribed)
R_ENG = 2.0               # engagement radius [m]
T_RECHARGE = 1.0          # minimum pulse interval [s]
DT_RADAR = 0.05           # radar interval [s]
SIG_R = 1.0               # range noise [m]
SIG_TH = math.radians(0.02)   # elevation noise [rad]
V_TGT_MAX = 150.0
A_TGT_MAX = 20.0
TH_MIN = math.radians(-10.0)
TH_MAX = math.radians(85.0)
TH_MARGIN = math.radians(0.3)

# plant (theta/V = 1000 / (s (320 s^2 + 3350 s + 1.0015e6)))
FF_C3, FF_C2, FF_C1 = 0.32, 3.35, 1001.5   # exact plant inverse: V = c3 y''' + c2 y'' + c1 y'
V_LIM = 600.0

CFG = dict(
    order=4,              # states per axis: 2 = CV, 3 = Singer accel, 4 = Markov jerk
    sigma_top=20.0,       # std of the top-derivative Markov process (jerk [m/s^3] for order 4)
    tau_top=0.4,          # its correlation time [s]
    q_cv=200.0,           # white-acceleration intensity for order 2 [m^2/s^3]
    tau_fast=0.02,        # extra reference-model pole for orders < 4 (relative degree 3)
    kp=4000.0,            # feedback PI gain [V/rad]
    ti=1.0,               # integral time [s]
    tau_f=0.03,           # feedback roll-off [s]
    p_dev=20.0,           # deviation-dynamics pole (triple) of the reference generator [rad/s]
    omega_ref_max=0.58,   # reference-model rate clamp (acquisition only) [rad/s]
    alpha_ref_max=2.5,    # reference-model acceleration clamp (acquisition only) [rad/s^2]
    use_reid=True,        # maneuver-class re-identification after measurement loss
    use_ambush=False,     # ambush hold during long measurement loss
    p_sure=0.99,          # tracking: fire immediately above this hit probability
    p_fire=0.90,          # tracking: arm above this; fire at the peak or after t_wait
    t_wait=0.30,          # tracking: longest wait for a better shot once armed [s]
    p_min=0.95,           # coasting: fire only while certainty is still near tracking quality
    q_cons=0.80,          # minimum fraction of maneuver classes inside R_ENG
    n_upd_min=8,          # radar updates before the track may fire
    chi2_jump=40.0,       # innovation gate for track-jump detection (2 dof)
    # research options (docs/research_lemma_homotopy.md); defaults reproduce the final design
    ref_mode='deviation',  # 'deviation': T = 1/D_dev with deviation feedback K (final design)
                           # 'exact_taylor': T = 1/A_n(tp s), no K (needs order >= ref_n + 1)
                           # 'moment_matched': T = 1/Q(tp s), Q = A_3 + q4 x^4 + q5 x^5
                           # 'blend': P = (1 - blend_a) A_3 + blend_a d3 (s + p_dev)^3
                           # 'contact': P = ref_poly, agreeing with A_3 to order k (see below)
    ref_n=3,              # Taylor degree n of A_n in 'exact_taylor' mode (3, 4; 5 is unstable)
    mm_q=(0.7414 / 24.0, 0.4697 / 120.0),   # q4, q5 of Q (contact-order study, m = 3, f = 2)
    blend_a=0.5,          # 'blend': weight of the deviation-feedback polynomial (0 = exact, 1 = final)
    ref_poly=(1.0, 0.5, 1.0 / 12.0, 1.0 / 216.0),   # 'contact' P(s), ascending; default (1 + s/6)^3:
                          # contact order 1 at the decay bound N/tp (research_pole_shift.py);
                          # contact order 2: (1, 0.5, 0.125, 1/64, 1/1024) = ((s^2 + 8s + 32)/32)^2
    ex_band_in=0.003,     # exact modes: enter the linear model when |y_d - theta_hat| < this [rad]
    ex_band_out=0.015,    # ... and leave it (saturated acquisition cascade) above this [rad]
    ex_rate_band=0.02,    # ... entry also needs |y_d' - theta_hat'| < this [rad/s]
    phit_mode='mixture',  # 'mixture': 17 constant-acceleration offsets on sigma (final design)
                          # 'homotopy': event-signature fire planner (T-MPC mapping)
    h_c=0.75,             # homotopy: consistency weight of the previous fire plan (T-MPC Eq. 12)
    h_sig_on=10.0,        # homotopy: prior spread of a maneuver onset [m/s^2]
    h_look=True,          # homotopy: plan over the remaining window (False: fire-now P only)
)
DEFAULT_CFG = dict(CFG)        # the final design (analysis scripts may change CFG)

_S = None


def configure(**kw):
    """Change design options (used only by analysis/ablation scripts) and reset the state."""
    CFG.update(kw)
    reset()


def reset():
    global _S
    _S = None


def get_debug():
    return None if _S is None else _S['dbg']


# ----------------------------------------------------------------------------- helpers
def _poly_mul(a, b):
    """Multiply polynomials given as ascending coefficient lists."""
    out = [0.0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def _norm_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _model(n):
    """Continuous per-axis model: chain of integrators, Markov top derivative."""
    A = np.zeros((n, n))
    for k in range(n - 1):
        A[k, k + 1] = 1.0
    if n == 2:
        q = CFG['q_cv']
    else:
        alpha = 1.0 / CFG['tau_top']
        A[n - 1, n - 1] = -alpha
        q = 2.0 * alpha * CFG['sigma_top'] ** 2
    return A, q


def _disc_target(S, dt):
    """Discrete transition F and process noise Q for the 2-axis target model (Van Loan)."""
    key = int(round(dt * 1e9))
    hit = S['disc'].get(key)
    if hit is not None:
        return hit
    n = S['n']
    A, q = S['A'], S['q']
    Qc = np.zeros((n, n))
    Qc[n - 1, n - 1] = q
    M = np.zeros((2 * n, 2 * n))
    M[:n, :n] = -A
    M[:n, n:] = Qc
    M[n:, n:] = A.T
    E = expm(M * dt)
    F1 = E[n:, n:].T
    Q1 = F1 @ E[:n, n:]
    Q1 = 0.5 * (Q1 + Q1.T)
    F = np.kron(np.eye(2), F1)
    Q = np.kron(np.eye(2), Q1)
    if len(S['disc']) < 4000:
        S['disc'][key] = (F, Q)
    return F, Q


def _ref_jerk(S, m, derivs):
    """Third derivative commanded by the reference generator.

    Linear law (small deviations):  D(D) y_d = r_pred - K(D)(y_d - theta_hat), with
    K = D_dev - D and D_dev(s) = d3 (s + p)^3, i.e.
        y''' = phi - 3p d2 - 3p^2 d1 - p^3 d0,     d_j = y^(j) - theta_hat^(j),
        phi  = (r_pred - sum_{j<3} D_j theta_hat^(j)) / d3   (= theta_hat''' for order 4).
    Realised as a cascade (gains p/3, p, 3p) whose intermediate rate/acceleration demands
    are saturated, so large acquisition slews respect the drive limits without chatter.
    """
    D = S['D']
    r_pred = sum(S['A_m'][k] * derivs[k] for k in range(len(S['A_m'])))
    phi = (r_pred - D[0] * derivs[0] - D[1] * derivs[1] - D[2] * derivs[2]) / D[3]
    return _cascade(m, derivs, phi)


def _cascade(m, derivs, phi):
    """Saturated cascade (gains p/3, p, 3p) toward the nowcast; returns y'''."""
    p = CFG['p_dev']
    y0, y1, y2 = m[0], m[1], m[2]
    om, al = CFG['omega_ref_max'], CFG['alpha_ref_max']
    y1_des = min(max(derivs[1] - (p / 3.0) * (y0 - derivs[0]), -om), om)
    y2_des = min(max(derivs[2] - p * (y1 - y1_des), -al), al)
    return phi - 3.0 * p * (y2 - y2_des)


def _ref_top(S, m, derivs, lin):
    """Top derivative y^(N) of the reference generator (N = len(m)).

    'deviation' mode: _ref_jerk. Exact modes: the linear model Q(D) y_d = r_pred with
    r_pred = A_n(D) theta_hat and no deviation feedback, so y_d = A_n/Q theta_hat. 'blend':
    P_a(D) y_d = (1 - a) r_pred + a D_dev(D) theta_hat with P_a = (1 - a) A_3 + a D_dev, the
    weighted average of the exact-Taylor and deviation-feedback laws. In general the input
    polynomial is S['Num'] and the characteristic polynomial S['D']. The linear model is used
    inside a small-deviation band (hysteresis, state in the one-element list `lin`); outside it
    the generator reduces to the saturated third-order acquisition cascade: higher states are
    slaved (y''' = cascade command, y'''' = ... = 0; `m` is modified in place).
    """
    if CFG['ref_mode'] == 'deviation':
        return _ref_jerk(S, m, derivs)
    N = len(m)
    d0, d1 = abs(m[0] - derivs[0]), abs(m[1] - derivs[1])
    if lin[0] and d0 > CFG['ex_band_out']:
        lin[0] = False
    elif not lin[0] and d0 < CFG['ex_band_in'] and d1 < CFG['ex_rate_band']:
        lin[0] = True
    if lin[0]:
        D, Num = S['D'], S['Num']
        u = sum(Num[k] * derivs[k] for k in range(len(Num)))
        return (u - sum(D[j] * m[j] for j in range(N))) / D[N]
    cmd = _cascade(m, derivs, derivs[3])
    if N == 3:
        return cmd
    m[3] = cmd
    for j in range(4, N):
        m[j] = 0.0
    return 0.0


def _taylor(derivs, s):
    """Advance an angle-derivative tuple by s seconds (constant top derivative)."""
    if len(derivs) == 4:
        d0, d1, d2, d3 = derivs
        return (d0 + d1 * s + d2 * s * s / 2.0 + d3 * s ** 3 / 6.0,
                d1 + d2 * s + d3 * s * s / 2.0,
                d2 + d3 * s,
                d3)
    n = len(derivs)
    return tuple(sum(derivs[j + k] * s ** k / math.factorial(k) for k in range(n - j)) for j in range(n))


def _converted(r, th):
    """Unbiased converted measurement and its covariance (Lerro & Bar-Shalom / Mo et al.)."""
    lam = math.exp(-SIG_TH ** 2 / 2.0)
    lam2 = math.exp(-2.0 * SIG_TH ** 2)
    c, s = math.cos(th), math.sin(th)
    z = np.array([r * c, r * s]) / lam
    k = (lam ** -2 - 2.0) * r * r
    h = 0.5 * (r * r + SIG_R ** 2)
    R11 = k * c * c + h * (1.0 + lam2 * math.cos(2 * th))
    R22 = k * s * s + h * (1.0 - lam2 * math.cos(2 * th))
    R12 = k * c * s + h * lam2 * math.sin(2 * th)
    return z, np.array([[R11, R12], [R12, R22]])


def _los(x, n, nd=3):
    """Range, LOS angle and its first nd (3 or 4) time derivatives from the Cartesian state."""
    zs = []
    for k in range(5):
        zs.append(complex(x[k], x[n + k]) if k < n else 0j)
    z0 = zs[0]
    w1, w2, w3, w4 = zs[1] / z0, zs[2] / z0, zs[3] / z0, zs[4] / z0
    th0 = math.atan2(z0.imag, z0.real)
    th1 = w1.imag
    th2 = (w2 - w1 * w1).imag
    th3 = (w3 - 3.0 * w1 * w2 + 2.0 * w1 ** 3).imag
    if nd == 3:
        return abs(z0), th0, th1, th2, th3
    th4 = (w4 - 4.0 * w1 * w3 - 3.0 * w2 * w2 + 12.0 * w1 * w1 * w2 - 6.0 * w1 ** 4).imag
    return abs(z0), th0, th1, th2, th3, th4


def _project_envelope(x, n):
    """Project velocity/acceleration estimates onto the admissible target envelope."""
    for k, lim in ((1, 1.05 * V_TGT_MAX), (2, 1.15 * A_TGT_MAX), (3, 400.0), (4, 2.0e4)):
        if k >= n:
            break
        mag = math.hypot(x[k], x[n + k])
        if mag > lim:
            x[k] *= lim / mag
            x[n + k] *= lim / mag
    return x


def _classes():
    """Maneuver classes: nominal estimate plus envelope-bounded acceleration hypotheses."""
    acc = [(0.0, 0.0)]
    for mag in (0.5 * A_TGT_MAX, A_TGT_MAX):
        for k in range(8):
            ang = k * math.pi / 4.0
            acc.append((mag * math.cos(ang), mag * math.sin(ang)))
    return acc


# ----------------------------------------------------------------------------- init
def _init(theta, theta_dot):
    n = CFG['order']
    A, q = _model(n)
    if CFG['ref_mode'] == 'deviation':
        m = n - 1                               # Taylor-advance degree
        A_m = [TP ** k / math.factorial(k) for k in range(m + 1)]
        D = A_m                                 # prefilter denominator, relative degree 3
        for _ in range(3 - m):
            D = _poly_mul(D, [1.0, CFG['tau_fast']])
    else:
        m = CFG['ref_n'] if CFG['ref_mode'] == 'exact_taylor' else 3
        if n < m + 1:
            raise ValueError(f"ref_mode {CFG['ref_mode']} with n = {m} needs order >= {m + 1}")
        A_m = [TP ** k / math.factorial(k) for k in range(m + 1)]
        D = list(A_m)                           # Q(tp s), ascending powers of s
        if CFG['ref_mode'] == 'moment_matched':
            D += [qk * TP ** (4 + k) for k, qk in enumerate(CFG['mm_q'])]
    Num = list(A_m)
    if CFG['ref_mode'] == 'contact':
        D = [float(c) for c in CFG['ref_poly']]
        Num = D[:4]                             # P y_d = A_3 theta_hat + (P - A_3) theta_hat; no
                                                # theta_hat'''' from an order-4 estimator
    elif CFG['ref_mode'] == 'blend':
        p, a = CFG['p_dev'], CFG['blend_a']
        D_dev = [A_m[3] * c for c in (p ** 3, 3.0 * p * p, 3.0 * p, 1.0)]   # d3 (s + p)^3
        D = [(1.0 - a) * q + a * d for q, d in zip(A_m, D_dev)]
        Num = list(D)
    n_ref = max(3, len(D) - 1)                  # reference-generator states y_d .. y_d^(N-1)
    nd = max(3, len(A_m) - 1)                   # LOS derivatives used by the predictor
    H = np.zeros((2, 2 * n))
    H[0, 0] = 1.0
    H[1, n] = 1.0
    return dict(
        n=n, A=A, q=q, H=H, A_m=A_m, D=D, Num=Num, nd=nd, disc={},
        classes=_classes(),
        t_prev=0.0, first=True,
        # estimator
        trk=False, x=None, P=None, t_upd=-1.0, n_upd=0, last_meas=(None, None), n_jump=0,
        n_proj=0,
        # reference model / feedback
        m_ref=[theta, theta_dot] + [0.0] * (n_ref - 2), derivs_prev=None, lin=[False],
        xi=0.0, xf=0.0, sat=0,
        # fire control
        t_fire=-1e9, p_prev=0.0, t_armed=None, p_armed_max=0.0, ambush=None,
        h_plan=None, n_wait=0,
        dbg={},
    )


# ----------------------------------------------------------------------------- estimator
def _kf_update(S, z, R):
    H = S['H']
    P = S['P']
    nu = z - H @ S['x']
    Sm = H @ P @ H.T + R
    K = P @ H.T @ np.linalg.inv(Sm)
    S['x'] = S['x'] + K @ nu
    I_KH = np.eye(P.shape[0]) - K @ H
    S['P'] = I_KH @ P @ I_KH.T + K @ R @ K.T
    x_raw = S['x'].copy()
    S['x'] = _project_envelope(S['x'], S['n'])
    if not np.array_equal(x_raw, S['x']):
        S['n_proj'] += 1


def _track_init(S, z, R, t):
    n = S['n']
    x = np.zeros(2 * n)
    x[0], x[n] = z[0], z[1]
    P = np.zeros((2 * n, 2 * n))
    P[np.ix_([0, n], [0, n])] = R
    sig0 = [None, V_TGT_MAX, A_TGT_MAX, CFG['sigma_top'] if n == 4 else 60.0, CFG['sigma_top']]
    for k in range(1, n):
        P[k, k] = P[n + k, n + k] = sig0[k] ** 2
    S.update(trk=True, x=x, P=P, t_upd=t, n_upd=1, n_jump=0)


def _reidentify(S, z, R, F, Q):
    """Maneuver-class re-identification after measurement loss.

    Each class holds a constant acceleration (envelope-bounded) through the gap; the nominal
    class is the model propagation. The class with the best prior-weighted likelihood seeds
    the state before the normal update. Prior weights play the role of the consistency
    weight in T-MPC (de Groot et al., Eq. 12): the nominal hypothesis is preferred unless the
    data clearly favour a maneuver.
    """
    n = S['n']
    x_b = S['x']
    Pp = F @ S['P'] @ F.T + Q
    H = S['H']
    Sm = H @ Pp @ H.T + R
    Sinv = np.linalg.inv(Sm)
    n_cls = len(S['classes'])
    best_cost, best_x = None, None
    for i, (ax, ay) in enumerate(S['classes']):
        xb = x_b.copy()
        if i == 0:
            prior = 0.5
        else:
            prior = 0.5 / (n_cls - 1)
            if n >= 3:
                xb[2], xb[n + 2] = ax, ay
                for k in range(3, n):
                    xb[k] = xb[n + k] = 0.0
            else:
                continue
        xi = F @ xb
        nu = z - H @ xi
        cost = 0.5 * float(nu @ Sinv @ nu) - math.log(prior)
        if best_cost is None or cost < best_cost:
            best_cost, best_x = cost, xi
    S['x'] = best_x
    S['P'] = Pp


def _estimator_step(S, t, r_meas, th_meas, valid):
    if not valid or r_meas is None or not (r_meas == r_meas) or not (th_meas == th_meas):
        return
    if (r_meas, th_meas) == S['last_meas']:
        return                                   # held duplicate, not a new measurement
    S['last_meas'] = (r_meas, th_meas)
    z, R = _converted(r_meas, th_meas)
    if not S['trk']:
        _track_init(S, z, R, t)
        return
    gap = t - S['t_upd']
    if gap <= 0.0:
        return
    F, Q = _disc_target(S, gap)
    x_pred = F @ S['x']
    P_pred = F @ S['P'] @ F.T + Q
    H = S['H']
    nu = z - H @ x_pred
    Sm = H @ P_pred @ H.T + R
    d2 = float(nu @ np.linalg.solve(Sm, nu))
    if d2 > CFG['chi2_jump']:
        S['n_jump'] += 1
        if S['n_jump'] >= 2:                     # persistent jump: new target, restart track
            _track_init(S, z, R, t)
        return
    S['n_jump'] = 0
    if gap > 1.5 * DT_RADAR and CFG['use_reid']:
        _reidentify(S, z, R, F, Q)
    else:
        S['x'], S['P'] = x_pred, P_pred
    _kf_update(S, z, R)
    S['t_upd'] = t
    S['n_upd'] += 1


# ----------------------------------------------------------------------------- fire control
def _hit_probability(S, x_now, P_now, tau, theta_turret):
    """Gaussian mixture over maneuver classes of the beam/target cross-range offset."""
    n = S['n']
    a_hat = (x_now[2], x_now[n + 2]) if n >= 3 else (0.0, 0.0)
    px, py = x_now[0], x_now[n]
    rng = math.hypot(px, py)
    th = math.atan2(py, px)
    nvec = np.array([-math.sin(th), math.cos(th)])
    Ppos = P_now[np.ix_([0, n], [0, n])]
    sig = math.sqrt(max(float(nvec @ Ppos @ nvec), 1e-6))
    n_cls = len(S['classes'])
    p_hit, inside = 0.0, 0
    for i, (ax, ay) in enumerate(S['classes']):
        w = 0.5 if i == 0 else 0.5 / (n_cls - 1)
        if i == 0:
            qx, qy = px, py
        else:
            qx = px + 0.5 * (ax - a_hat[0]) * tau * tau
            qy = py + 0.5 * (ay - a_hat[1]) * tau * tau
        d = math.hypot(qx, qy) * math.sin(math.atan2(qy, qx) - theta_turret)
        p_hit += w * (_norm_cdf((R_ENG - d) / sig) - _norm_cdf((-R_ENG - d) / sig))
        if abs(d) <= R_ENG:
            inside += 1
    return p_hit, inside / n_cls, sig, rng


# ----------------------------------------------------------------------------- homotopy planner
# T-MPC mapping (docs/research_lemma_homotopy.md). State-time plane (t, e), e = signed
# cross-range offset of the target from the turret's predicted beam path. Obstacles are the
# candidate firing events O_j = {t_j} x [-R_ENG, R_ENG]. A target future is a graph e(t), so its
# winding around each O_j is fixed by the side on which it passes the event: its homotopy class
# (H-signature) is h_j in {+1 above, 0 through = hit, -1 below}. The decision is a fire plan
# (one event); its cost is the posterior mass of futures that do not pass through that event,
# and plans are selected with the consistency rule of T-MPC Eq. 12.
H_ABOVE, H_THROUGH, H_BELOW = 1, 0, -1
H_SHAPES = ('constant', 'onset', 'reverse')


def _h_signature(e):
    """Homotopy signature of offsets e[..., j] at the firing events."""
    e = np.asarray(e)
    return np.where(e > R_ENG, H_ABOVE, np.where(e < -R_ENG, H_BELOW, H_THROUGH))


def _h_hypotheses(S, x, P):
    """Guidance hypotheses: admissible maneuver futures about the estimate x.

    Nominal (prior 0.5) plus 8 directions x 3 shapes (constant from the last update, onset at
    mid-horizon, reversing at mid-horizon), each at the largest magnitude keeping
    |a_hat + da| <= A_TGT_MAX. Posterior weight = prior x likelihood of da under the
    acceleration uncertainty P_a plus a maneuver-onset spread. Returns w (K,), da (K, 2), shape (K,).
    """
    n = S['n']
    if n >= 3:
        a = np.array([x[2], x[n + 2]])
        Sig = P[np.ix_([2, n + 2], [2, n + 2])] + CFG['h_sig_on'] ** 2 * np.eye(2)
    else:
        a = np.zeros(2)
        Sig = CFG['h_sig_on'] ** 2 * np.eye(2)
    Sinv = np.linalg.inv(Sig)
    a2 = float(a @ a)

    def cap(au):
        return max(-au + math.sqrt(max(au * au - a2 + A_TGT_MAX ** 2, 0.0)), 0.0)

    w, da, shp = [0.5], [(0.0, 0.0)], [0]
    n_man = 8 * len(H_SHAPES)
    for k in range(8):
        u = np.array([math.cos(k * math.pi / 4.0), math.sin(k * math.pi / 4.0)])
        au = float(a @ u)
        quad = float(u @ Sinv @ u)
        for sh in range(len(H_SHAPES)):
            lam = min(cap(au), cap(-au)) if sh == 2 else cap(au)
            w.append(0.5 / n_man * math.exp(-0.5 * lam * lam * quad))
            da.append((lam * u[0], lam * u[1]))
            shp.append(sh)
    w = np.array(w)
    return w / w.sum(), np.array(da), np.array(shp)


def _h_disp(s, horizon, shp):
    """Displacement per unit acceleration after s seconds, per hypothesis shape (vectorised)."""
    s1 = 0.5 * horizon
    const = 0.5 * s * s
    onset = 0.5 * max(s - s1, 0.0) ** 2
    rev = const if s <= s1 else 0.5 * s1 * s1 + s1 * (s - s1) - 0.5 * (s - s1) ** 2
    return np.where(shp == 0, const, np.where(shp == 1, onset, rev))


def _beam_path(S, theta_turret, taus):
    """Predicted beam angle at t + tau: reference generator rolled forward (10 ms steps) on the
    extrapolated nowcast, offset by the present turret error."""
    beam = [theta_turret]
    if len(taus) == 1:
        return beam
    h = 0.01
    m = list(S['m_ref'])
    lin = [S['lin'][0]]
    d0 = S['derivs_prev']
    off = theta_turret - m[0]
    n_step = int(round(taus[-1] / h))
    for k in range(1, n_step + 1):
        top = _ref_top(S, m, _taylor(d0, (k - 1) * h), lin)
        m[-1] += h * top
        for j in range(len(m) - 2, -1, -1):
            m[j] += h * m[j + 1]
        if k % 5 == 0:
            beam.append(m[0] + off)
    return beam[:len(taus)]


def _h_plan(S, t, theta_turret, horizon):
    """Probability that the target passes through each firing event t + tau_j (tau_0 = 0: now).

    Hypotheses start at the last update (gap g), so the maneuver spread over the gap and the
    look-ahead is carried by the hypotheses; each keeps a Gaussian position spread from the
    noise-free propagation F P F' (Rao-Blackwellised, no process noise: no offsets stacked on
    sigma). Returns taus, P_through (J,), signatures (K, J), weights (K,), sigma_0, range.
    """
    n = S['n']
    g = max(t - S['t_upd'], 0.0)
    n_ev = int(round(horizon / 0.05))
    taus = [0.05 * j for j in range(n_ev + 1)]
    beam = _beam_path(S, theta_turret, taus)
    w, da, shp = _h_hypotheses(S, S['x'], S['P'])
    H = g + taus[-1]
    E = np.zeros((len(w), len(taus)))
    P_thr = np.zeros(len(taus))
    sig0 = rng0 = 0.0
    for j, tau in enumerate(taus):
        s = g + tau
        F, _ = _disc_target(S, s)
        xm = F @ S['x']
        Pm = F @ S['P'] @ F.T
        px, py = xm[0], xm[n]
        th = math.atan2(py, px)
        nvec = np.array([-math.sin(th), math.cos(th)])
        sig = math.sqrt(max(float(nvec @ Pm[np.ix_([0, n], [0, n])] @ nvec), 1e-6))
        disp = _h_disp(s, H, shp)
        qx = px + da[:, 0] * disp
        qy = py + da[:, 1] * disp
        e = qy * math.cos(beam[j]) - qx * math.sin(beam[j])
        E[:, j] = e
        p = 0.5 * (erf((R_ENG - e) / (sig * math.sqrt(2.0))) - erf((-R_ENG - e) / (sig * math.sqrt(2.0))))
        P_thr[j] = float(w @ p)
        if j == 0:
            sig0, rng0 = sig, math.hypot(px, py)
    return taus, P_thr, _h_signature(E), w, sig0, rng0


def _ambush_angle(S, x_upd, gap, th_dot):
    """Angle ahead of the target that every maneuver class crosses within the next second."""
    n = S['n']
    sgn = 1.0 if th_dot >= 0.0 else -1.0
    t_settle, window = 0.3, 1.0
    lo, hi = -1e9, 1e9
    for i, (ax, ay) in enumerate(S['classes']):
        for tau, which in ((gap + t_settle, 'lo'), (gap + window, 'hi')):
            F, _ = _disc_target(S, tau)
            xb = x_upd.copy()
            if i > 0 and n >= 3:
                xb[2], xb[n + 2] = ax, ay
                for k in range(3, n):
                    xb[k] = xb[n + k] = 0.0
            xp = F @ xb
            th = sgn * math.atan2(xp[n], xp[0])
            if which == 'lo':
                lo = max(lo, th)
            else:
                hi = min(hi, th)
    if hi > lo:
        return sgn * 0.5 * (lo + hi)
    return None


# ----------------------------------------------------------------------------- main entry
def turret_controller(t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret):
    global _S
    if _S is None or t < _S['t_prev'] - 1e-12 or t == 0.0:
        _S = _init(theta_turret, theta_dot_turret)
    S = _S
    dt = 0.0 if S['first'] else t - S['t_prev']
    S['first'] = False
    n = S['n']

    # 1. advance the reference generator over the last interval; inside it, the nowcast
    #    derivatives evolve by their own Taylor model (sub-steps <= 1 ms)
    if dt > 0.0 and S['derivs_prev'] is not None:
        n_sub = max(1, int(math.ceil(dt / 1e-3 - 1e-9)))
        h = dt / n_sub
        m = list(S['m_ref'])
        for k in range(n_sub):
            top = _ref_top(S, m, _taylor(S['derivs_prev'], k * h), S['lin'])
            m[-1] += h * top
            for j in range(len(m) - 2, -1, -1):
                m[j] += h * m[j + 1]
        S['m_ref'] = m

    # 2. estimation
    _estimator_step(S, t, r_meas, theta_meas, meas_valid)

    # 3. nowcast + tp-prediction
    mode = 'acquire'
    y, y1, y2 = S['m_ref'][:3]
    zeros = (0.0,) * S['nd']
    derivs = (y,) + zeros                        # before track: hold the current reference
    r_pred = y
    p_hit, cons, sig, rng_hat, th_hat = 0.0, 0.0, float('nan'), float('nan'), float('nan')
    if S['trk']:
        tau = max(t - S['t_upd'], 0.0)
        F, Q = _disc_target(S, tau)
        x_now = F @ S['x']
        P_now = F @ S['P'] @ F.T + Q
        los = _los(x_now, n, S['nd'])
        rng_hat, th_hat = los[0], los[1]
        derivs = tuple(los[1:])
        mode = 'track'
        if tau > 1.5 * DT_RADAR:
            mode = 'coast'
            if CFG['use_ambush'] and tau > 0.25:
                if S['ambush'] is None:
                    S['ambush'] = _ambush_angle(S, S['x'], tau, los[2])
                if S['ambush'] is not None:
                    derivs = (S['ambush'],) + zeros
                    mode = 'ambush'
        else:
            S['ambush'] = None
        if CFG['phit_mode'] == 'homotopy':
            _, P_thr, sgn, w_h, sig, _ = _h_plan(S, t, theta_turret, 0.0)
            p_hit = float(P_thr[0])
            cons = float(w_h @ (sgn[:, 0] == H_THROUGH))   # posterior mass of 'through' futures
        else:
            p_hit, cons, sig, _ = _hit_probability(S, x_now, P_now, tau, theta_turret)
    if not (TH_MIN + TH_MARGIN <= derivs[0] <= TH_MAX - TH_MARGIN):
        derivs = (min(max(derivs[0], TH_MIN + TH_MARGIN), TH_MAX - TH_MARGIN),) + zeros
    # predicted elevation reference: Taylor advance of the nowcast by tp
    r_pred = sum(S['A_m'][k] * derivs[k] for k in range(len(S['A_m'])))
    S['derivs_prev'] = derivs

    # 4-5. reference generator + 2-DOF control: exact plant-inverse feed-forward + PI feedback
    if len(S['m_ref']) == 3:
        y3 = _ref_top(S, S['m_ref'], derivs, S['lin'])
    else:
        y3 = S['m_ref'][3]
    u_ff = FF_C3 * y3 + FF_C2 * y2 + FF_C1 * y1
    r_ref = r_pred
    e = y - theta_turret
    if dt > 0.0:
        if not (S['sat'] != 0 and S['sat'] * e > 0):
            S['xi'] += e * dt                    # conditional integration (anti-windup)
        a = math.exp(-dt / CFG['tau_f'])
        S['xf'] = a * S['xf'] + (1.0 - a) * CFG['kp'] * (e + S['xi'] / CFG['ti'])
    u = u_ff + S['xf']
    if u > V_LIM:
        V, S['sat'] = V_LIM, 1
    elif u < -V_LIM:
        V, S['sat'] = -V_LIM, -1
    else:
        V, S['sat'] = u, 0

    # 6. fire control: fire-now-or-wait over the next recharge window
    fire = 0
    ready = (t - S['t_fire'] >= T_RECHARGE + 1e-6) and S['trk'] and S['n_upd'] >= CFG['n_upd_min']
    if ready and cons >= CFG['q_cons'] and CFG['phit_mode'] == 'homotopy' and CFG['h_look']:
        # fire plans = firing events within the remaining window (one recharge period from first
        # arming: waiting longer would forfeit a pulse); cost J_j = 1 - P(through event j);
        # T-MPC Eq. 12: the previously selected plan, re-identified by its absolute fire time,
        # has its cost scaled by c < 1
        floor = CFG['p_min'] if mode in ('coast', 'ambush') else CFG['p_fire']
        if mode == 'track' and p_hit >= CFG['p_sure']:
            fire = 1
        elif p_hit >= floor:
            if S['t_armed'] is None:
                S['t_armed'] = t
            W = T_RECHARGE - (t - S['t_armed'])
            if W <= DT_RADAR:
                fire = 1
            else:
                taus, P_thr, _, _, _, _ = _h_plan(S, t, theta_turret, min(W, T_RECHARGE))
                cost = 1.0 - P_thr
                if S['h_plan'] is not None:
                    jp = int(np.argmin([abs(t + tj - S['h_plan']) for tj in taus]))
                    if abs(t + taus[jp] - S['h_plan']) <= 0.03:
                        cost[jp] *= CFG['h_c']
                j = int(np.argmin(cost))
                if j == 0:
                    fire = 1
                else:
                    S['h_plan'] = t + taus[j]
                    S['n_wait'] += 1
        else:
            S['t_armed'] = None
            S['h_plan'] = None
    elif ready and cons >= CFG['q_cons']:
        if mode in ('coast', 'ambush'):
            # no new data: certainty can only decay, so shoot while it is still high
            if p_hit >= CFG['p_min'] and p_hit < S['p_prev']:
                fire = 1
        elif p_hit >= CFG['p_sure']:
            fire = 1
        elif p_hit >= CFG['p_fire']:
            if S['t_armed'] is None:
                S['t_armed'], S['p_armed_max'] = t, p_hit
            S['p_armed_max'] = max(S['p_armed_max'], p_hit)
            if p_hit < S['p_armed_max'] - 0.02 or t - S['t_armed'] >= CFG['t_wait']:
                fire = 1                         # peak passed, or waited long enough
        else:
            S['t_armed'] = None
    if fire:
        S['t_fire'] = t
        S['t_armed'] = None
        S['h_plan'] = None
    S['p_prev'] = p_hit

    S['t_prev'] = t
    S['dbg'] = dict(t=t, mode=mode, th_hat=th_hat, r_hat=rng_hat, r_ref=r_ref, y_d=y,
                    u_ff=u_ff, u_fb=S['xf'], V=V, p_hit=p_hit, cons=cons, sig=sig,
                    n_proj=S['n_proj'], lin=S['lin'][0], n_wait=S['n_wait'],
                    x_est=(S['x'][0], S['x'][n]) if S['trk'] else (float('nan'),) * 2)
    return V, fire
