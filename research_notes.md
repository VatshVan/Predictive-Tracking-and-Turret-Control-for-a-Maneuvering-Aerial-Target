# ME 319 Assignment 2: Research Notes (Mumbai 2044 / V-44)

Deadline: **Nov 10, 2026, 11:59 PM**. Deliverables: ≤6-page IEEE paper + 1-page AI appendix, public GitHub repo, peer-evaluation form.

---

## 1. Hard constraints from the spec that change the plan

| Spec says | Consequence |
|---|---|
| "The controller shall be formulated as a **conventional linear time-invariant (LTI) compensator** … supported by analytical control-system considerations rather than numerical tuning alone." (§4) | **MPPI cannot be the turret feedback controller.** G_c(s) must be LTI (lead/lag/PID/2-DOF), with root locus, Bode, and margins. Sampling ideas have to live upstream, in prediction, aim-point selection or fire decision. |
| Evaluation uses fixed **MATLAB signatures** `target_motion_generator`, `turret_controller`, `turret_plant`; repo must contain `turret_controller.m` and `turret_plant.m` (§3, §9) | **Final Evaluation is 5/15 points** and runs through the instructor's MATLAB script. Python is fine for prototyping, but the submitted controller must be `.m`. MATLAB R2025a is installed on this machine, so I'd write it in MATLAB from the start. |
| Controller is a stateless-looking function called each step | All filter/controller state lives in `persistent` variables inside `turret_controller.m`. We don't know the evaluation step size, so the controller should measure its own dt from `t`. |
| "designed once … without retuning for individual cases" | One parameter set for A/B/C and the hidden cases. Tune for the *envelope* (|v|≤150, |a|≤20, 800–2000 m, blackout ≤0.5 s at t≥1 s), not the scenarios. |
| "Preference is given to designs that are simple, well-motivated and analytically supported" | Originality should come from *justified adaptation*, not stacking algorithms. Every block needs a one-line "why this, for this problem" argument. |
| Rubric: Originality 2, Estimation and Prediction 2, Modelling 1.5, Controller 1.5, Sim 1.5, **Final Eval 5**, Writing 1.5 | Hit performance on unseen cases dominates. The fire-decision logic is as important as tracking. |

---

## 2. Turret plant: derived numbers (verified in Python)

Reflect the motor to the turret side (gear ratio N = 100):

- J_eq = J_turret + N²·J_m = 1200 + 10⁴·0.2 = **3200 kg·m²**
- Torque at turret = N·K_t·i = 1000·i; back-EMF = K_e·N·θ̇ = 1000·θ̇

```
L di/dt   = V − R i − N K_e θ̇
J_eq θ̈   = N K_t i − b θ̇

θ(s)/V(s) = N K_t / ( s [ (L s + R)(J_eq s + b) + N² K_t K_e ] )
          = 1000 / ( s (320 s² + 3350 s + 1.0015e6) )
```

- Complex pole pair: **ω_n ≈ 55.9 rad/s (8.9 Hz), ζ ≈ 0.094**. The plant is lightly damped and electromechanically resonant because the large back-EMF term N²K_tK_e couples the L/R electrical lag (τ_e = 0.1 s) with the very fast mechanical pole (J R / N²K_tK_e = 3.2 ms).
- Plus a free integrator. Type-1 plant.
- Velocity DC gain ≈ 1.0×10⁻³ (rad/s)/V. Tracking a 0.19 rad/s LOS rate (150 m/s at 800 m) needs ≈ **190 V** steady-state. No voltage limit is stated. Ask whether there is one; if not, state an assumed limit in the paper.
- Implications for G_c: a notch or a lead network is needed around 8.9 Hz, and the reference must be **smooth**. Step changes from 20 Hz radar updates or the blackout recovery will ring the resonance. The controller must run faster than 20 Hz. The resonance (8.9 Hz) sits just under the radar Nyquist (10 Hz).

---

## 3. Error budget: the numbers that drive the whole design

| Quantity | Value | Meaning |
|---|---|---|
| Hit tolerance R_eng | 2 m | **1.0 mrad at 2000 m**, 2.5 mrad at 800 m |
| Radar σ_θ | 0.02° = **0.35 mrad** | Cross-range σ ≈ 0.7 m at 2 km per raw sample |
| Radar σ_r | 1 m | |
| Unknown-acceleration drift over t_p = 0.5 s | ½·20·0.5² = **2.5 m** | Worst case > R_eng. Pure prediction 0.5 s ahead *cannot* guarantee a hit against an adversarial maneuver. |
| Same drift over one radar interval (50 ms) | 0.025 m | Negligible. Short-horizon prediction is easy. |
| Blackout drift (0.5 s) | up to 2.5 m (plus velocity error × 0.5 s) | During loss, P(hit) genuinely drops. Fire logic should know this. |

**Key observation (good paper material):** *laser propagation time is negligible*. A hit is judged at the instant of the pulse, between the beam and the target's **current** position. Unlike a gun, no physical lead angle is needed. Prediction exists only to beat the **turret's closed-loop lag** and to schedule the shot. Two consistent ways to read the spec's t_p = 0.5 s requirement:

1. **Lag compensation.** The reference is θ̂(t + τ), where τ is matched to the closed-loop phase lag, so the turret's *actual* angle coincides with the target's *current* angle. t_p sets the predictor horizon, and the controller uses the predicted trajectory segment as feed-forward.
2. **Ambush / pre-positioning** (the 1944 picture in Fig. 1). Drive the turret to θ̂(t + t_p), let it settle, and fire when the estimated target crosses the beam.

We need to choose deliberately and say why. Option 1 gives continuous firing opportunities. Option 2 gives a well-settled turret. This interpretation choice is a real originality point; most teams will just point at θ(t+0.5) and miss.

---

## 4. Target-state estimation: EKF vs UKF vs PF vs CMKF

Measurements are (r, θ) at 20 Hz, polar, in a 2-D plane. State in Cartesian (x, y, ẋ, ẏ, [ẍ, ÿ]).

- **Nonlinearity is tiny here.** The relevant measure is cross-range error r·σ_θ ≈ 0.7 m versus σ_r = 1 m. Both are small compared with range, and r·σ_θ² ≈ 2.4×10⁻⁴ m. The polar→Cartesian conversion bias that Lerro & Bar-Shalom warn about only matters when r·σ_θ²/σ_r is not small. Here it is ~10⁻⁴, so **EKF, UKF and converted-measurement KF (CMKF) will perform almost identically.**
- **Particle filter: not justified.** The posterior is essentially Gaussian, the state is 4–6-D, and the rubric favours simple analytically supported designs. A PF adds cost and noise for no gain. One sentence in the paper explaining why it was rejected is worth more than using it.
- **Recommendation:** a **converted-measurement linear KF** (Lerro & Bar-Shalom debiased conversion, with covariance computed at the predicted state per Bordonaro/Willett/Bar-Shalom) **or** an EKF in mixed coordinates. The CMKF keeps the filter linear, which makes the gain/bandwidth analysis LTI-friendly and analysable with the same tools as the controller. EKF is the "default" choice. Either is defensible. I lean CMKF because the steady-state gain can be analysed (α-β-γ equivalence), which ties estimation to control-theory content.
- Range-rate is **not** measured, so velocity comes purely from the dynamics model.

From your LiDAR/GPS background: the change is that radar error is anisotropic and polar, tight in angle and looser in range, and its Cartesian ellipse rotates with θ. Recomputing R_cart every step from the current θ handles this.

## 5. Motion model and maneuvers

| Model | Fit to this problem |
|---|---|
| CV (constant velocity) | Scenario A only. Lags in B. |
| CA (constant acceleration) | Handles the steady part of a maneuver, but is noisy. |
| **Singer** (exp.-correlated accel, time constant τ_m, variance tied to a_max) | Natural fit: the spec *gives* a_max = 20 m/s². Setting σ_a from a_max makes Q *derived from the threat envelope, not hand-tuned*. That is exactly the "analytically supported" argument. |
| Current Statistical Model (Zhou) | Adaptive mean acceleration with Rayleigh bounded by a_max. Also envelope-driven. |
| Coordinated turn (CT) | Good for turning flight. The V-44 maneuvers are "bounded in x and y", not necessarily turns. |
| **IMM** (e.g. CV + Singer/CA, Markov switching) | The standard answer for maneuvering targets (Mazor et al. 1998; Li & Jilkov 2005). Pitre/Jilkov/Li found IMM gives the best cost/accuracy trade-off. It is **prior art**, so it needs adaptation. |

**Adaptation ideas that make it ours (pick 1–2, not all):**
- **Envelope-constrained estimate/prediction.** After each update, project the estimate onto the admissible set (|v| ≤ 150, |a| ≤ 20, 800 ≤ r ≤ 2000). The spec guarantees these bounds, so a prediction that violates them is known to be wrong. This is cheap, provably helps, and uses information most filters ignore.
- **Q derived from a_max** (Singer) instead of tuned. Justify τ_m from "smooth bounded maneuver".
- **IMM with a non-maneuver mode and a maneuver mode** whose process noise is set by the a_max envelope. Mode probability μ_maneuver also feeds the fire decision (see §7).

## 6. Radar loss (Scenario C)

- During loss: run predict-only (coast). Covariance grows naturally. Freeze IMM mode probabilities (or let them drift via the Markov matrix). This is standard.
- The parts that go beyond standard:
  - **Bound the coasted state** with the envelope projection. Covariance growth can also be capped by the reachable-set bound (½·a_max·t² position), since the true uncertainty cannot exceed what physics allows.
  - **Recovery transient.** The first measurement after loss can produce a large innovation. Don't gate it out. The KF gain will be large because P grew, so the estimate jumps. **Rate-limit or smooth the reference** at recovery so the jump doesn't excite the 8.9 Hz turret mode. The spec explicitly asks to analyse "transient response following measurement recovery".
  - **Fire inhibit** while the predicted miss-uncertainty exceeds R_eng.
- `meas_valid` tells us when data is missing, so loss detection is given, not estimated.

## 7. Where your MPPI / sampling idea fits (and where the T-MPC paper helps)

Since G_c must be LTI, put the sampling where it is legitimate and actually useful: **the fire-control decision and aim-point selection.**

**Sampling-based engagement decision.** Each control step:
1. Draw M trajectories of target motion over the next short window from the filter posterior. Sample accelerations as smooth random processes bounded by a_max (MPPI-style rollouts of the *target*, not the turret).
2. For candidate fire times (now, +Δ, …) evaluate the beam–target miss distance using the turret's predicted angle from the known closed-loop model.
3. Estimated P(hit) = fraction of rollouts with miss ≤ R_eng. **Fire on the 0→1 transition when P(hit) ≥ threshold and recharge allows.** Since only the first hit counts for interception time and pulses are limited to 1 per second, a wasted pulse costs ≥1 s.

This is the honest home for MPPI's "simulate many futures, weight by cost" idea. It is original in this setting and easy to justify.

**What transfers from T-MPC (de Groot et al., 2024):**
- T-MPC runs **several distinct hypotheses in parallel** (homotopy classes) and picks the lowest-cost one. The homotopy/obstacle machinery itself does **not** transfer: there are no obstacles. What maps is **distinct maneuver hypotheses**: climb/dive/accelerate/brake as separate predicted futures, i.e. the IMM modes or bounded-maneuver extremes. Choose the aim point that minimises the worst-case or expected miss across them (a minimax aim).
- Their **consistency rule** (Eq. 12, weight c_i < 1 on the previously chosen plan) maps directly to **hysteresis in aim-point / mode selection**. That matters here because switching the reference excites the lightly damped turret resonance. They used c_i = 0.75 and showed both extremes are worse. That is a nice precedent to cite.
- Cite it as the inspiration for "parallel hypotheses + consistent selection". Don't claim we're doing T-MPC.

## 8. Feedback controller G_c(s)

- Plant: integrator × lightly damped 2nd-order (ω_n 55.9 rad/s, ζ 0.094).
- Reference is a ramp/parabola in angle (LOS rate and acceleration), so a type-1 plant gives a constant ramp error. Options: add an integrator (type-2 loop, zero ramp error), or a **2-DOF structure with feed-forward from the predicted LOS rate/accel**. The predictor already outputs θ̇̂ and θ̈̂, which makes feed-forward natural. A feed-forward + feedback structure is still LTI.
- Candidate: lead-lag or PID + notch/low-pass near 8.9 Hz, crossover well below the resonance (≈1–2 Hz), designed with root locus / Bode and showing GM/PM. Noise from 0.35 mrad radar must not pass through to the motor voltage. Discuss bandwidth vs noise.
- Respect −10°…+85° elevation limits (reference saturation + anti-windup if integral action).
- The controller is designed in continuous time and then discretised (Tustin) at the controller rate. Discuss discretisation.

## 9. What real systems do (for the "prior art" paragraph)

- WWII Kerrison/M9: analog computation of lead from a constant-velocity assumption.
- Modern gun fire-control computers use a **Kalman target-state estimator** for position/velocity/acceleration, from which lead angles are computed (DTIC anti-aircraft fire-control report). Patents note the constant-velocity/linear assumption fails on maneuvering or weaving targets.
- Naval gun fire control has used **IMM with CV + Singer + CT** models.
- High-energy-laser C-UAS chains: wide-FOV acquisition → narrow-FOV fine tracking (fast steering mirrors) → aim-point selection → **aim-point maintenance / dwell**. Our problem abstracts to a single pulse with instantaneous effect, so dwell is not modelled. Say so in the assumptions.

## 10. Open questions for you

1. **Voltage limit on V_m?** Not in the spec. ~190 V is needed for the worst-case LOS rate.
2. **Evaluation step size**: will their script call the controller at 20 Hz, or faster (e.g. 1 kHz with `meas_valid=0` between radar samples)? This changes the controller discretisation a lot. Check the course page/Piazza, or ask the TA.
3. **Initial turret angle and target start position.** The first ~1 s is acquisition. Interception time is scored, so how fast we slew and settle matters.
4. **Team size and split.** 1–5 students. If there are teammates, the natural split is estimator/predictor, controller design and analysis, fire logic + simulation harness, paper.
5. **t_p interpretation** (§3). I'd like us to decide this together because it shapes everything downstream.

---

## 11. Round-2 decisions and proposals (Nov 2)

**Team:** solo. **Deadline:** Nov 10, 2026 (per spec). Papers are in `papers/` (see `papers/README.md`).

### Workflow: Python first, MATLAB at the end
- Develop and analyse in Python. Write `turret_controller` / `turret_plant` / `target_motion_generator` in a **MATLAB-portable subset**: same signatures, all state in one dict (maps to `persistent` struct), plain numpy linear algebra, no classes/closures in the controller path.
- Port to `.m` at the end. **Parity test:** log Python inputs/outputs for every scenario, replay them through the `.m` files via `matlab -batch`, assert outputs match to ~1e-9. MATLAB R2025a is installed, so this is automatic and no manual MATLAB work is needed.

### Step sizes (for README)
- **Plant integration:** forward-Euler on `turret_plant` is **unstable for dt > 3.34 ms** (plant poles −5.23 ± 55.7j). README: integrate the plant at **dt ≤ 1 ms** (or ode45).
- **Controller:** the controller computes its own dt from `t`, so any rate works. The real constraint is a **minimum** rate: it must be called at ≥ 20 Hz (radar rate) or it drops measurements. A pure-P loop with K_v = 2 was checked under ZOH: stable at 10, 20, 40, 50, 100 and 1000 Hz. Recommendation: plant at 1 kHz, controller at 1 kHz or any integer multiple of 20 Hz.

### Voltage limit (assumed, stated in paper)
- Tracking need is set by back-EMF: V ≈ 1000·θ̇. Worst LOS rate 0.1875 rad/s and LOS accel 0.060 rad/s² need ≈ 188 V.
- Proposal: **|V_m| ≤ 600 V**, i.e. a rectified 440 VAC 3-phase naval/industrial bus (1.35×440 ≈ 594 V DC). That gives 3.2× headroom over worst-case tracking and ≈ 34°/s slew, which matters for acquisition and interception time. Implement saturation + anti-windup.

### Hidden-scenario stress suite (real-world-inspired guesses)
| ID | Name | What it stresses |
|---|---|---|
| D1 | **Stealth break**: blackout at a later onset (e.g. t = 3 s) and the target starts a max-a maneuver *inside* the blackout | Exactly the V-44 story: maneuver when unseen. Coasting + recovery. |
| D2 | **Weave/jink**: sinusoidal ±a_max in x and y, period ~1–2 s | Worst case for a 0.5 s predictor |
| D3 | **Close-in overhead pass**: r → 800 m, high crossing speed, elevation toward +85° | Max LOS rate/accel, elevation limit |
| D4 | **Repeated blackouts** (spec §2 says "intervals") | Repeated coasting |
| D5 | **Combined**: weave + late blackout + close range | Everything at once |
| MC | **Monte Carlo** random admissible trajectories (random smooth bounded accel, random blackout onset ≥1 s, duration ≤0.5 s) | The real defence against unknown cases; report hit rate and interception-time distribution |

### t_p integration: proposal (pending discussion)
- **Design rule:** make the closed turret loop behave like a **t_p-second delay at low frequency**, T(s) ≈ e^{−s t_p} ≈ 1 − s t_p + s² t_p²/2. Feeding the reference θ̂(t + t_p) then gives θ_turret(t) ≈ θ_target(t), with the prediction exactly cancelling the lag.
  - First-order match: **velocity error constant K_v = 1/t_p = 2 s⁻¹.**
  - Second-order match: second error coefficient = t_p²/2 (achievable with lead-lag or a 2-DOF prefilter). Cancels LOS-acceleration lag.
- Check (pure P, K_p ≈ 2003 V/rad gives K_v = 2): crossover ≈ 2 rad/s, PM ≈ 90°, GM ≈ 14.4 dB (limited by the 8.9 Hz resonance). It is quiet and robust, and the prediction "buys" zero lag without needing bandwidth.
- **Fire decision uses the nowcast** (filter estimate at time t versus the measured turret angle), not the t_p prediction, because the laser is instantaneous.
- **Where the T-MPC idea enters:** (1) a bank of distinct *maneuver-class* hypotheses from the a_max envelope, (2) beam side-of-line as the topological invariant for firing, (3) re-identifying the matching class at blackout recovery, (4) Eq. 12-style hysteresis in class selection.

---

## 12. Round-3 decisions (Nov 2)

### Step size
Build first, then measure. README will state a **validated controller-rate range [20 Hz, f_max]**, with f_max found empirically (filter convergence, numerical conditioning, run time per call). Plant integration ≤ 1 ms regardless.

### Power stage add-on (our own assumption layer)
Source: 230 V 50 Hz single-phase household supply (India). Chain:
EMI filter → bridge rectifier (≈325 V peak) → **PFC boost converter to a 650–700 V DC link** → DC-link capacitor bank → **4-quadrant H-bridge** (SiC MOSFET/IGBT, ~20 kHz PWM) → motor, giving **|V_m| ≤ 600 V**. Brake chopper for regenerated energy.

| Operating point | i [A] | V [V] | P [W] |
|---|---|---|---|
| Worst-case tracking (0.1875 rad/s, 0.06 rad/s²) | 0.47 | 188 | 89 |
| Slew cruise 0.6 rad/s | 0.90 | 601 | 541 |
| Slew accel 1 rad/s² at 0.6 rad/s | 4.10 | 604 | 2477 |
| Slew accel 2 rad/s² at 0.6 rad/s | 7.30 | 607 | 4433 |

- A household 230 V / 16 A circuit gives 3.68 kW. If acceleration is capped at ~1 rad/s² (current limit ~4 A), the drive runs off household supply. Peaks above that come from the DC-link capacitors.
- Braking energy at 0.6 rad/s = ½·3200·0.36 = 576 J. Absorbing it with a 650→750 V rise needs C ≈ 8.2 mF; alternatively use a brake chopper resistor.
- The laser's own power supply is out of scope (state this in the assumptions).

### t_p integration: end-to-end LTI cascade design
- **Time-shift property of LTI systems:** a reference advanced by t_p (e^{+s t_p}) through a closed loop T(s) gives turret angle Θ = T(s)·e^{s t_p}·Θ_T. Perfect pointing needs T(s)e^{s t_p} = 1.
- Treat the whole chain (steady-state KF → t_p predictor → closed loop) as one LTI transfer H(s) from true LOS angle to turret angle, in the linearised LOS frame (small-perturbation, as the spec allows). The predictor is a truncated Taylor advance A_n(s) = Σ_{k≤n}(s t_p)^k/k!, built from KF derivative estimates.
- **Design target:** the error series of 1 − H(s) vanishes up to order m, i.e. end-to-end "type-m" pointing.
  - Order 3 is mandatory: position, velocity and acceleration, all bounded by the spec.
  - Order 4 (jerk) is matched or traded against noise: each extra order needs one more estimated derivative, and noise grows roughly as (ω_KF t_p)^k. The spec bounds a but not jerk.
- The choice of m comes from an **error-budget plot**: truncation error vs. noise-induced error vs. order. This is analytical paper material.
- Use a 2-DOF structure. Feedback C(s) handles robustness (margins, 8.9 Hz resonance). Prefilter/feed-forward F(s) does the coefficient matching. RHP zeros (Padé-like) are allowed in F because it sits outside the loop.

### Firing policy
- **Spec constraint:** ≥ 1.0 s between pulses; continuous operation is *not* available. CIWS-style bursts (Phalanx: 3000–4500 rds/min, 60/100-round bursts) do not apply. Max ~10 pulses in 10 s. There is no ammunition cost (laser), so the cost of a miss is **1 s of lockout** plus the "engagement accuracy" score.
- **No kill feedback** reaches the controller, so we can't do true shoot-look-shoot. Keep engaging after every recharge for the whole run.
- **Rule (fire-now vs wait, an optimal-stopping problem):** fire when
  1. **nowcast** P_hit(now) ≥ p_min (current estimate vs measured turret angle), AND
  2. **prediction consensus:** maneuver-class hypotheses agree the target is inside R_eng, AND
  3. no higher P_hit is predicted within the next 1 s window, discounted by a **blackout hazard**. Since a stealth blackout can start any time after 1 s, holding a good shot is risky, so lean towards firing early.
- **Ambush mode** (fallback when the classes diverge, e.g. after a blackout or a hard maneuver): put the beam at an angle that **every** class trajectory crosses within the window, between θ_now and min_i θ_i(t+τ) in the direction of LOS motion. That gives a guaranteed crossing in the topological (side-of-beam) sense. Fire one pulse when the nowcast confirms the crossing, then return to tracking.

---

## 13. Round-4 decisions (Nov 2): locked for build

- **Power stage:** written up in `docs/power_stage.md`. 230 V outlet → PFC boost → 700 V DC link → H-bridge, |V_m| ≤ 600 V, ≈2.5 kW design point, 11.5 A input (< 16 A outlet), 2 mF link + brake chopper.
- **Error order: go to 4 (jerk).**
  - Reasoning: A/B/C exercise position, velocity and acceleration, so the hidden case plausibly exercises jerk: sudden maneuver onsets, and the post-blackout "jump" where the target appears to have switched acceleration while unseen.
  - Estimator therefore needs a jerk state. Candidate: the **Mehrotra & Mahapatra (1997) jerk model**, an exponentially correlated jerk, the Singer idea pushed one derivative up.
  - The order-vs-error budget plot (truncation error vs noise amplification ~(ω_KF·t_p)^k) is still produced, to justify order 4 numerically and show the noise cost.
- **Firing look-ahead window = 1 s only** (one recharge period). Looking further adds nothing because a fire decision only blocks the next 1 s.
- **Firing judged from radar data only** (no kill feedback). Keep engaging after every recharge.
- **Multiple V-44s:** the interface carries a single (r, θ) track, but the story says swarms. Robustness consequence: a sudden **track discontinuity** (radar locks onto another V-44) must be detected via an innovation gate (χ² test on the innovation, beyond blackout-recovery tolerance). Then re-initialise the filter instead of dragging the old track across. Keep firing logic per-track.
- **Ambush mode** confirmed as the fallback when maneuver classes disagree.

---

## 14. Build log and evidence (Nov 2–3)

### Final architecture (Python `python/guardian/`, MATLAB `matlab/`, parity verified to 4e-9 V)
1. **Estimator:** converted-measurement KF (unbiased conversion), order-4 Markov-jerk model (σ_j = 20 m/s³, τ_j = 0.4 s), envelope projection, χ² track-jump re-init, maneuver-class re-identification after loss.
2. **Predictor:** exact LOS derivatives θ̂…θ̂''' via complex log-derivatives (z'/z chain); r = A₃(D)θ̂ = θ̂(t+t_p) Taylor.
3. **Reference generator:** (A₃ + K)(D) y_d = r + K(D) θ̂. The A₃ lag cancels the t_p lead, so y_d → θ̂(t) exactly. K = d₃(s+p)³ − A₃ places the deviation poles at a triple −20 rad/s. Realised as a cascade (gains p/3, p, 3p) with rate/accel saturation (0.58 rad/s, 2.5 rad/s²), active only during acquisition.
4. **Control:** exact plant-inverse feed-forward V = 0.32 y''' + 3.35 y'' + 1001.5 y', plus PI + roll-off C(s) = 4000(1 + 1/s)/(0.03s + 1). ω_c = 4.1 rad/s, PM 68°, GM 16 dB; the same within 1° and 2 dB at 20 Hz ZOH.
5. **Fire control:** Gaussian mixture over 17 maneuver classes. Tracking: fire at P ≥ 0.999, or arm at P ≥ 0.90 and fire at the peak (−0.02) or after 0.3 s. Coasting: fire only while P ≥ 0.95 and decaying.

### Analytical results for the paper
- Taylor sections of e^{s t_p} are Hurwitz only for degree ≤ 4. Their minimum damping falls 1.0 → 0.71 → 0.36 → 0.11 for degree 1–4. A predictor-matched lag model therefore exists only up to order 4, and the bare A₃ lag rings (ζ = 0.36). This is why deviation feedback K is needed.
- In the LTI limit, prediction + matched lag ≡ nowcast tracking (time-shift property). t_p sets the reference model, not the end-to-end transfer.

### Decisions backed by Monte Carlo (300–500 random admissible engagements, 100 Hz)
| Change | Evidence | Decision |
|---|---|---|
| Bare 1/A₃ reference | 5.3 m pointing RMS (46 m acquisition overshoot, ζ = 0.36) | Rejected → deviation feedback |
| State-clamped reference at p ≥ 8 | 14–90 A current spikes, chatter | Rejected → cascaded saturation |
| p_dev 5 → 20 | D6 RMS 1.59 → 0.61 m, first hit 2.1 → 1.0 s | Accepted (p = 30: noisier, no gain) |
| Fire at first P ≥ 0.9 | Edge shots missed | Rejected → fire-now-or-wait rule |
| Order 2 (CV) | 99.3 % success, 77 % time in R_eng | Rejected |
| Order 3 vs order 4 | Order 3 ~10 % better on smooth MC; order 4 better on jerk/blackout cases (D4, D6, B), and order 3 missed in D4 at 20 Hz | **Order 4** (robust worst case) |
| σ_j 40 → 20 | first-hit max 2.86 → 2.05 s | Accepted |
| Class re-identification | No measurable MC difference vs plain KF recovery | Kept (harmless, motivated), reported honestly |
| Ambush hold | Pointing p95 0.8 → 4.1 m, no hit gain | **Rejected** (kept as option, off) |
| NIS-inflated firing σ | No measurable change | **Rejected** (simplicity) |
| Coast threshold p_min 0.75 → 0.95 | All 9 MC misses were coast shots 0.2–0.43 s into a loss; 9 → 1 misses, no time cost | Accepted |
| Instant-fire p_sure 0.99 → 0.999 | Removes the D6 jerk-reversal miss (fired at P = 0.9901, missed by 2 cm); MC first-hit mean +0.04 s | Accepted |

---

## Sources
- Spec: `Assignment_2_ME_319_Control_Systems_S3.pdf`
- de Groot, Ferranti, Gavrila, Alonso-Mora, *Topology-Driven Parallel Trajectory Optimization in Dynamic Environments*, arXiv:2401.06021v2 (2024): `2401.06021v2.pdf`
- Lerro & Bar-Shalom, *Tracking with debiased consistent converted measurements versus EKF*, IEEE TAES 29(3), 1993: https://data.zhenzhi.com/details/1/988196848
- Bordonaro, Willett, Bar-Shalom, converted-measurement vs EKF results: https://www.bohrium.com/paper-details/consistent-linear-tracker-with-converted-range-bearing-and-range-rate-measurements/813116169740353540-2511
- Polar→Cartesian conversion bias illustration: https://code.eng.buffalo.edu/tracking/unscented_filter/polar2cart.html
- Mazor, Averbuch, Bar-Shalom, Dayan, *IMM methods in target tracking: a survey*, IEEE TAES 34(1), 1998: https://www.bohrium.com/paper-details/interacting-multiple-model-methods-in-target-tracking-a-survey/812459701454241793-105711
- Li & Jilkov, *Overview of multiple-model methods for maneuvering target tracking*, SPIE 2004: https://biomedicaloptics.spiedigitallibrary.org/conference-proceedings-of-spie/5204/0000/Overview-of-multiple-model-methods-for-maneuvering-target-tracking/10.1117/12.511164.full
- Pitre, Jilkov, Li, *Comparative study of MM algorithms*, SPIE 2005: https://electronicimaging.spiedigitallibrary.org/conference-proceedings-of-spie/5809/0000/A-comparative-study-of-multiple-model-algorithms-for-maneuvering-target/10.1117/12.609681.full
- *A Survey of Manoeuvring Target Tracking Methods* (Singer, CSM, coasting): https://arxiv.org/pdf/1503.07828
- MathWorks Singer acceleration model: https://www.mathworks.com/help/fusion/ref/fusion.tracker.transition.singeraccelerationmodel.html
- Current Statistical Model summary: https://rfessentials.com/resources/rf-glossary/current-statistical-model/
- DTIC anti-aircraft gun fire-control (Kalman target-state estimator): https://apps.dtic.mil/sti/pdfs/ADA063704.pdf
- Non-linear prediction for gun fire control (patent US4794235): https://patents.google.com/patent/US4794235
- IMM for warship gun fire control (CV/Singer/CTR): https://koreascience.kr/article/JAKO200706942461766.pdf
- NPS HEL counter-drone tracking chain: https://www.eurekalert.org/news-releases/1074343
- Missing measurements + VSMM-UKF (Lv, Cai, Dong 2022): https://doaj.org/article/0be728b911d647b58687d9424ed3236f
- Williams et al., *Information-theoretic MPC / MPPI* (2017), for the sampling idea (to be cited when we write)
- Phalanx rate of fire / burst length: https://navweaps.com/Weapons/WNUS_Phalanx.php ; https://twz.com/sea/phalanx-ciws-costs-3500-per-second-in-ammo-to-fire
- Laser weapon dwell / magazine limits: https://www.laserwars.net/p/laser-weapon-infinite-magazine-myth ; https://eos-aus.com/defence/high-energy-laser-weapon
- Shoot-look-shoot vs salvo: Glazebrook & Washburn, https://apps.dtic.mil/sti/tr/pdf/ADA487449.pdf ; Pryluk, Shima, Golan, *Shoot-Shoot-Look for an Air Defense System*, https://cris.iucc.ac.il/en/publications/shoot-shoot-look-for-an-air-defense-system/
