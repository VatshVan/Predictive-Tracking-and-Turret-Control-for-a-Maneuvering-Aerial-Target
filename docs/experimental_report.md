# Guardian Fire-Control System: Experimental Report

ME 319 Assignment 2 (Mumbai 2044). This report documents every experiment run on the final design: how test engagements were generated, how performance was measured, which statistical tests were applied, and what the evidence supports and does not support. It is the evidence base for the paper.

All numbers are reproducible from `python/` (commands in §9). Figures are in `docs/figures/`.

---

## 1. System under test

| Block | Final configuration |
|---|---|
| Estimator | Converted-measurement Kalman filter (unbiased polar→Cartesian conversion), Markov-jerk model per axis (position, velocity, acceleration, jerk), σ_j = 20 m/s³, τ_j = 0.4 s; envelope projection (\|v\| ≤ 157.5 m/s, \|a\| ≤ 23 m/s²); χ² track-jump gate 40; maneuver-class re-identification after measurement loss |
| Predictor | Exact LOS derivatives θ̂, θ̂′, θ̂″, θ̂‴ from the Cartesian state; Taylor advance r = A₃(D)θ̂, A₃(s) = 1 + t_p s + (t_p s)²/2 + (t_p s)³/6, t_p = 0.5 s |
| Reference generator | (A₃ + K)(D) y_d = r + K(D)θ̂ with A₃ + K = (t_p³/6)(s + 20)³; cascaded realisation with rate/acceleration saturation 0.58 rad/s, 2.5 rad/s² |
| Feedback | V = G_p⁻¹(D) y_d + C(s)(y_d − θ), C(s) = 4000 (1 + 1/s)/(0.03 s + 1), \|V\| ≤ 600 V with conditional-integration anti-windup |
| Fire control | Gaussian mixture over 17 maneuver classes. Tracking: fire if P_hit ≥ 0.999, or arm at P_hit ≥ 0.90 and fire at the peak (drop of 0.02) or after 0.3 s. Coasting: fire only while P_hit ≥ 0.95 and decaying. Class consensus ≥ 0.80; ≥ 8 radar updates; 1 s recharge |

The plant is the first-principles elevation drive θ/V = 1000 / [s(320 s² + 3350 s + 1.0015·10⁶)], with poles at 0 and −5.23 ± 55.70j (8.90 Hz, ζ = 0.094).

## 2. Simulation protocol

### 2.1 Call pattern (mirrors the evaluation interface)
At each controller step t_k = k / f_ctrl:
1. `[r, θ, valid] = target_motion_generator(t_k, scenario)`
2. `[V, fire] = turret_controller(t_k, r, θ, valid, θ_turret, θ̇_turret)`
3. A laser pulse is **judged at t_k** if `fire` rises 0→1 and at least 1.0 s has passed since the previous pulse. Hit if the minimum distance from the target to the beam ray (origin, direction θ_turret(t_k)) is ≤ R_eng = 2 m.
4. The plant is integrated from t_k to t_{k+1} with V held (zero-order hold), using RK4 at ≤ 1 ms. Mechanical stops at −10° and +85° clamp the state.

Unless stated otherwise: f_ctrl = 100 Hz, the turret starts at rest at 0°, and t_end = 10 s.

### 2.2 Radar model
Measurements exist only at ticks t = 0.05 k. Noise is additive Gaussian: σ_r = 1 m, σ_θ = 0.02°. Noise is indexed by tick number, so a scenario sees the same noise sequence at every controller rate. During blackouts, and between ticks, the generator returns NaN with `meas_valid = 0`.

### 2.3 Metrics

| Metric | Definition |
|---|---|
| Interception success | At least one pulse with miss ≤ 2 m |
| First-hit time | Time of the first successful pulse (the spec's interception time) |
| Pulse hit rate | Hits / pulses issued |
| Pointing RMS | RMS of cross-range pointing error r·sin(θ_turret − θ_true) over t ≥ 2 s |
| Time in R_eng | Fraction of t ≥ 2 s with \|pointing error\| ≤ 2 m |
| Estimation error | r·sin(θ̂ − θ_true), the nowcast's cross-range error |
| Peak voltage | max \|V_m\| |

---

## 3. Test engagements

### 3.1 Standardized scenarios (from the specification)

| ID | Start (x, y) [m] | Initial velocity [m/s] | Maneuver | Measurement loss |
|---|---|---|---|---|
| A | (1700, 700) | (−120, 0) | none | none |
| B | (1700, 700) | (−120, 0) | from t = 1 s: a_x = 8 sin(2π·0.15(t−1)), a_y = 14 sin(2π·0.30(t−1)) with a 0.5 s quintic onset | none |
| C | as A | as A | none | [1.0, 1.5) s |

### 3.2 Stress scenarios (our own, modelling plausible hidden cases)

| ID | Name | Construction | What it stresses |
|---|---|---|---|
| D1 | Stealth break | Pull-up (9, 17) m/s² starting at 3.1 s, **inside** a blackout on [3.0, 3.5) s | Maneuver while unobserved; recovery |
| D2 | Weave | a = (13 sin(2πt/1.6), 14 cos(2πt/1.2)), \|a\| ≈ 19 m/s² | Worst case for a 0.5 s predictor |
| D3 | Close-in pass | r → 853 m, elevation 30° → 83° | Peak LOS rate and acceleration, near the +85° limit |
| D4 | Repeated loss | Scenario B with four blackouts (1.0, 3.2, 5.5, 7.5 s) | Repeated coasting |
| D5 | Combined | Weave + blackout [4.0, 4.5) s, start at 1476 m | Everything at once |
| D6 | Jerk test | ±18 m/s² vertical reversals every 1.2 s, 0.15 s transitions (peak jerk 450 m/s³) | Order-4 (jerk) estimation |

All scenarios are checked programmatically against the envelope (800 ≤ r ≤ 2000 m, \|v\| ≤ 150 m/s, \|a\| ≤ 20 m/s²) and rejected if violated.

### 3.3 Monte Carlo engagements (`MC<seed>`)
Each seed produces one random admissible engagement:

| Quantity | Distribution |
|---|---|
| Initial range | U(1500, 2000) m |
| Initial elevation | U(8°, 55°) |
| Initial speed | U(70, 150) m/s |
| Heading | Inbound (π + θ₀) ± U(−40°, 40°) |
| Maneuver onset | U(0.5, 3.0) s, 0.5 s quintic ramp |
| Maneuver shape | Per axis, a sum of 3 sinusoids with frequency U(0.1, 1.0) Hz, phase U(0, 2π), amplitude U(0.2, 1.0) |
| Break maneuver (probability 0.5) | Smooth pulse with start U(1, 7) s, length U(0.5, 2) s, ramp U(0.15, 0.6) s, random direction |
| Acceleration scale | The whole profile is scaled to a peak of U(0.4, 1.0) × 0.97 × a_max |
| Blackouts | One with onset U(1.0, 8.5) s and duration U(0.1, 0.5) s. With probability 0.3, a second one starting ≥ 0.5 s after the first ends (before 9.5 s) |
| Admissibility | Rejection sampling (≤ 500 draws): envelope as above, elevation within [−8°, 83°] |
| Noise | Independent stream, seed + 10007 |

Figure `mc_paths` shows 80 sampled trajectories and 12 acceleration profiles. Figure `mc_features` shows the realised design space over the 1000 baseline seeds:
- 1284 blackout intervals: 716 runs with one, 284 with two;
- peak \|a\| from 7.0 to 19.4 m/s²;
- peak jerk from 17 to 169 m/s³;
- peak LOS rate up to 0.171 rad/s;
- initial slew from 8° to 55°.

**Limitation of this generator.** The random maneuvers are smooth sinusoid sums, so their jerk (≤ 169 m/s³) is far below the stress scenario D6 (450 m/s³). The noise-ensemble experiment (§6) covers that gap.

---

## 4. Statistical methodology

| Question | Method |
|---|---|
| Interception probability | Wilson score 95 % interval on k/n |
| Pulse hit rate | Pooled rate; 95 % CI by **cluster bootstrap over engagements** (10 000 resamples), because pulses within one engagement are correlated. A naive Wilson CI is also stored |
| First-hit time, pointing RMS, time in R_eng | Percentile bootstrap 95 % CIs (10 000 resamples) |
| Variant vs final design | **Paired design**: identical seeds 0–499, so each pair shares trajectory, blackouts and radar noise |
| — first-hit time, pointing RMS | Wilcoxon signed-rank (Pratt zero handling) + bootstrap CI of the mean paired difference |
| — missed pulses | Two-sided sign-flip permutation test on per-run miss-count differences (20 000 permutations) |
| — interception success | Exact McNemar test (binomial on discordant pairs) |
| Multiplicity | Holm–Bonferroni within each metric family (11 comparisons) |
| Censoring | A run with no hit is assigned t_end = 10 s for paired first-hit comparisons (it never occurred for the final design) |
| Estimator consistency | Normalised error z = (estimation error)/σ̂, σ̂ the filter's cross-range standard deviation; E[z²] and the fraction with \|z\| < 2 compared with N(0, 1) |
| P_hit calibration | Reliability diagram over 270 000 controller samples, Wilson CI per bin, Brier score |
| Difficulty dependence | Spearman rank correlation |

All random draws in the analysis use fixed seeds (bootstrap: 2044; permutation: 7).

---

## 5. Results: final design

### 5.1 Standardized and stress scenarios (single realisation, 1 kHz)

| Scenario | First hit [s] | Pulses hit | Pointing RMS [m] | Time in R_eng |
|---|---|---|---|---|
| A | 1.00 | 9/9 | 0.32 | 100 % |
| B | 1.01 | 9/9 | 0.53 | 100 % |
| C | 0.98 | 10/10 | 0.41 | 100 % |
| D1 | 0.72 | 10/10 | 0.40 | 100 % |
| D2 | 1.03 | 9/9 | 0.99 | 96.3 % |
| D3 | 1.28 | 9/9 | 0.32 | 100 % |
| D4 | 0.98 | 9/9 | 0.59 | 100 % |
| D5 | 1.17 | 9/9 | 0.69 | 100 % |
| D6 | 1.10 | 9/9 | 0.78 | 98.1 % |

Time histories: `scenario_<ID>` figures. Peak voltage 583 V in every case.

### 5.2 Monte Carlo, 1000 engagements (100 Hz)

| Metric | Value | 95 % CI |
|---|---|---|
| Interceptions | 1000 / 1000 | [99.62 %, 100 %] (Wilson); failure probability < 0.38 % |
| Pulse hit rate | 99.96 % (4 misses / 9138) | [99.91 %, 99.99 %] (cluster bootstrap) |
| First-hit mean | 1.260 s | [1.235, 1.286] s |
| First-hit p95 / max | 1.90 s / 2.30 s | |
| Pointing RMS median | 0.559 m | |
| Time in R_eng | 99.40 % | |
| Peak voltage | 589 V | (limit 600 V) |

**Miss distances** (`mc_miss_hist`): median 0.37 m, p90 0.90 m, p99 1.45 m, max 2.26 m. All four misses were tracking-mode shots that missed by 0.01–0.26 m.

**Interception time is acquisition-limited** (`mc_first_hit`):
- First-hit time correlates with the initial slew \|θ_T(0) − θ(0)\| at Spearman ρ = 0.94.
- It is essentially uncorrelated with maneuver intensity (ρ = 0.05 with RMS acceleration), peak jerk (ρ = 0.03) and blackout total (ρ = 0.06, p = 0.05).
- Engagements needing < 10° of slew are intercepted at 0.52 s on average.
- The scatter lies on a line parallel to the pure-slew bound θ/0.58 rad/s, offset by settling and the 8-update firing gate.
- In short, the time to the first shot is set by turret kinematics, not by the estimator or the target.

**Pointing accuracy degrades gracefully with difficulty** (`mc_difficulty`). Pointing RMS correlates with RMS target acceleration (ρ = 0.45, p ≈ 10⁻⁵¹), peak acceleration (ρ = 0.41), peak jerk (ρ = 0.30) and total blackout time (ρ = 0.22). The 95th-percentile RMS stays at 0.83 m.

### 5.3 Pointing error over time (`mc_error_bands`, 300 logged engagements)
- After acquisition (t > 2 s), the median \|pointing error\| stays around 0.3–0.4 m.
- The 90th percentile stays below 1 m.
- The 99th percentile touches 2 m only at blackout recoveries.

### 5.4 Behaviour across measurement loss (`mc_blackout_aligned`, 276 aligned blackouts)
- **Before onset:** median pointing error ≈ 0.4 m, estimation error ≈ 0.25 m, filter σ ≈ 0.35 m.
- **During the loss:** the filter σ grows (median ≈ 1.0 m at 0.3 s). Estimation error grows more slowly (median ≈ 0.5 m), and pointing error lags it (median peak ≈ 0.8 m at ≈ 0.35 s).
- **Worst runs:** the 90th-percentile band reaches ≈ 2 m at the end of the longest losses. This is consistent with the ½·a_max·t² = 2.5 m worst-case drift at 0.5 s.
- **On recovery:** the first measurement collapses σ immediately. Pointing returns to its pre-loss level within ≈ 0.4 s, set by the reference generator's triple pole at −20 rad/s.

### 5.5 Estimator consistency (`mc_consistency`)

| Mode | E[z²] (ideal 1.0) | \|z\| < 2 (ideal 95.4 %) |
|---|---|---|
| Tracking (t ≥ 2 s) | **1.008** | 95.2 % |
| Coasting (measurement loss) | 0.919 | 96.5 % |

The filter is statistically consistent while tracking, and slightly conservative while coasting. That conservatism is the intended behaviour, because coasting shots are the risky ones.

### 5.6 Hit-probability calibration (`mc_calibration`)
- Over about 270 000 controller samples, observed frequency tracks predicted P_hit across all bins (Brier score 0.0045).
- Mid-range bins are conservative: predicted 0.62 → observed 0.69; predicted 0.75 → observed 0.83.
- At the 9138 fire instants the mean predicted P_hit was 0.9991, implying 7.8 expected misses. Only 4 occurred.
- So the firing decision is an honest probability statement, and errs on the safe side.

---

## 6. Results: design decisions (paired ablations)

### 6.1 Monte Carlo ablations (500 paired seeds, Holm-adjusted)

| Variant (one change from final) | Δ first hit [s] (95 % CI), p | Δ pointing RMS [m] (95 % CI), p | Misses final → variant, p | Verdict |
|---|---|---|---|---|
| Order 2 (constant velocity) | +0.172 (0.128, 0.226), 2·10⁻⁴⁶ | +0.878 (0.814, 0.939), 10⁻⁷⁴ | 2 → 77, 5.5·10⁻⁴; **2 interceptions lost** | Strongly worse |
| Deviation pole p = 10 | +0.418 (0.401, 0.435), 10⁻⁸¹ | +0.206 (0.191, 0.222), 10⁻⁸² | 2 → 13, 0.037 | Strongly worse |
| Ambush hold enabled | +0.001 (−0.001, 0.003), 1 | +0.822 (0.717, 0.932), 6·10⁻⁵¹ | 2 → 1, 1 | Worse pointing, no hit gain: **rejected** |
| Controller at 20 Hz | +0.369 (0.346, 0.393), 6·10⁻⁷⁸ | +0.369 (0.348, 0.391), 10⁻⁸² | 2 → 4, 1 | Slower, still 500/500 |
| Jerk noise σ_j = 40 | +0.021 (0.011, 0.031), 3·10⁻⁷ | +0.022 (0.015, 0.029), 3·10⁻¹⁰ | 2 → 2, 1 | Slightly worse |
| Deviation pole p = 30 | +0.179 (0.159, 0.200), 2·10⁻³² | −0.021 (−0.037, −0.001), 3·10⁻⁴⁹ | 2 → 1, 1 | Slower acquisition for marginal pointing gain |
| Instant-fire threshold 0.99 | **−0.047** (−0.058, −0.037), 2·10⁻³⁶ | 0 | 2 → 2, 1 | Faster on smooth targets (see §6.2) |
| Coast threshold 0.75 | −0.002 (−0.004, 0.000), 0.63 | 0 | 2 → 8, 0.62 | More misses; not significant at this n |
| No class re-identification | 0.000, 1 | 0.000, 0.77 | 2 → 2, 1 | **No measurable effect** |
| Order 3 (Singer acceleration) | −0.016 (−0.026, −0.005), 0.66 | **−0.061** (−0.066, −0.056), 5·10⁻⁷⁰ | 2 → 0, 1 | Better on smooth targets (see §6.2) |
| Turret starts at 60° | −0.149 (−0.216, −0.080), 1.8·10⁻⁴ | −0.008, 0.002 | 2 → 1, 1 | Shorter slew → earlier hit |

Figure: `mc_ablation_forest`.

### 6.2 Stress-scenario noise ensembles

<<STRESS>>

### 6.3 What the evidence supports

<<VERDICT>>

---

## 7. Controller rate, timing and implementation parity

**Rate study** (24 cases per rate: A–D6 and 15 Monte Carlo seeds; `rate_study`):

| f_ctrl [Hz] | Intercepted | Pulses hit | First hit mean / max [s] | Median pointing RMS [m] |
|---|---|---|---|---|
| 20 | 24/24 | 200/200 | 1.70 / 2.65 | 0.76 |
| 40 | 24/24 | 210/211 | 1.49 / 2.25 | 0.63 |
| 60 | 24/24 | 212/212 | 1.34 / 2.10 | 0.56 |
| 100 | 24/24 | 217/217 | 1.26 / 2.05 | 0.54 |
| 200 | 24/24 | 220/220 | 1.21 / 1.90 | 0.53 |
| 500 | 24/24 | 221/221 | 1.22 / 1.91 | 0.52 |
| 1000 | 24/24 | 221/221 | 1.23 / 1.91 | 0.52 |
| 2000 | 24/24 | 221/221 | 1.23 / 1.91 | 0.52 |
| 5000 | 24/24 | 221/221 | 1.23 / 1.92 | 0.52 |

- Performance saturates at about 200 Hz.
- Forward-Euler plant integration at 1 ms (controller 100 Hz) works: Scenario B first hit at 0.95 s. Forward Euler is unstable above 3.34 ms.

**Feedback margins** (exact discrete C(z) with ZOH plant):

| Rate | ω_c [rad/s] | PM [°] | GM [dB] |
|---|---|---|---|
| continuous | 4.10 | 68.5 | 16.1 |
| 100 Hz | 4.12 | 68.5 | 16.1 |
| 50 Hz | 4.14 | 68.4 | 16.3 |
| 20 Hz | 4.19 | 67.4 | 18.6 |

**Execution time** (scenario B, 1 kHz, 10 001 calls, quiet machine): `turret_controller.m` mean 0.51 ms, median 0.33 ms; Python mean 0.26 ms. Real-time use is practical up to about 1 kHz.

**Parity.** 12 logged runs (A, B, C, D1, D4, D6 at 20 Hz and 1 kHz) were replayed open-loop through `turret_controller.m`. The maximum \|ΔV\| is 3·10⁻⁹ V, with zero fire-command differences.

---

## 8. Threats to validity

1. **Generator realism.** The Monte Carlo maneuvers are smooth sinusoid sums with a 50 % chance of a break pulse. Their jerk is bounded well below a bang-bang adversary. The stress ensembles test sharper cases, but no generator can stand in for the instructor's hidden cases.
2. **Single noise model.** Radar noise is white Gaussian as specified. Glint, bias and range-dependent noise are not modelled.
3. **Plant exactness.** The feed-forward inverts the same plant the simulator integrates. Robustness to parameter error rests on the feedback margins (PM 68°, GM 16 dB), not on simulation.
4. **Laser model.** Pulses are instantaneous rays, per the spec. Dwell and beam divergence are not modelled.
5. **Small-effect power.** With 2 baseline misses in 500 paired runs, differences of a few misses (e.g. coast threshold 0.75: 8 misses) cannot reach significance. Miss-count conclusions are directional only.
6. **Rate-study scale.** 24 cases per rate support "no failures", not fine-grained comparisons between neighbouring rates.

---

## 9. Reproducing every number

```bash
cd python
python -m pytest -q tests                     # 40 tests
python scripts/design_analysis.py             # plant, margins, Taylor-polynomial roots
python scripts/plot_scenarios.py              # §5.1 table + scenario figures
bash scripts/run_ablations.sh                 # §5.2, §6.1 Monte Carlo campaign
python scripts/mc_ensemble.py --n 300         # §5.3-5.6 logged ensemble
python scripts/statistics.py                  # CIs and paired tests -> results/stats_report.md
python scripts/stress_ensemble.py --reps 50   # §6.2 -> results/stress_report.md
python scripts/rate_study.py                  # §7
python scripts/plot_experiments.py            # Monte Carlo figures
python scripts/export_parity.py               # then in MATLAB: parity_check, time_controller
```
