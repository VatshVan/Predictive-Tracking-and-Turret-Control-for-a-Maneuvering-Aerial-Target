# Guardian: Predictive Tracking and Turret Control for a Maneuvering Aerial Target

ME 319 Control Systems (Autumn 2026, IIT Bombay), Assignment 2: *Mumbai 2044*.

A predictive-feedback fire-control system for the elevation axis of a laser turret. It estimates a maneuvering V-44 target from noisy 20 Hz radar (range, elevation), predicts its elevation t_p = 0.5 s ahead, steers the turret with an LTI two-degree-of-freedom controller, and decides when to fire the 1 Hz pulsed laser. It keeps working through radar blackouts.

## Architecture

```
radar (r, θ) ─► converted-measurement KF ─► LOS nowcast θ̂, θ̂', θ̂'', θ̂''' ─► Taylor advance A₃(s): r = θ̂(t+t_p)
                (Markov-jerk model,           │                                 │
                 envelope projection,         │                                 ▼
                 class re-identification)     │                  reference generator (A₃+K)(D) y_d = r + K(D) θ̂
                                              │                                 │ y_d, y_d', y_d'', y_d'''
                                              ▼                                 ▼
                         maneuver-class hit probability          V = Gp⁻¹(D) y_d  +  C(s)(y_d − θ)  ─► turret
                         + consensus + fire-now-or-wait ─► fire_cmd       (feed-forward)   (PI + roll-off)
```

The core identity: in the linear time-invariant limit the reference generator's lag `1/A₃(s)` exactly cancels the predictor's lead `A₃(s) ≈ e^{s t_p}`, so the turret tracks the target's *current* elevation (the laser is instantaneous) while the controller only ever sees the t_p-predicted reference. The closed loop is therefore invariant to t_p, and the predicted reference adds no information beyond the estimated derivatives; the useful prediction is the filter's motion model (coasting through blackouts, feeding the plant inverse and the firing look-ahead). The deviation feedback `K(s)` places the transient poles (triple pole at −20 rad/s) without disturbing that cancellation.

## Repository layout

| Path | Contents |
|---|---|
| `python/guardian/turret_controller.py` | Estimator, predictor, reference generator, feedback, fire control |
| `python/guardian/turret_plant.py` | First-principles elevation drive, ±600 V supply limit |
| `python/guardian/target_motion_generator.py` | Scenarios A, B, C, stress scenarios D1–D6, Monte Carlo family `MC<seed>` |
| `python/guardian/sim.py` | Closed-loop harness mirroring the evaluation call pattern, pulse judging |
| `python/scripts/` | Experiments (`run_all_experiments.sh` runs all), statistics (`mc_statistics.py`), figures, paper-number export |
| `python/tests/` | pytest suite |
| `matlab/` | MATLAB port: `turret_controller.m`, `turret_plant.m`, `target_motion_generator.m`, run scripts |
| `docs/power_stage.md` | 230 V outlet → ±600 V drive design (the voltage-limit assumption) |
| `docs/Guardian_Experimental_Report.pdf` | Detailed experimental report (36 pages, all figures, statistics, final-evaluation analysis); source in `docs/report/` |
| `docs/research_lemma_homotopy.md` | Order-4 cap (contact order), alternative reference laws, homotopy-guided adversarial testing |
| `python/guardian/homotopy.py` | Beam-relative homotopy classes and the `HG<seed>` adversarial scenarios |
| `paper/` | IEEE paper source (`main.tex`; numbers generated into `numbers.tex`) |
| `research_notes.md` | Design research and decisions |
| `papers/README.md` | References (PDFs not redistributed) |

## Reproducing the results (Python)

```bash
cd python
pip install -r requirements.txt
python -m pytest -q tests
python scripts/plot_scenarios.py
python scripts/design_analysis.py
python scripts/monte_carlo.py --n 300 --f-ctrl 100 --tag baseline
python scripts/rate_study.py
```

Figures and Monte Carlo outputs are written to `python/results/`.

## Step size and call rate (for the evaluation script)

| Item | Requirement |
|---|---|
| Plant integration step | **≤ 1 ms** (or `ode45`). Forward Euler on `turret_plant` is unstable for steps above 3.34 ms because of the 8.9 Hz, ζ = 0.094 drive resonance. |
| Controller call rate | **Recommended 100 Hz to 1 kHz; validated from 20 Hz to 5 kHz.** Use an integer multiple of the 20 Hz radar rate so every radar sample lands on a call. The controller measures its own time step from `t`, so no rate needs configuring. Below 20 Hz radar samples are lost. `turret_controller.m` costs about 0.5 ms per call (MATLAB R2025a, median 0.33 ms), which caps real-time use near 1 kHz. |
| Reset between runs | Automatic: the controller re-initialises when `t` decreases or `t == 0`. In MATLAB, `clear turret_controller` also works. |
| Initial turret state | Any; the controller starts from the turret state at its first call. |

## Running in MATLAB

```matlab
cd matlab
run_scenarios          % A, B, C, D1-D6 at 100 Hz: metrics table + figures in matlab/figures/
o = run_engagement('B', 1000);   % one engagement at a chosen controller rate
time_controller        % per-call execution time
```

`parity_check` replays logged Python runs (create them with `python scripts/export_parity.py`) through `turret_controller.m`. The two implementations agree to within 4e-9 V, with identical fire commands.

## Results

Standardized and stress scenarios (1 kHz controller, turret starting at 0°):

| Scenario | First hit [s] | Pulses hit | Pointing RMS after 2 s [m] |
|---|---|---|---|
| A nominal | 1.00 | 9/9 | 0.32 |
| B maneuvering | 1.01 | 9/9 | 0.53 |
| C measurement loss | 0.98 | 10/10 | 0.41 |
| D1 stealth break | 0.72 | 10/10 | 0.40 |
| D2 weave | 1.03 | 9/9 | 0.99 |
| D3 close-in pass | 1.28 | 9/9 | 0.32 |
| D4 repeated loss | 0.98 | 9/9 | 0.59 |
| D5 combined | 1.17 | 9/9 | 0.69 |
| D6 jerk test | 1.10 | 8/9 | 0.78 |

The single D6 miss (2 cm outside R_eng, right after a ±18 m/s² reversal) belongs to this noise realisation only: across 50 noise realisations of every scenario, the design intercepted 450/450 with 99.90 % of pulses hitting.

Monte Carlo over random admissible engagements (smooth bounded maneuvers, random loss onset ≥ 1 s and duration ≤ 0.5 s), with 95 % confidence intervals:

| Set | Intercepted | Pulse hit rate | First hit mean [CI] / p95 / max [s] |
|---|---|---|---|
| 1000 runs, 100 Hz | 1000/1000 (failure < 0.38 %) | 99.96 % [99.91, 99.99] | 1.21 [1.19, 1.24] / 1.83 / 2.14 |
| 500 runs, 20 Hz | 500/500 | 99.81 % | 1.62 / 2.20 / 3.25 |
| 500 runs, 100 Hz, turret starting at 60° | 500/500 | 99.94 % | 1.07 / 1.70 / 2.64 |

Only the 1000 independent runs enter the failure bound. The 20 Hz and 60° sets reuse the same seeds.

Beyond the random engagements:
- **Shot-timed adversary.** On 1000 white-box engagements where the target reverses its maneuver just before each of the controller's shots (`HG<seed>`), every engagement was intercepted, with 97.3 % of pulses hitting. The hit probability is overconfident there (99.7 % predicted).
- **Plant error.** With R, L, J_eq and b each off by up to ±20 %, the design intercepted 1000/1000 with unchanged accuracy. With K_t and K_e also perturbed it still intercepted 1000/1000 (99.82 % of pulses hit), but median pointing RMS rose to 1.69 m: the feed-forward's velocity term is ≈ N·K_e. Friction is not modelled.
- **Alternative reference laws.** Laws that keep the Taylor prediction inside the loop (exact Taylor n = 3, 4, moment-matched, contact order 1–2) miss less per shot under attack. They point 0.4–4.4 m worse (still 1.6–3.2× worse with no initial slew) and deliver fewer hits, so the deviation-feedback design stays (`docs/research_lemma_homotopy.md`).

Feedback loop: ω_c = 4.1 rad/s, phase margin 68°, gain margin 16 dB (18.6 dB at 20 Hz ZOH).

The full experimental record (generation of test paths, statistical methods, paired ablations, stress ensembles, estimator consistency and hit-probability calibration) is in [docs/Guardian_Experimental_Report.pdf](docs/Guardian_Experimental_Report.pdf). `bash python/scripts/run_all_experiments.sh` regenerates every number and figure.

## Assumptions beyond the specification

- Motor voltage limited to **±600 V** (rectified-and-boosted 230 V outlet; see `docs/power_stage.md`).
- `meas_valid = 1` only when a new radar sample exists; repeated identical samples are ignored.
- Mechanical stops at −10° and +85° clamp the turret state in the simulation harness.
