# Contact order, pole placement and beam-relative homotopy

Research extension of the Guardian controller. Every mode below sits behind a `CFG` flag in
`python/guardian/turret_controller.py`. The defaults are unchanged and reproduce the final design
bit for bit (checked against `results/mc_S_final.json`). No CFG default, `matlab/`, `paper/` or
`docs/report/` file was changed.

Raw numbers are in `python/results/research/` (regenerate with the scripts in §8):
`contact_order.md`, `pole_shift.md`, `research_report.md`, `deconfound/deconfound_report.md`.

## 1. Why this study

A self-review of the final design found two gaps:

1. **The prediction is cancelled.** The reference generator solves
   D(D) y_d = r − K(D)(y_d − θ̂) with r = A₃(D) θ̂ and A₃ + K = D_dev = (t_p³/6)(s + 20)³. That
   polynomial agrees with e^{s t_p} only at s = 0 (contact order 0). The order-4 cap on the
   Taylor advance therefore never binds in the loop.
2. **"Homotopy" was nominal.** The hit probability used 17 hand-picked constant-acceleration
   offsets stacked on a σ that already contains the maneuver noise.

## 2. Contact order (`scripts/research_contact_order.py`)

Let Q(x) = Σ q_k x^k with q_k = 1/k! for k ≤ m, so Q − e^x = O(x^{m+1}) (contact order m), plus f
free higher coefficients.

**Fixed minors.** The reciprocal polynomial's Hurwitz minors Δ_j with 2j − 1 ≤ m are fixed by
the contact order:

| Δ₁ | Δ₂ | Δ₃ | Δ₄ | Δ₅ | Δ₆ |
|---|---|---|---|---|---|
| 1 | 1/3 | 1/45 | 1/4725 | 1/4465125 | 1/46414974375 |

All are positive, so the necessary conditions never rule out a stable completion. They do
shrink super-exponentially.

**Best minimum root damping ζ_min** (differential evolution, 3 restarts):

| m | f=0 (Taylor) | f=1 | f=2 | f=3 | f=4 | f=5 | f=6 |
|---|---|---|---|---|---|---|---|
| 3 | +0.362 | +0.362 | +0.436 | +0.436 | +0.436 | +0.436 | +0.436 |
| 4 | +0.107 | +0.107 | +0.222 | +0.222 | +0.222 | +0.222 | +0.271 |
| 5 | −0.076 | −0.076 | +0.056 | +0.056 | +0.136 | +0.136 | +0.136 |
| 6 | −0.212 | −0.212 | −0.074 | −0.074 | −0.073 | +0.014 | +0.014 |
| 7 | −0.316 | −0.316 | −0.177 | −0.177 | −0.086 | −0.086 | −0.086 |
| 8 | −0.397 | −0.397 | −0.261 | −0.261 | −0.170 | −0.170 | −0.170 |

- Pure Taylor is Hurwitz exactly for n ≤ 4.
- Free coefficients help in pairs; an odd extra term is driven to zero.
- Order 5 needs 2 extra terms. Order 6 needs 5 and is then barely stable.
- Orders 7–8 found no Hurwitz completion with up to 6 extra terms.

**Contact order is effectively capped at 4–5 in a stable loop.**

## 3. The reference generator is a deviation-pole choice

In the linear band, every reference mode in this study has the form P(D) y_d = N(D) θ̂ with
N = P on the estimated derivatives. The deviation therefore obeys

  **P(D)(y_d − θ̂) = 0**, driven only by the estimator's correction step at each 20 Hz radar update.

"Putting the prediction in the loop" is algebraically the choice P = A₃(t_p s), i.e. the
deviation poles become the roots of the Taylor polynomial. The final design is the choice
P = d₃(s + 20)³.

**Decay bound** (`scripts/research_pole_shift.py`). If P keeps first-order prediction fidelity
(p₀ = 1, p₁ = t_p), then Σ 1/s_i = −t_p. A pole with decay rate σ has |Re(1/s_i)| ≤ 1/σ. Hence

  **σ_min ≤ N / t_p**: at most 6 rad/s for a cubic, versus 20 rad/s in the final design.

| fidelity kept | design | slowest pole | ζ |
|---|---|---|---|
| order 3 (exact Taylor) | A₃(t_p s) | 1.40 rad/s | 0.362 |
| order 3, degree 5 | moment-matched Q = A₃ + q₄x⁴ + q₅x⁵ | 2.89 | 0.436 |
| order 3, degree 9 (best found) | 6 free terms | 4.43 | 0.458 |
| **order 2** | **((s² + 8s + 32)/32)², poles −4 ± 4j (double)** | **4.00** | **0.707** |
| **order 1** | **(1 + s/6)³, the bound attained** | **6.00** | **1.000** |
| none (final design) | d₃(s + 20)³ | 20.0 | 1.000 |

The contact-order-2 design has degree 4. Its θ̂'''' input term is dropped, because the order-4
estimator does not estimate it.

**Weighted average** (`ref_mode='blend'`): P_α = (1 − α) A₃ + α d₃(s + 20)³.

| α | 0 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 1 |
|---|---|---|---|---|---|---|---|
| slowest pole [rad/s] | 1.40 | 1.10 | 1.95 | 4.06 | 6.35 | 10.8 | 20 |
| p₁/p₀ (t_p = 0.5 is faithful) | 0.500 | 0.186 | 0.168 | 0.158 | 0.155 | 0.152 | 0.150 |

The blend moves the poles left, but it gives up prediction fidelity almost immediately: at
α = 0.05 the velocity term is already 0.19 instead of 0.5. The contact-order designs above are
the principled way to move the poles left while keeping the prediction terms exact.

## 4. Homotopy, defined the T-MPC way

T-MPC (de Groot et al. 2024) works in state-time space. There, two trajectories are homotopic iff
they can be deformed into each other without crossing an obstacle (Def. 1), and classes are
compared by H-signature or winding numbers (App. A).

**Mapping to the turret** (`H_THROUGH`, `_h_plan` in the controller; `guardian/homotopy.py`):
- **Plane.** Work in the plane (t, e), where e is the target's signed cross-range offset from the
  turret's predicted beam path.
- **Obstacles.** These are the firing events, i.e. the hit tubes O_k = {t_k} × [−R_eng, R_eng].
- **Classes.** A target path is a graph e(t), so its winding number around each O_k is fixed by
  the side on which it passes. Two paths are homotopic iff they pass every event on the same
  side. The class of a path is its H-signature h_k ∈ {+1 above, 0 through = hit, −1 below}.
- **Plans and cost.** The shooter's plans are firing events, and a plan's cost is
  1 − P(the target passes through that event). The previously selected plan, re-identified by
  its absolute time, gets T-MPC's Eq. 12 consistency weight c = 0.75.

**Fire planner** (`phit_mode='homotopy'`):
- **Hypotheses.** 25 maneuver hypotheses start at the last radar update: nominal, plus 8
  directions × {constant, onset at mid-horizon, reversing at mid-horizon}, each at the largest
  admissible |â + Δa| ≤ 20 m/s².
- **Weights.** Posterior weights are a prior times N(Δa; 0, P_a + (10 m/s²)² I).
- **Position spread.** Each hypothesis keeps a Gaussian spread from the noise-free F P Fᵀ
  (Rao-Blackwellised), so maneuver offsets are no longer stacked on σ.
- **Window.** Events are spaced 50 ms apart over the remaining window: one recharge period from
  first arming. The beam path comes from rolling the reference generator forward.

**HG adversary** (closed loop, `HG<seed>`, class = seed % 5):
- **Why closed loop.** Inside the envelope (|a| ≤ 20 m/s², blackouts ≤ 0.5 s), the beam follows
  the radar after a blackout. A target can therefore leave the tube at a firing event only
  through estimator lag.
- **Search.**
  1. Find the final design's firing events on a base engagement.
  2. Reverse the normal-to-LOS acceleration at full envelope (quintic, finite jerk) a lead δ
     before each event, toward the side the class signature requires.
  3. Choose δ ∈ {0.1, 0.2, 0.35, 0.5, 0.7} s to maximize the final design's mean miss.
  4. Re-run once at the fixed point. Envelope violators are rejected.
- **Classes.** The five classes are the signatures (+…+), (−…−), (+…−), (−…+) and (+−+…).
- **Caching.** Cases are cached, so every variant faces identical trajectories.
  `HX<target>_<seed>` re-tunes the same search against another design (§7).

During generation, the realized signature had every pulse "through" in 633/1000 cases. Only 116
realized the intended class, even with full knowledge of the shot timing. Where the final design
does miss, the side follows the attack: "above" seeds miss 38 above vs 4 below, and "below"
seeds miss 1 above vs 31 below.

**What this adversary is and is not.** The homotopy classes do steer the search: the side of
the misses follows the intended class. But the adversary is weak, and most engagements still
end with every pulse "through". The HG family is therefore worst-case evidence about where
misses come from. It is not a Monte Carlo that reliably finds failures, and it does not
estimate a failure probability.

## 5. Results: MC0–999 (random engagements; the failure bound)

1000 independent seeds; the controller runs at 100 Hz and the plant uses RK4 at 1 ms. Every
configuration intercepts **1000/1000**. The final design's failure probability is
**≤ 0.38 %** (one-sided 97.5 % Wilson bound).

| Design | Pointing RMS median | Misses / pulses | Time in R_eng | Δ RMS vs final [95 % CI] | misses final → design (p, Holm) |
|---|---|---|---|---|---|
| **Final design** | **0.559 m** | **4 / 9197** | **99.40 %** | — | — |
| Blend α = 0.5 | 0.579 | 5 / 9203 | 99.33 % | +0.023 [0.022, 0.023] | 4 → 5 (1) |
| Blend α = 0.3 | 0.608 | 6 / 9203 | 99.22 % | +0.053 [0.052, 0.054] | 4 → 6 (1) |
| Blend α = 0.2 | 0.648 | 7 / 9194 | 99.02 % | +0.096 [0.094, 0.097] | 4 → 7 (1) |
| Blend α = 0.1 | 0.817 | 17 / 9125 | 97.23 % | +0.276 [0.269, 0.284] | 4 → 17 (0.07) |
| Contact order 1 | 0.985 | 24 / 8931 | 94.04 % | +0.453 [0.437, 0.469] | 4 → 24 (0.004) |
| Blend α = 0.05 | 1.249 | 23 / 8900 | 89.78 % | +0.788 [0.759, 0.816] | 4 → 23 (0.004) |
| Contact order 2 | 1.821 | 18 / 8461 | 81.02 % | +1.664 [1.593, 1.738] | 4 → 18 (0.04) |
| Exact Taylor n = 3 | 2.817 | 23 / 8056 | 73.33 % | +2.494 [2.415, 2.573] | 4 → 23 (0.004) |
| Moment-matched | 2.941 | 33 / 7937 | 66.01 % | +2.733 [2.644, 2.822] | 4 → 33 (0.0007) |
| Exact Taylor n = 4 | 4.907 | 23 / 7463 | 44.73 % | +4.407 [4.315, 4.499] | 4 → 23 (0.004) |
| Homotopy fire planner | 0.559 | 11 / 9228 | 99.40 % | 0 | 4 → 11 (0.4) |
| No re-identification | 0.558 | 3 / 9197 | 99.40 % | 0 | 4 → 3 (1) |
| Final design, plant ±20 % | 0.557 | 6 / 9190 | 99.40 % | 0 | 4 → 6 (1) |

- **First-hit time.** All blends and contact order 1 hit 3–9 ms earlier (significant). Exact4,
  contact2 and moment-matched hit 30–83 ms later. The homotopy planner hits 19 ms earlier.
- **Homotopy planner activity.** It waited at least once in 18 % of runs.
- **Envelope projection.** It fires in 44 % of runs for every order-4 design. The order-5
  estimator (exact4) is projected about 36 times per run: its acceleration estimate overshoots
  1.15 × 20 m/s².
- **exact4 tuning.** exact4 runs with a retuned estimator, `sigma_top=150, tau_top=0.2` (snap
  model; see `scripts/run_research.sh`). Every other design uses the default order-4 estimator.
  This tuning was not optimized, so the exact4 numbers are pessimistic.

## 6. Results: HG0–999 (closed-loop adversary tuned against the final design)

Worst-case evidence, kept separate from the failure bound. Every configuration intercepts **1000/1000**.

| Design | Misses / pulses | Mean miss | Pointing RMS median | misses final → design (p, Holm) |
|---|---|---|---|---|
| Final design | 244 / 9145 (2.67 %) | ~1.0 m | 1.209 m | — |
| Blend α = 0.2 / 0.3 / 0.5 | 322 / 290 / 262 | ~1.0 m | 1.41 / 1.32 / 1.26 | worse or n.s. |
| Homotopy fire planner | 230 / 9174 | ~1.0 m | 1.209 | 244 → 230 (0.68) |
| No re-identification | 245 / 9146 | ~1.0 m | 1.211 | 244 → 245 (1) |
| Blend α = 0.05 | 187 / 8573 | 0.88 m | 2.485 | 244 → 187 (0.04) |
| Moment-matched | 86 / 7862 | 0.73 m | 4.780 | 244 → 86 (0.0006) |
| **Contact order 1** | **93 / 8608** | **0.68 m** | 1.903 | **244 → 93 (0.0006)** |
| Contact order 2 | 89 / 8205 | 0.70 m | 3.092 | 244 → 89 (0.0006) |
| Exact Taylor n = 3 | 74 / 7876 | 0.70 m | 3.834 | 244 → 74 (0.0006) |
| Exact Taylor n = 4 | 29 / 7488 | 0.82 m | 6.263 | 244 → 29 (0.0006) |

Under this adversary the **per-pulse miss rate** ranking flips: the slow-pole designs miss
60–88 % fewer pulses. The likely mechanism: the 20 rad/s deviation loop chases every estimator
correction, so the lag the adversary induces is converted into beam motion right at the shot.
Slow poles smooth these corrections out.

The miss rate is not the whole picture. The slow-pole designs fire fewer pulses, because their
beam is inside R_eng far less of the time (time in R_eng 26–65 % vs 93 %) and the hit-probability
gate holds the shot back. Counted as **hits delivered**, the final design is still best:

| Design | Hits delivered (pulses − misses) | Time in R_eng |
|---|---|---|
| Final design | **8901** | 93.4 % |
| Homotopy fire planner | 8944 | 93.4 % |
| Contact order 1 | 8515 | 64.7 % |
| Blend α = 0.05 | 8386 | 51.9 % |
| Contact order 2 | 8116 | 46.8 % |
| Exact Taylor n = 3 | 7802 | 55.2 % |
| Moment-matched | 7776 | 29.4 % |
| Exact Taylor n = 4 | 7459 | 26.2 % |

The slow designs look robust because the firing gate refuses doubtful shots, not because they
point better. Time to first hit, the quantity that decides an interception, is the same within
about 50 ms for all of them.

The adversary in this table was tuned against the final design only. §7 checks whether the
advantage survives an adversary tuned against each design itself.

## 7. De-confounding: attacker × defender

The same closed-loop search was re-tuned against contact order 1 and exact Taylor n = 3, on
seeds 0–199. Each seed keeps the base engagement and radar noise of HG<seed>.
`scripts/research_deconfound.py` produces the attacker × defender matrix:

| adversary tuned on ↓ / design → | final | contact order 1 | exact Taylor n = 3 |
|---|---|---|---|
| final | **49 / 1832 (2.67 %)**, 0.96 m | 20 / 1729 (1.16 %), 0.67 m | 14 / 1580 (0.89 %), 0.72 m |
| contact order 1 | 36 / 1829 (1.97 %), 0.84 m | **24 / 1723 (1.39 %)**, 0.67 m | 15 / 1576 (0.95 %), 0.72 m |
| exact Taylor n = 3 | 40 / 1826 (2.19 %), 0.86 m | 20 / 1637 (1.22 %), 0.70 m | **14 / 1502 (0.93 %)**, 0.72 m |

Cells show misses / pulses and the mean miss distance. All 200/200 engagements were intercepted
in every cell. Within a row, cells share trajectories and are paired with the final design:

| adversary | design | misses final → design | p (Holm, 6 tests) |
|---|---|---|---|
| final | contact order 1 | 49 → 20 | 0.004 |
| final | exact Taylor n = 3 | 49 → 14 | 0.0003 |
| contact order 1 | contact order 1 | 36 → 24 | 0.15 |
| contact order 1 | exact Taylor n = 3 | 36 → 15 | 0.008 |
| exact Taylor n = 3 | contact order 1 | 40 → 20 | 0.016 |
| exact Taylor n = 3 | exact Taylor n = 3 | 40 → 14 | 0.0003 |

Findings:
- **Each design is hit hardest by its own adversary.** That is the diagonal of the first table.
  Its worst-case miss rates are 2.67 % for the final design, 1.39 % for contact order 1 and
  0.93–0.95 % for exact Taylor n = 3.
- **The robustness advantage survives, but it is smaller than §6 suggested.**
  - Exact Taylor n = 3 beats the final design under every adversary (p ≤ 0.008).
  - Contact order 1 cuts the worst case roughly in half. Under the adversary tuned against it, its
    advantage over the final design is not significant (36 → 24, p = 0.15).
- **The adversary is weak in absolute terms.** Tuned against contact1 or exact3, it realized its
  intended class in only 6/200 and 12/200 cases. Every engagement is still intercepted.
- **In hits delivered, the final design wins every cell.**

  | adversary | final | contact order 1 | exact Taylor n = 3 |
  |---|---|---|---|
  | final | **1783** | 1709 | 1566 |
  | contact order 1 | **1793** | 1699 | 1561 |
  | exact Taylor n = 3 | **1786** | 1617 | 1488 |

## 7a. Assessment

1. **The order-4 cap is real, but it binds only on the transient.** Contact order is
   effectively capped at 4–5 in a stable loop. Exact Taylor stops at n = 4, and ζ = 0.107 there.
   In every mode the deviation obeys P(D)(y_d − θ̂) = 0 (§3), so in steady state the loop tracks
   the nowcast and the prediction cancels regardless of P. The cap restricts which deviation
   poles may be chosen when P = A_n. The one design where it truly binds (exact Taylor) is clearly
   worse on the random engagements: 2.8 m vs 0.56 m pointing RMS.
2. **Prediction in the loop is a choice of deviation poles,** and fidelity bounds them:
   σ ≤ N/t_p. The final design gives up fidelity to get 20 rad/s poles. Contact order 1 attains
   the bound (6 rad/s) with the velocity prediction kept exact.
3. **The pole choice trades pointing accuracy for a lower per-pulse miss rate under attack, not
   for more hits.**

   | Design | MC pointing RMS | MC misses | worst-case adversarial miss rate | hits delivered under its own adversary (200 seeds) |
   |---|---|---|---|---|
   | final | 0.56 m | 4 / 9197 | 2.67 % | **1783** |
   | contact order 1 | 0.99 m | 24 / 8931 | 1.39 % | 1699 |
   | exact Taylor n = 3 | 2.82 m | 23 / 8056 | 0.95 % | 1488 |

   All three intercept every one of the 1000 + 1000 + 600 engagements. On the admissible random
   engagements, the final design is clearly best. Under attack it misses more often per pulse,
   but it still delivers the most hits and hits first at the same time. The slow designs' lower
   miss rate comes from shooting less. **The final design stays the default.**
4. **Blends are a smooth dial, but not a useful one.** From α ≈ 0.2 they behave like the final
   design, including its exposure to the adversary.
5. **The homotopy fire planner is now principled,** with event signatures and the Eq. 12 plan
   consistency, but it does not beat the mixture rule (MC 4 → 11, HG 244 → 230, both n.s.). The
   firing decision is not the bottleneck; the deviation poles are.
6. **Robustness checks.** A ±20 % plant error and removing class re-identification have no
   measurable effect.

## 8. Reproduce

```bash
cd python
python scripts/research_contact_order.py
python scripts/research_pole_shift.py
sh scripts/run_research.sh          # HG cache (~40 min) + MC/HG campaign, 12 cores
sh scripts/run_research2.sh         # contact-order modes
sh scripts/run_deconfound.sh        # adversary re-tuned against contact1 / exact3 (200 seeds)
python scripts/research_eval.py
python scripts/research_deconfound.py
python -m pytest -q
```
