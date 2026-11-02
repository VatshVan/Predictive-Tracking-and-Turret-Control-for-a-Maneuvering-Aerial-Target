# Elevation-Drive Power Stage (design assumption)

The assignment does not specify a supply or a voltage limit for V_m. We assume the turret drive is fed from a **standard 230 V, 50 Hz single-phase outlet** (IS 1293, 16 A) and derive the converter that produces the **±600 V** motor voltage used as the saturation limit in the controller. The laser's own power supply is outside the scope of this work.

## 1. Why ±600 V

The drive is back-EMF dominated (turret-side constant N·K_e = 1000 V·s/rad):

```
V_m = R·i + N·K_e·θ̇        i = (J_eq·θ̈ + b·θ̇) / (N·K_t)
J_eq = 3200 kg·m²,  b = 1500 N·m·s/rad,  N·K_t = 1000 N·m/A,  R = 1 Ω
```

| Operating point | θ̇ [rad/s] | θ̈ [rad/s²] | i [A] | V_m [V] | P [W] |
|---|---|---|---|---|---|
| Worst-case tracking (v = 150 m/s at r = 800 m) | 0.1875 | 0.060 | 0.47 | 188 | 89 |
| Slew cruise | 0.60 | 0 | 0.90 | 601 | 541 |
| Slew + acceleration (design point) | 0.60 | 1.0 | 4.10 | 604 | 2477 |
| Slew + hard acceleration | 0.60 | 2.0 | 7.30 | 607 | 4433 |

- Worst-case LOS rate: v_max / r_min = 150/800 = 0.1875 rad/s.
- Worst-case LOS acceleration: a_max/r_min + v_max²/r_min² = 0.060 rad/s².
- **±600 V gives 3.2× headroom over worst-case tracking** and a slew rate of ≈ 0.6 rad/s (34°/s). That lets the turret acquire from an arbitrary initial angle within the 10 s engagement.

## 2. Converter chain

```
230 V AC ─► EMI filter ─► bridge rectifier ─► PFC boost ─► DC link 700 V ─► H-bridge (4-quadrant PWM) ─► motor
             (CISPR)       (≈325 V peak)        (PF≈0.99)   C_dc + brake chopper    |V_m| ≤ 600 V
```

| Stage | Function | Rating (design) |
|---|---|---|
| EMI filter | Conducted-emission compliance | 16 A, 250 V AC |
| Bridge rectifier | AC → pulsating DC, 325 V peak | 600 V, 25 A |
| PFC boost converter | Boost to 700 V DC; sinusoidal input current (IEC 61000-3-2 harmonic limits) | 3 kW, 1200 V SiC MOSFET + SiC diode, ~65 kHz |
| DC-link capacitor | Buffers acceleration peaks and absorbs regenerated energy | see §3 |
| Brake chopper | Dumps regenerated energy if the link exceeds 760 V | 1200 V switch + 50 Ω resistor |
| H-bridge | Bidirectional voltage and current (drive + regenerative braking); bipolar PWM ~20 kHz ≫ L/R corner (1.6 Hz) and the 8.9 Hz mechanical mode | 4 × 1200 V / 20 A SiC MOSFET |

Why 700 V on the link: the H-bridge can only produce |V_m| < V_dc. That leaves ~100 V for switch drops and PWM duty-cycle limits (max duty ≈ 0.9). 1200 V devices then have a > 1.5× voltage margin.

## 3. Sizing checks

**Input current at the design point (2.5 kW):**
I_in = P / (V_ac · η · PF) = 2477 / (230 · 0.95 · 0.99) ≈ **11.5 A < 16 A outlet**. Worst-case continuous tracking draws < 0.5 A. At this design point the drive runs from a household socket.

**Regenerated energy when stopping from full slew:**
E = ½·J_eq·θ̇² = ½ · 3200 · 0.6² = **576 J**.
- Absorbed by capacitance alone with a 700 → 760 V rise: C = 2E / (760² − 700²) = 1152 / 87 600 ≈ **13 mF** (bulky).
- **Chosen:** C_dc = 2 mF plus a brake chopper. The capacitor holds the link during transients and the chopper dissipates the rest. Peak chopper power ≈ 760²/50 ≈ 11.6 kW for < 0.1 s.

**Single-phase 100 Hz ripple** at 2.5 kW with ±5 % ripple (ΔV = 35 V):
C ≥ P / (2π · 100 Hz · V_dc · ΔV) = 2477 / (2π · 100 · 700 · 35) ≈ **0.16 mF**. Covered by the 2 mF bank.

**Current limit:** the H-bridge current limit is set at **10 A**, which is 2.4× the design point. The controller only commands voltage, so the ±600 V saturation is enforced in the controller with anti-windup. The current limit is a hardware protection layer, not part of the control law.

## 4. What enters the control design

- Actuator saturation: **|V_m| ≤ 600 V** (hard limit in `turret_controller`, with anti-windup).
- PWM at ~20 kHz is far above every plant mode, so the bridge is modelled as an ideal voltage source (average model). Its delay of ~25 µs is negligible against the 8.9 Hz mode.
- Elevation limits −10° … +85° are handled in the reference generator, independent of the power stage.

## References
- N. Mohan, T. M. Undeland, W. P. Robbins, *Power Electronics: Converters, Applications, and Design*, 3rd ed., Wiley, 2003 (rectifier, PFC boost, H-bridge, DC link).
- IEC 61000-3-2, *Limits for harmonic current emissions (equipment input current ≤ 16 A per phase)*.
- IS 1293, *Plugs and socket-outlets of rated voltage up to and including 250 volts and rated current up to and including 16 amperes*.
