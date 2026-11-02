"""Specification constants shared by the simulation harness (Sections 2-3 of the assignment).

The controller keeps its own copy of the constants it needs, so it can be ported to a
single self-contained MATLAB file.
"""
import numpy as np

# Engagement
T_END = 10.0            # standardized engagement duration [s]
TP = 0.50               # prescribed prediction horizon [s]
R_ENG = 2.0             # engagement radius [m]
T_RECHARGE = 1.0        # minimum time between laser pulses [s]

# Radar
DT_RADAR = 0.05         # 20 Hz measurement interval [s]
SIGMA_R = 1.0           # range noise std [m]
SIGMA_TH = np.deg2rad(0.02)  # elevation noise std [rad]

# Target envelope
R_MIN, R_MAX = 800.0, 2000.0
V_MAX_TGT = 150.0
A_MAX_TGT = 20.0
T_BO_MAX = 0.50         # maximum measurement-loss duration [s]
T_BOO_MIN = 1.00        # earliest measurement-loss onset [s]

# Turret
J_TURRET = 1200.0       # [kg m^2]
J_MOTOR = 0.20          # [kg m^2]
N_GEAR = 100.0
B_TURRET = 1500.0       # turret-side viscous damping [N m s/rad]
K_T = 10.0              # torque constant [N m/A]
K_E = 10.0              # back-emf constant [V s/rad]
R_ARM = 1.0             # armature resistance [ohm]
L_ARM = 0.10            # armature inductance [H]
J_EQ = J_TURRET + N_GEAR**2 * J_MOTOR   # 3200 kg m^2, reflected to the turret
TH_MIN = np.deg2rad(-10.0)
TH_MAX = np.deg2rad(85.0)

# Power stage assumption (docs/power_stage.md)
V_SUPPLY_MAX = 600.0    # [V]
