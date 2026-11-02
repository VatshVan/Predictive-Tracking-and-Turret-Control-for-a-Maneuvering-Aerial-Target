"""Log controller inputs/outputs from closed-loop Python runs for the MATLAB parity check.

python scripts/export_parity.py
Writes matlab/parity/<scenario>_<rate>Hz.csv with columns
t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret, V_m, fire_cmd.
MATLAB (matlab/parity_check.m) replays the inputs through turret_controller.m and
compares the outputs.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from guardian import turret_controller as tc  # noqa: E402
from guardian.sim import run  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'matlab', 'parity')
CASES = [(s, f) for s in ('A', 'B', 'C', 'D1', 'D4', 'D6') for f in (20, 1000)]


def main():
    os.makedirs(OUT, exist_ok=True)
    for s, f in CASES:
        rows = []

        def ctrl(t, r, th, valid, tht, thd):
            V, fire = tc.turret_controller(t, r, th, valid, tht, thd)
            rows.append((t, r, th, valid, tht, thd, V, fire))
            return V, fire

        tc.reset()
        m = run(s, f_ctrl=f, controller=ctrl, log=False)['metrics']
        path = os.path.join(OUT, f'{s}_{f}Hz.csv')
        np.savetxt(path, np.array(rows, dtype=float), delimiter=',', fmt='%.17g',
                   header='t,r_meas,theta_meas,meas_valid,theta_turret,theta_dot_turret,V_m,fire_cmd',
                   comments='')
        print(f'{path}: {len(rows)} calls, {m["n_hits"]}/{m["n_pulses"]} hits')


if __name__ == '__main__':
    main()
