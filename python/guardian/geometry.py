"""Engagement geometry used by the harness to judge laser pulses."""
import numpy as np


def beam_miss_distance(theta_beam, x, y):
    """Minimum distance from the target point (x, y) to a ray from the origin at theta_beam."""
    ux, uy = np.cos(theta_beam), np.sin(theta_beam)
    along = x * ux + y * uy
    if along < 0.0:
        return float(np.hypot(x, y))
    return float(abs(x * uy - y * ux))
