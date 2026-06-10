import numpy as np
from scipy.optimize import bisect

def find_critical_load(L, E, A, r, c, e, sigma_allow):

    def f(P):
        angle = (L / (2 * r)) * np.sqrt(P / (E * A))

        sigma_max = (P / A) * (
            1 + (e * c / r**2) * (1 / np.cos(angle))
        )

        return sigma_max - sigma_allow

    p_min = 1e-5
    p_max = sigma_allow * A

    while f(p_max) < 0:
        p_max *= 2

    return float(bisect(f, p_min, p_max))
