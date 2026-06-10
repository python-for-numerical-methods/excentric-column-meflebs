import numpy as np
from scipy.optimize import bisect

def find_critical_load(L, E, A, r, c, e, sigma_allow):

    def sigma_max(P):
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
        sec_theta = 1.0 / np.cos(theta)

        return (P / A) * (
            1.0 + (e * c / r**2) * sec_theta
        )

    def f(P):
        return sigma_max(P) - sigma_allow

    p_min = 1e-8
    p_max = sigma_allow * A

    while f(p_max) < 0:
        p_max *= 2

    return bisect(f, p_min, p_max)
