import numpy as np
from scipy.optimize import newton


def find_critical_load(L, E, A, r, c, e, sigma_allow):

    def f(P):
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
        sigma = (P / A) * (
            1 + (e * c / r**2) * (1 / np.cos(theta))
        )
        return sigma - sigma_allow

    P0 = sigma_allow * A
    return newton(f, P0)
