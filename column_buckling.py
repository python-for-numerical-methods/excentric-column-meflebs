import numpy as np
from scipy.optimize import bisect


def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: length [mm]
    E: modulus of elasticity [MPa = N/mm²]
    A: area [mm²]
    r: radius of gyration [mm]
    c: distance to extreme fiber [mm]
    e: eccentricity [mm]
    sigma_allow: allowable stress [MPa]

    Returns:
        Critical load P [N]
    """

    def sigma_max(P):
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
        sec_theta = 1.0 / np.cos(theta)

        return (P / A) * (
            1 + (e * c / r**2) * sec_theta
        )

    def f(P):
        return sigma_max(P) - sigma_allow

    # Euler load gives a safe upper bound before sec() becomes singular
    P_euler = (np.pi**2 * E * A * r**2) / (L**2)

    lower = 0.0
    upper = 0.99 * P_euler

    return bisect(f, lower, upper, xtol=1e-6)
