import math
from scipy.optimize import bisect


def find_critical_load(L, E, A, r, c, e, sigma_allow):

    euler_limit = E * A * ((math.pi * r) / L) ** 2

    def secant_equation(P):

        alpha = (L / (2.0 * r)) * math.sqrt(P / (E * A))

        cos_alpha = math.cos(alpha)

        if abs(cos_alpha) < 1e-12:
            return float("inf")

        sigma_max = (P / A) * (
            1.0 + (e * c) / (r**2 * cos_alpha)
        )

        return sigma_max - sigma_allow

    return bisect(
        secant_equation,
        0.0,
        0.9999 * euler_limit,
        xtol=1e-6
    )
