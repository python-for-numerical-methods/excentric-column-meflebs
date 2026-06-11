import numpy as np
from scipy.optimize import brentq
def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa
    Return: העומס P בניוטון (float)
    """
    def f(P):
        """f(P) = sigma_max(P) - sigma_allow"""
        if P <= 0:
            return -sigma_allow
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
        sec_term = 1.0 / np.cos(theta)
        sigma_max = (P / A) * (1.0 + (e * c) / (r ** 2) * sec_term)
        return sigma_max - sigma_allow
    # Upper bound: Euler buckling load (where sec argument -> pi/2)
    P_euler = (np.pi ** 2) * E * A * (r ** 2) / (L ** 2)
    # Ensure we are safely below the singularity
    P_low = 0.0
    P_high = 0.9999 * P_euler
    # Edge case: no eccentricity → direct solution
    if e == 0:
        return float(sigma_allow * A)
    # If f(P_high) is still negative, reduce upper bound until sign changes
    if f(P_high) < 0:
        # This can happen only for extremely small eccentricities or very low allowable stress
        # In such case the solution is very close to sigma_allow * A and below P_high
        P_high = sigma_allow * A
        if f(P_high) < 0:
            # Binary search for a valid upper bound
            factor = 2.0
            while f(P_high) < 0 and P_high < P_euler:
                P_high *= factor
                if P_high >= P_euler:
                    P_high = 0.9999 * P_euler
                    break
    P_critical = brentq(f, P_low, P_high, xtol=1e-9, rtol=1e-9)
    return float(P_critical)
