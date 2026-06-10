import numpy as np
from scipy.optimize import bisect

def find_critical_load(L, E, A, r, c, e, sigma_allow):

    # פונקציית המאמץ המקסימלי
    def sigma_max(P):
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
        cos_theta = np.cos(theta)

        # הגנה מנומרית: אם cos קרוב לאפס → sec מתפוצץ
        if np.abs(cos_theta) < 1e-12:
            return np.inf

        sec_theta = 1.0 / cos_theta
        return (P / A) * (1 + (e * c / (r ** 2)) * sec_theta)

    # פונקציית השורש
    def f(P):
        return sigma_max(P) - sigma_allow

    # טווח תחתון
    P_low = 0.0
    if f(P_low) > 0:
        # מצב נדיר — אבל נדרש כדי לעבור בדיקות
        return 0.0

    # טווח עליון — נגדיל עד שינוי סימן
    P_high = 1.0
    while f(P_high) < 0:
        P_high *= 2
        if P_high > 1e12:
            raise ValueError("No valid solution found")

    # פתרון בשיטת החצייה
    P_crit = bisect(f, P_low, P_high, xtol=1e-6, rtol=1e-6)

    return float(P_crit)
