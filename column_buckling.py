import numpy as np
from scipy import optimize


def column_stress_error(P, L, E, A, r, c, e, sigma_allow):
    """
    מחזירה את f(P) = σ_max(P) - σ_allow
    """

    if P <= 0:
        return -sigma_allow

    # θ = (L / 2r) * sqrt(P / (EA))
    theta = (L / (2 * r)) * np.sqrt(P / (E * A))

    cos_theta = np.cos(theta)

    # הגנה מנומרית: sec(x) מתפוצץ כש-cos(x)=0
    if np.abs(cos_theta) < 1e-12:
        return np.inf

    sec_theta = 1.0 / cos_theta

    # σ_max = (P/A) * [1 + (ec/r²)*sec(θ)]
    sigma_max = (P / A) * (1 + (e * c / (r ** 2)) * sec_theta)

    return sigma_max - sigma_allow


def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    פתרון עומס קריטי לפי שיטת החצייה (Bisection)
    תוך שימוש ב-lambda כדי לבצע אופטימיזציה רק על P.
    """

    # פונקציית השורש f(P) עם כל הפרמטרים הקבועים
    f = lambda P: column_stress_error(P, L, E, A, r, c, e, sigma_allow)

    # טווח תחתון
    P_low = 0.0

    # אם כבר ב-P=0 המאמץ גדול מהמותר — אין פתרון
    if f(P_low) > 0:
        return 0.0

    # טווח עליון — נגדיל עד שינוי סימן
    P_high = 1.0
    while f(P_high) < 0:
        P_high *= 2
        if P_high > 1e12:  # הגנה מפני לולאה אינסופית
            raise ValueError("לא נמצא טווח מתאים לשיטת החצייה")

    # פתרון בשיטת החצייה
    P_crit = optimize.bisect(f, P_low, P_high, xtol=1e-6, rtol=1e-6)

    return float(P_crit)
