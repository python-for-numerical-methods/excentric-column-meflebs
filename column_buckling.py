import numpy as np
from scipy import optimize


def column_stress_error(P, L, E, A, r, c, e, sigma_allow):
    """
    פונקציית עזר: מחשבת f(P) = sigma_max(P) - sigma_allow
    כאשר f(P) = 0, מצאנו את העומס הקריטי.
    """
    arg = (L / (2.0 * r)) * np.sqrt(P / (E * A))
    sec_val = 1.0 / np.cos(arg)
    sigma_max = (P / A) * (1.0 + (e * c / (r ** 2)) * sec_val)
    return sigma_max - sigma_allow


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
    # עומס אוילר - הגבול העליון התיאורטי
    P_euler = (np.pi ** 2 * E * A * r ** 2) / (L ** 2)

    # גבולות לשיטת החצייה
    P_low = 1e-6
    P_high = 0.999 * P_euler

    # מציאת השורש בשיטת החצייה
    P_critical = optimize.bisect(
        lambda P: column_stress_error(P, L, E, A, r, c, e, sigma_allow),
        P_low,
        P_high
    )

    return float(P_critical)
