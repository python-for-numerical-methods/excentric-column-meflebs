import numpy as np
from scipy import optimize


def column_stress_error(P, L, E, A, r, c, e, sigma_allow):
    """
    מחשב את ההפרש בין המאמץ המקסימלי (לפי נוסחת הסקנט) לבין המאמץ המותר.
    כאשר הפונקציה מחזירה 0, P הוא העומס הקריטי.
    """
    arg = (L / (2.0 * r)) * np.sqrt(P / (E * A))
    sec_val = 1.0 / np.cos(arg)
    sigma_max = (P / A) * (1.0 + (e * c / (r ** 2)) * sec_val)
    return sigma_max - sigma_allow


def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    מוצא את העומס הקריטי P שבו המאמץ המקסימלי שווה למאמץ המותר.

    L: אורך העמוד (מ"מ)
    E: מודול אלסטיות (MPa)
    A: שטח חתך (ממ"ר)
    r: רדיוס אינרציה (מ"מ)
    c: מרחק לסיב קיצוני (מ"מ)
    e: אקסצנטריות (מ"מ)
    sigma_allow: מאמץ מותר (MPa)

    Return: P_critical (N)
    """
    # ניחוש התחלתי - חצי מעומס אוילר
    P_euler = (np.pi ** 2 * E * A * r ** 2) / (L ** 2)
    P_guess = 0.5 * P_euler

    P_critical = optimize.newton(
        lambda P: column_stress_error(P, L, E, A, r, c, e, sigma_allow),
        P_guess
    )

    return float(P_critical)
