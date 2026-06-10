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
    # ניחוש התחלתי - מבוסס על מאמץ מותר כפול שטח חתך (כ-50% מעומס אוילר)
    P_euler = (np.pi ** 2 * E * A * r ** 2) / (L ** 2)
    P_guess = 0.5 * P_euler

    # מציאת השורש בשיטת ניוטון-רפסון
    P_critical = optimize.newton(
        lambda P: column_stress_error(P, L, E, A, r, c, e, sigma_allow),
        P_guess
    )

    return float(P_critical)
