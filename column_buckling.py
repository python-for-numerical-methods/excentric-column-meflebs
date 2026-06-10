import numpy as np
from scipy.optimize import bisect


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
        # sigma_max(P) - sigma_allow = 0
        arg = (L / (2.0 * r)) * np.sqrt(P / (E * A))
        sec_val = 1.0 / np.cos(arg)
        sigma_max = (P / A) * (1.0 + (e * c / (r ** 2)) * sec_val)
        return sigma_max - sigma_allow

    # הגבול העליון הוא עומס אוילר - שם sec שואף לאינסוף
    # P_euler = pi^2 * E * I / L^2 = pi^2 * E * A * r^2 / L^2
    P_euler = (np.pi ** 2 * E * A * r ** 2) / (L ** 2)

    # גבול תחתון - ערך קטן מאוד
    P_low = 1e-3

    # גבול עליון - 99.9% מעומס אוילר כדי להימנע מאסימפטוטה
    P_high = 0.999 * P_euler

    # וידוא שיש שינוי סימן בתחום
    f_low = f(P_low)
    f_high = f(P_high)

    # אם f(P_low) > 0 אז גם בעומס מינימלי המאמץ כבר עובר את המותר
    # אם f(P_high) < 0 אז גם בעומס מקסימלי המאמץ לא מגיע למותר
    if f_low * f_high > 0:
        # ננסה להרחיב את הגבול העליון
        P_high = 0.9999 * P_euler
        f_high = f(P_high)

    P_critical = bisect(f, P_low, P_high, xtol=1e-6, rtol=1e-9)

    return float(P_critical)
