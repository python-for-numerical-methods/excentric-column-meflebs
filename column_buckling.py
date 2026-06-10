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

    Return:
        P (N)
    """

    def sigma_max(P):
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
@@ -27,14 +14,10 @@ def sigma_max(P):
    def f(P):
        return sigma_max(P) - sigma_allow

    # עומס אוילר משמש חסם עליון סביר
    P_euler = (np.pi**2 * E * A * r**2) / (L**2)

    # מחפשים תחום שבו יש החלפת סימן
    P_low = 0.0
    P_high = 0.99 * P_euler
    p_min = 1e-8
    p_max = sigma_allow * A

    while f(P_high) < 0:
        P_high *= 1.5
    while f(p_max) < 0:
        p_max *= 2

    return bisect(f, P_low, P_high, xtol=1e-6)
    return bisect(f, p_min, p_max)
