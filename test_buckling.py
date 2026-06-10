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
        angle = (L / (2 * r)) * np.sqrt(P / (E * A))

        sigma_max = (P / A) * (
            1 + (e * c / (r ** 2)) * (1 / np.cos(angle))
        )

        return sigma_max - sigma_allow

    # גבול תחתון
    p_min = 1e-8

    # התחלה עם עומס לחיצה ישיר
    p_max = sigma_allow * A

    # הרחבת הגבול העליון עד שיש שינוי סימן
    while f(p_max) <= 0:
        p_max *= 2

    # מציאת השורש
    return float(bisect(f, p_min, p_max))
