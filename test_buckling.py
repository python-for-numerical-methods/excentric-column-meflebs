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
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))

        # sec(theta) = 1/cos(theta)
        sec_theta = 1.0 / np.cos(theta)

        sigma_max = (P / A) * (
            1 + (e * c / r**2) * sec_theta
        )

        return sigma_max - sigma_allow

    # גבולות התחלתיים לחיפוש
    P_low = 1.0
    P_high = sigma_allow * A

    # הרחבת התחום במידת הצורך
    while f(P_low) * f(P_high) > 0:
        P_high *= 2

    return bisect(f, P_low, P_high)
