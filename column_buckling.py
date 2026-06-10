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

    # פונקציית המאמץ המקסימלי
    def sigma_max(P):
        if P <= 0:
            return -sigma_allow  # כדי שהשורש לא ייפול על 0

        theta = (L / (2 * r)) * np.sqrt(P / (E * A))

        # הגנה מפני sec בנקודות אסורות
        cos_theta = np.cos(theta)
        if np.abs(cos_theta) < 1e-12:
            return np.inf

        sec_theta = 1.0 / cos_theta

        return (P / A) * (1 + (e * c / (r ** 2)) * sec_theta)

    # פונקציית השורש
    def f(P):
        return sigma_max(P) - sigma_allow

    # טווח חיפוש: מתחילים נמוך ומגדילים עד שינוי סימן
    P_low = 1e-6
    P_high = 1.0

    # מגדילים את הגבול העליון עד שהמאמץ עובר את המותר
    while f(P_high) < 0:
        P_high *= 2
        if P_high > 1e12:  # הגנה מפני לולאה אינסופית
            raise ValueError("לא נמצא פתרון בתחום סביר")

    # פתרון בשיטת החצייה
    P_crit = bisect(f, P_low, P_high, xtol=1e-6, rtol=1e-6)

    return float(P_crit)
