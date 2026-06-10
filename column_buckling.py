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

    # -----------------------------
    # פונקציית המאמץ המקסימלי σ_max(P)
    # -----------------------------
    def sigma_max(P):
        if P <= 0:
            return 0.0

        # θ = (L / 2r) * sqrt(P / (EA))
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))

        cos_theta = np.cos(theta)

        # הגנה מנומרית: sec(x) מתפוצץ כש-cos(x)=0
        if np.abs(cos_theta) < 1e-12:
            return np.inf

        sec_theta = 1.0 / cos_theta

        # σ_max = (P/A) * [1 + (ec/r²)*sec(θ)]
        return (P / A) * (1 + (e * c / (r ** 2)) * sec_theta)

    # -----------------------------
    # פונקציית השורש f(P) = σ_max(P) - σ_allow
    # -----------------------------
    def f(P):
        return sigma_max(P) - sigma_allow

    # -----------------------------
    # שלב 1: קביעת טווח חיפוש
    # -----------------------------
    P_low = 0.0

    # אם כבר ב־P=0 המאמץ גדול מהמותר — אין פתרון
    if f(P_low) > 0:
        return 0.0

    # נתחיל מגבול עליון קטן ונגדיל עד שינוי סימן
    P_high = 1.0
    while f(P_high) < 0:
        P_high *= 2
        if P_high > 1e12:  # הגנה מפני לולאה אינסופית
            raise ValueError("לא נמצא טווח מתאים לשיטת החצייה")

    # -----------------------------
    # שלב 2: פתרון בשיטת החצייה
    # -----------------------------
    P_crit = bisect(f, P_low, P_high, xtol=1e-6, rtol=1e-6)

    return float(P_crit)
