# column_buckling.py

import numpy as np
from scipy.optimize import bisect

def _sigma_max(P, L, E, A, r, c, e):
    """
    חישוב המאמץ המקסימלי לפי נוסחת הסקנט עבור עומס נתון P.
    """
    if P <= 0:
        return 0.0

    # ארגומנט הסקנט
    theta = (L / (2.0 * r)) * np.sqrt(P / (E * A))

    # sec(x) = 1 / cos(x)
    cos_theta = np.cos(theta)

    # הגנה מפני חלוקה באפס (קרבה לנקודות אי-הגדרה של sec)
    if np.isclose(cos_theta, 0.0):
        # אם הגענו קרוב מדי לנקודת אי-הגדרה, נחזיר ערך מאוד גדול
        return np.inf

    sec_theta = 1.0 / cos_theta

    sigma_max = (P / A) * (1.0 + (e * c / (r ** 2)) * sec_theta)
    return sigma_max


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

    # פונקציית השורש: f(P) = sigma_max(P) - sigma_allow
    def f(P):
        return _sigma_max(P, L, E, A, r, c, e) - sigma_allow

    # נקודת התחלה תחתונה: כמעט אפס (לא 0 כדי להימנע מבעיות נומריות)
    P_low = 1e-6
    f_low = f(P_low)

    # נניח שעבור עומס קטן המאמץ קטן מהמאמץ המותר
    # אם לא – נזיז קצת את הנקודה
    if f_low > 0:
        # במקרה קצה – אם כבר כאן המאמץ גדול מהמותר, ננסה להקטין עוד
        P_low = 0.0
        f_low = f(P_low)

    # נבנה גבול עליון ע"י הגדלה הדרגתית עד שנקבל שינוי סימן
    P_high = 1.0
    f_high = f(P_high)

    # נגדיל את P_high עד ש-f(P_high) יהיה חיובי (כלומר sigma_max > sigma_allow)
    # או עד שנגיע לגבול עליון גדול מאוד
    max_iter = 60
    iter_count = 0
    while f_high <= 0 and iter_count < max_iter:
        P_high *= 2.0
        f_high = f(P_high)
        iter_count += 1

    if f_low * f_high > 0:
        # אם לא הצלחנו למצוא טווח עם שינוי סימן – נזרוק שגיאה ברורה
        raise ValueError("לא נמצא טווח מתאים לשיטת החצייה. בדקו את הנתונים.")

    # שימוש בשיטת החצייה למציאת השורש
    P_crit = bisect(f, P_low, P_high, xtol=1e-6, rtol=1e-6, maxiter=100)

    return float(P_crit)
