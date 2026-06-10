import numpy as np
from scipy.optimize import bisect


def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    מוצא את העומס המקסימלי P שיביא את העמוד למאמץ המותר.
    
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa
    
    Return: העומס P בניוטון (float)
    """
    
    def equation(P):
        arg = (L / (2 * r)) * np.sqrt(P / (E * A))
        secant_value = 1.0 / np.cos(arg)
        sigma_max = (P / A) * (1 + (e * c / r**2) * secant_value)
        return sigma_max - sigma_allow
    
    P_euler = (np.pi**2 * E * A) / (L / r)**2
    P_low = 1.0
    P_high = 0.99 * P_euler
    
    P_critical = bisect(equation, P_low, P_high, xtol=1e-6)
    
    return P_critical
