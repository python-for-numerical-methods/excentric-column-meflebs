import numpy as np
from scipy.optimize import brentq

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
    
    def sigma_max(P):
        if P <= 0:
            return 0
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
        # sec(x) = 1/cos(x), and theta must be < pi/2
        sec_term = 1 / np.cos(theta)
        return (P / A) * (1 + (e * c) / (r ** 2) * sec_term)
    
    def f(P):
        return sigma_max(P) - sigma_allow
    
    # Euler buckling load - the theoretical upper limit
    P_euler = np.pi**2 * E * A * r**2 / L**2
    
    # Upper bound should be safely below Euler load to keep cos(theta) > 0
    # theta = pi/2 at P_euler, so use a fraction of it
    P_low = 0.0
    P_high = 0.999 * P_euler
    
    # Verify signs
    if f(P_low) >= 0:
        return P_low
    if f(P_high) <= 0:
        # This shouldn't happen for well-posed problems
        # But if it does, the allowable stress might be very high
        # Try a lower upper bound
        pass
    
    try:
        P_critical = brentq(f, P_low, P_high, xtol=1e-6)
    except ValueError:
        # If signs are the same, try adjusting bounds
        # Binary search for an upper bound where f > 0
        test = P_high
        while test > 1e-6 and f(test) <= 0:
            test *= 0.5
        if test <= 1e-6:
            return 0.0
        P_critical = brentq(f, P_low, test, xtol=1e-6)
    
    return float(P_critical)
