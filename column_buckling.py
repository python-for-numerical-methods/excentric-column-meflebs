import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


def find_critical_load(L, E, A, r, c, e, sigma_allow):
  """
  L: אורך במ"מ
  E: מודול אלסטיות ב-MPa
  A: שטח חתך בממ"ר
  r: רדיוס אינרציה במ"מ
  c: מרחק לסיב קיצוני במ"מ
  e: אקסצנטריות במ"מ
  sigma_allow: מאמץ מותר ב-MPa

  Return:  P בניוטון (float)
  """העומס

  def column_stress_error_internal(P):
    arg = (L / (2 * r)) * np.sqrt(P / (E * A))
    sec_term = 1 / np.cos(arg)
    sigma_max = (P / A) * (1 + (e * c / r**2) * sec_term)
    return sigma_max - sigma_allow

  initial_guess_P = 500000
  P_critical = optimize.newton(column_stress_error_internal, initial_guess_P)
  return P_critical
