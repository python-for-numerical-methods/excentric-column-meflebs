import numpy as np
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

  Return: העומס P בניוטון (float)
  """

  def column_stress_error_internal(P):
    # Calculate the argument for np.cos
    arg = (L / (2 * r)) * np.sqrt(P / (E * A))

    # Calculate the secant term (1 / cos)
    # Note: Ensure the argument to np.cos is not too close to pi/2 + k*pi, as it can cause `OverflowError`.
    # For the purpose of root finding, scipy.optimize methods can often handle large values, but
    # extreme values might cause issues. A robust solution might include boundary checks.
    sec_term = 1 / np.cos(arg)

    # Calculate the maximum stress using the secant formula
    sigma_max = (P / A) * (1 + (e * c / r**2) * sec_term)

    # Return the difference from the allowable stress (we want this to be zero)
    return sigma_max - sigma_allow

  # A good initial guess for P is crucial for newton's method.
  # Based on the plot in the notebook, a value around 500,000 N seems reasonable.
  initial_guess_P = 500000

  # Use scipy.optimize.newton to find the root of the error function
  P_critical = optimize.newton(column_stress_error_internal, initial_guess_P)
  return P_critical
