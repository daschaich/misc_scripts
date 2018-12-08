#!/usr/bin/python3
import os
import sys
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammainc
# ------------------------------------------------------------------
# Fit susceptibility data to power law, A * x**B
# May make power_fit.py a bit redundant
# TODO: May eventually compute and print out error band for gnuplotting
#       But ignoring that for now

# Parse argument: the file to analyze (FORMAT: L dat err)
if len(sys.argv) < 2:
  print("Usage:", str(sys.argv[0]), "<file>")
  sys.exit(1)
filename = str(sys.argv[1])

if not os.path.isfile(filename):
  print("ERROR:", filename, "does not exist")
  sys.exit(1)

# errfunc will be minimized via least-squares optimization
# p_in are order-of-magnitude initial guesses
expfunc = lambda p, x: p[0] * np.power(x, p[1])
errfunc = lambda p, x, y, err: (expfunc(p, x) - y) / err
p_in = np.array([1.0, 1.0])

# Define corresponding Jacobian matrix
# Recall x^p = exp(p * log x)
def jac(p, x, y, err):
  J = np.empty((x.size, p.size), dtype = np.float)
  J[:, 0] = np.power(x, p[1])
  J[:, 1] = p[0] * np.log(x) * np.power(x, p[1])
  for i in range(p.size):
    J[:, i] /= err
  return J
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read, parse and fit data
# Assumed format: L dat err
LList = []
datList = []
errList = []
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  LList.append(float(temp[0]))
  datList.append(float(temp[1]))
  errList.append(float(temp[2]))
L = np.array(LList)
dat = np.array(datList)
err = np.array(errList)

# Demand that we have degrees of freedom
dof = len(L) - len(p_in)
if dof < 1:
  print("ERROR: dof > 0 required")
  sys.exit(1)

# Extract and save fit parameters
# and covariance matrix (J^T J)^{-1} where J is jacobian matrix
# method='lm' is Levenberg--Marquardt (can't handle bounds)
all_out = least_squares(errfunc, p_in, jac=jac, method='lm',
                        args=(L, dat, err))
p_out = all_out.x
tj = all_out.jac
cov = np.linalg.inv(np.dot(np.transpose(tj), tj))

if all_out.success < 0 or all_out.success > 4:
  print("ERROR: Fit failed with the following error message")
  print(errmsg)
  sys.exit(1)

print("Power: %.6g %.4g" % (p_out[1], np.sqrt(cov[1][1])))

# Compute chiSq and confidence level of fit
chiSq = ((errfunc(p_out, L, dat, err))**2).sum()
CL = 1.0 - gammainc(0.5 * dof, 0.5 * chiSq)
print("chiSq/dof = %.4g/%d = %.4g --> CL = %.4g" \
      % (chiSq, dof, chiSq / dof, CL))

# Format to copy+paste into gnuplot: fit and error function
# TODO: Ignore error bands for now
print("Fit: %.4g * x**%.4g" % (p_out[0], p_out[1]))
#print("Err: sqrt(%.4g - %.4g * x + %.4g * x**2)" \
#      % (cov[0][0], -1.0 * (cov[1][0] + cov[0][1]), cov[1][1]))
# ------------------------------------------------------------------
