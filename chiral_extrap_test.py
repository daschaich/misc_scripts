#!/usr/bin/python
import os
import sys
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammainc
# ------------------------------------------------------------------
# Perform chiral extrapolation (e.g., for Z_A data)

# Parse argument: the file to analyze
if len(sys.argv) < 2:
  print("Usage:", str(sys.argv[0]), "<file>")
  sys.exit(1)
filename = str(sys.argv[1])

if not os.path.isfile(filename):
  print("ERROR:", filename, "does not exist")
  sys.exit(1)

# errfunc will be minimized via least-squares optimization
func = lambda p, x: p[1] + p[0] * x
errfunc = lambda p, x, y, err: (func(p, x) - y) / err
p_in = np.array([0.1, 0.1])

# Define corresponding Jacobian matrix
def jac(p, x, y, err):
  J = np.empty((x.size, p.size), dtype = np.float)
  J[:, 0] = x / err
  J[:, 1] = 1.0 / err
  return J
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read, parse and fit data
# Assumed format: m dat err
mList = []
datList = []
errList = []
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  mList.append(float(temp[0]))
  datList.append(float(temp[1]))
  errList.append(float(temp[2]))
m = np.array(mList)
dat = np.array(datList)
err = np.array(errList)

# Require enough degrees of freedom
dof = len(m) - len(p_in)
if dof < 1:
  print("ERROR: dof > 0 required")
  sys.exit(1)

# Return fit parameters
# and covariance matrix (J^T J)^{-1} where J is jacobian matrix
all_out = least_squares(errfunc, p_in, #bounds=(-10.0, 10.0),
                        jac=jac, method='lm', args=(m, dat, err))
out = all_out.x
tj = all_out.jac
pcov = np.linalg.inv(np.dot(np.transpose(tj), tj))

if all_out.success < 0 or all_out.success > 4:
  print("WARNING: Fit failed with the following error message:")
  print(errmsg)

# Compute chiSq and print out intercept
chiSq = ((errfunc(out, m, dat, err))**2).sum()
CL = 1.0 - gammainc(0.5 * dof, 0.5 * chiSq)
intercept = out[1]
int_err = np.sqrt(pcov[1][1])     # Component of covariance matrix
print("0 %.8g %.4g # %.4g/%d = %.4g --> %.4g" \
      % (intercept, int_err, chiSq, dof, chiSq / dof, CL))
#print("%.6g %.4g" % (out[0], np.sqrt(pcov[0][0])))
# ------------------------------------------------------------------
