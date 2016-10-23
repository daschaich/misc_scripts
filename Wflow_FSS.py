#!/usr/bin/python
import os
import sys
import numpy as np
from scipy import optimize, special
# ------------------------------------------------------------------
# Perform fit to constant plus power for Wilson flow data

# Parse arguments: first is the file to analyze (FORMAT: L dat err)
# Second argument is minimum L to include in the fit
# Optional second argument fixes the constant
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<file> <Lmin> <fix gstar>"
  sys.exit(1)
filename = str(sys.argv[1])
Lmin = int(sys.argv[2])
if len(sys.argv) > 3:
  gstar = float(sys.argv[3])
else:
  gstar = -1.0      # Use as switch in fits below

if not os.path.isfile(filename):
  print "ERROR:", filename, "does not exist"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read, parse and fit data
# Assumed format: m dat err
LList = []
datList = []
errList = []
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  tL = float(temp[0])
  if tL >= Lmin:
    LList.append(tL)
    if gstar < 0:
      # Fit g^2(L) = gstar^2 + A*L**gamma
      datList.append(float(temp[1]))
    else:
      # Fit (g^2(L) - gstar^2) = A*L**gamma with gstar^2 fixed
      datList.append(float(temp[1]) - gstar)
    errList.append(float(temp[2]))
L = np.array(LList)
dat = np.array(datList)
err = np.array(errList)

if gstar < 0:
  # Fit g^2(L) = gstar^2 + A*L**gamma
  fitfunc = lambda p, x: p[0] + p[1] * np.power(x, p[2])
  errfunc = lambda p, x, y, err: (fitfunc(p, x) - y) / err
  p_in = np.array([1.0, 1.0, 1.0], dtype = np.float)
else:
  # Fit (g^2(L) - gstar^2) = A*L**gamma with gstar^2 fixed
  fitfunc = lambda p, x: p[0] * np.power(x, p[1])
  errfunc = lambda p, x, y, err: (fitfunc(p, x) - y) / err
  p_in = np.array([1.0, 1.0], dtype = np.float)


# Define corresponding Jacobian matrix
def jac(p, x, y, err):
  J = np.empty((x.size, p.size))
  if gstar < 0:
    J[:, 0] = 1.0 / err
    J[:, 1] = np.power(x, p[2]) / err
    J[:, 2] = p[1] * np.log(x) * np.power(x, p[2]) / err
  else:
    J[:, 0] = np.power(x, p[1]) / err
    J[:, 1] = p[0] * np.log(x) * np.power(x, p[1]) / err
  return J

# For now, demand that we have degrees of freedom
dof = len(L) - len(p_in)
if dof <= 0:
  print "ERROR: dof > 0 required for now"
  sys.exit(1)

all_out = optimize.least_squares(errfunc, p_in,
                                 jac=jac, method='lm', args=(L, dat, err))
out = all_out.x
tj = all_out.jac
pcov = np.linalg.inv(np.dot(np.transpose(tj), tj))

if all_out.success < 0 or all_out.success > 4:
  print "WARNING: Fit failed with the following error message:"
  print errmsg

# Compute chiSq and print exponent -- the last component of the output
chiSq = ((errfunc(out, L, dat, err))**2).sum()
CL = 1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq)
print "gamma %.6g %.4g" % (out[-1], np.sqrt(pcov[-1][-1]))
if gstar < 0:
  # Print estimated location of IRFP
  IRFP = float(out[0])
  IRFP_err = np.sqrt(pcov[0][0])
  print "gstar %.6g %.4g" % (IRFP, IRFP_err)
print "chiSq/dof = %.4g/%d = %.4g --> %.4g" % (chiSq, dof, chiSq / dof, CL)
# ------------------------------------------------------------------
