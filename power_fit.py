#!/usr/bin/python
import os
import sys
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Perform fit to power law (e.g., for M_P data)

# Parse arguments: first is the file to analyze (FORMAT: m dat err)
# Optional second argument fixes the power to 1 / (1 + gamma_m)
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file> <fix gamma_m>"
  sys.exit(1)
filename = str(sys.argv[1])
if len(sys.argv) > 2:
  power = 1.0 / (1.0 + float(sys.argv[2]))
  # For MP^2
#  power = 2. / (1. + float(sys.argv[2]))
  # For topological susceptibility
#  power = 4. / (1. + float(sys.argv[2]))
  # For more general fits, just read in power itself
#  power = float(sys.argv[2])
else:
  power = 999

if not os.path.isfile(filename):
  print "ERROR:", filename, "does not exist"
  sys.exit(1)
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

if power == 999:
  # Fit to A*x**B
  fitfunc = lambda p, x: p[0] * np.power(x, p[1])
  errfunc = lambda p, x, y, err: (fitfunc(p, x) - y) / err
  p_in = [1., 1.]
  dof = len(m) - 2
else:
  # Fit to A*x**power
  fitfunc = lambda p, x: p[0] * np.power(x, power)
  errfunc = lambda p, x, y, err: (fitfunc(p, x) - y) / err
  p_in = [1.]
  dof = len(m) - 1

# For now, demand that we have degrees of freedom
if dof == 0:
  print "ERROR: dof > 0 required for now"
  sys.exit(1)

all_out = optimize.leastsq(errfunc, p_in[:], args=(m, dat, err),
                           full_output = 1)
p_out = all_out[0]
covar = all_out[1]

print "coeff %.6g %.4g" % (p_out[0], np.sqrt(covar[0][0]))
if power == 999:
  x = float(p_out[1])
  x_err = np.sqrt(covar[1][1])
  gamma = (1.0 - x) / x       # Power = 1 / (1 + gamma_m)
  ga_err = x_err / x**2
  # For MP^2
#  gamma = (2.0 - x) / x      # Power = 2 / (1 + gamma_m)
#  ga_err = 2.0 * x_err / x**2
  # For topological susceptibility
#  gamma = (4.0 - x) / x      # Power = 4 / (1 + gamma_m)
#  ga_err = 4.0 * x_err / x**2
  print "power %.6g %.4g --> %.6g %.4g" \
        % (x, x_err, gamma, ga_err)

chiSq = ((errfunc(p_out, m, dat, err))**2).sum()
print "chiSq/dof = %.4g/%d = %.4g" % (chiSq, dof, chiSq / dof)
# ------------------------------------------------------------------
