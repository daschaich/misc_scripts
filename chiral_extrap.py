#!/usr/bin/python
import os
import sys
import numpy as np
from scipy import special
# ------------------------------------------------------------------
# Perform chiral extrapolation (e.g., for Z_A data)

# Parse arguments: first is the file to analyze (FORMAT: m dat err)
# Optional second argument is order of polynomial fit
if len(sys.argv) < 2:
  print("Usage:", str(sys.argv[0]), "<file> <poly degree (default=1)>")
  sys.exit(1)
filename = str(sys.argv[1])
if len(sys.argv) > 2:
  degree = int(sys.argv[2])
else:
  degree = 1

if not os.path.isfile(filename):
  print("ERROR:", filename, "does not exist")
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read, parse and fit data
# Assumed format: m dat err
mList = []
datList = []
wList = []
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  mList.append(float(temp[0]))
  datList.append(float(temp[1]))
  wList.append(np.power(float(temp[2]), -1))    # Squared in fit...
m = np.array(mList)
dat = np.array(datList)
weight = np.array(wList)

# Handle case of no degrees of freedom
dof = len(m) - degree - 1
if dof == 0:
  out = np.polyfit(m, dat, degree, full=False, w=weight, cov=False)
  print("%.6g (no dof)" % out[degree])
  sys.exit(0)

# Return polynomial coefficient values and covariance matrix
# but not diagnostic information and singular value decomposition
# Intercept is last polynomial coefficient
out, cov = np.polyfit(m, dat, degree, full=False, w=weight, cov=True)

# This is annoying: polyfit scales the covariance matrix by chiSq_dof
# so that the weights are considered relative
# That is, the errors' overall scale doesn't affect the final uncertainty
# For the record, gnuplot and Mathematica do the same thing...
# We need to divide the full covariance matrix by chiSq_dof to correct this,
# or equivalently divide the final error by sqrt(chiSq_dof),
# since the absolute scale of our uncertainties is meaningful
# !!! The easy way to check this is to scale all the errors by (e.g.) 1000x
# !!! and see whether the final uncertainty increases, as it should
fit = np.poly1d(out)
chi = (dat - fit(m)) * weight
chiSq = (chi**2).sum()
chiSq_dof = chiSq / float(dof)
cov /= chiSq_dof
CL = 1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq)
intercept = out[degree]
int_err = np.sqrt(cov[degree][degree])    # Component of covar matrix
print("0 %.8g %.4g # %.4g/%d = %.4g --> %.4g" \
      % (intercept, int_err, chiSq, dof, chiSq_dof, CL))
#print("%.6g %.4g" % (out[0], np.sqrt(cov[0][0])))
# ------------------------------------------------------------------
