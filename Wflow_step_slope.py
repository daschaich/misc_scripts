#!/usr/bin/python
import os
import sys
import numpy as np
from scipy import special
# ------------------------------------------------------------------
# Extract slope of given step-scaling function in given range of g^2

# Parse arguments: first is the file to analyze, second is scale change s
# Third is how far from the zero to include in the fit
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<file> <s> <range_size>"
  sys.exit(1)
filename = str(sys.argv[1])
s = float(sys.argv[2])
halfrange = float(sys.argv[3])

# Just do linear fit
degree = 1

if not os.path.isfile(filename):
  print "ERROR:", filename, "does not exist"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read, parse and fit data
# Format: gc^2 beta_s stat extrap interp optimize tot
# gc^2 is always incremented in steps of 0.01
# Take a first pass through to identify the IRFP
# Stick it in between the gc^2 where the sign of beta_s changes
old_beta = 1.0
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  max_gSq = float(temp[0])
  new_beta = float(temp[1])
  if old_beta > 0.0 and new_beta < 0.0:
    IRFP = max_gSq - 0.005
  old_beta = new_beta

# Check that requested fit range remains within the available data
if IRFP + halfrange > max_gSq:
  print "Requested fit range overflows the available data"
  print "The IRFP is at", IRFP, "while max_gSq =", max_gSq
  sys.exit(0)

# Now accumulate data in given range
gSqList = []
datList = []
wList = []
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  gc = float(temp[0])
  if gc > IRFP - halfrange and gc < IRFP + halfrange:
    gSqList.append(gc)
    datList.append(float(temp[1]))
    wList.append(np.power(float(temp[6]), -1))    # Squared in fit...
gSq = np.array(gSqList)
dat = np.array(datList)
weight = np.array(wList)

# Require dof > 0
dof = len(gSq) - degree - 1
if dof <= 0:
  print "Might as well expand the fit range to get some dof..."
  sys.exit(0)

print "Fitting %d points within %.4g of the IRFP at %.4g," % \
      (len(gSq), halfrange, IRFP),
print "from %.4g to %.4g" % (gSq[0], gSq[-1])

# Return polynomial coefficient values and covariance matrix
# but not diagnostic information and singular value decomposition
# Intercept is last polynomial coefficient
out, cov = np.polyfit(gSq, dat, degree, full=False, w=weight, cov=True)

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
chi = (dat - fit(gSq)) * weight
chiSq = (chi**2).sum()
chiSq_dof = chiSq / float(dof)
cov /= chiSq_dof
CL = 1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq)
slope = out[0]
err = np.sqrt(cov[0][0])          # Component of covar matrix
print "0 %.6g %.4g # %.4g/%d = %.4g --> %.4g," \
      % (slope, err, chiSq, dof, chiSq / dof, CL),

# Convert to gamma_g^* = log(1 + slope) / log(s)
# Print positive value for convenience
gamma = np.log(1.0 + 2.0 * slope * np.log(s)) / np.log(s)
print "gamma_g^* = %.4g" % (-1.0 * gamma)
# ------------------------------------------------------------------
