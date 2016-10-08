#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from scipy import optimize
from scipy import special
# ------------------------------------------------------------------
# Fit Wilson flow coupling as function of beta_F ("b")
#   1 / g_c^2 = (b / 12) * (c2 + c3*b + c4*b^2) / (1 + c0*b + c1*b^2)
#   --> g_c^2 = (1 + c0*b + c1*b^2) / (c2*b + c3*b^2 + c4*b^3)
# or
#   --> g_c^2 = (1 + c0*b + c1*b^2) / (c2*b + c3*b^2 + c4*b^3 + c5*b^4)
# Data files contain c=0.2, 0.25, 0.3 and 0.35

# Parse arguments: the Wilson flow parameter c,
# the file containing the data to fit,
# and the form of the rational function
# Different t-shifts tau correspond to different input files
# The data files are produced by average_Wflow.py, which already includes
# the perturbative finite-volume + zero-mode corrections
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<c> <file> <fit_form>"
  sys.exit(1)
c_tag = str(sys.argv[1]).rstrip('0')    # Strip trailing zeroes
infile = str(sys.argv[2])
fit_form = int(sys.argv[3])

# Set c_index and err_index based on c_tag read in
if c_tag == "0.2":
  c_index = 1
  err_index = 2
elif c_tag == "0.25":
  c_index = 3
  err_index = 4
elif c_tag == "0.3":
  c_index = 5
  err_index = 6
elif c_tag == "0.35":
  c_index = 7
  err_index = 8
else:
  print "Error: only c=0.2, 0.25, 0.3 and 0.35 set up,",
  print "while", c_tag, "was read in"
  sys.exit(1)

# Make sure we're calling this from the right place
if not os.path.isfile(infile):
  print "ERROR:", infile, "not found"
  sys.exit(1)

# errfunc will be minimized via least-squares optimization
# Small p_in seem to help the greedy algorithm find the right minimum
if fit_form == 11:
  func = lambda p, x: (1.0 + x * p[0]) / (x * (p[2] + x * p[3]))
  p_in = [0.01, 0.01, 0.01, 0.01, 0.01]
elif fit_form == 22:
  func = lambda p, x: (1.0 + x * (p[0] + x * p[1])) \
                    / (x * (p[2] + x * (p[3] + x * p[4])))
  p_in = [0.01, 0.01, 0.01, 0.01, 0.01]
elif fit_form == 23:
  func = lambda p, x: (1.0 + x * (p[0] + x * p[1])) \
                    / (x * (p[2] + x * (p[3] + x * (p[4] + x * p[5]))))
  p_in = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
  if c_tag == "0.2" or c_tag == "0.25":
    p_in = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
elif fit_form == 33:
  func = lambda p, x: (1.0 + x * (p[0] + x * (p[1] + x * p[2]))) \
                    / (x * (p[3] + x * (p[4] + x * (p[5] + x * p[6]))))
  p_in = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
  if 'L24' in infile or 'L32' in infile or 'L36' in infile:
    p_in = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
elif fit_form == 4:
  func = lambda p, x: p[0] + (p[1] + (p[2] + (p[3] + p[4] / x) / x) / x) / x
  p_in = [0.01, 0.01, 0.01, 0.01, 0.01]
else:
  print "Error: only (2, 2) and (2, 3) rational functions set up,",
  print "while", str(fit_form), "was read in"
  sys.exit(1)
errfunc = lambda p, x, y, err: (func(p, x) - y) / err
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct data to fit
xList = []
datList = []
errList = []
for line in open(infile):
  if line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  xList.append(float(temp[0]))
  datList.append(float(temp[c_index]))
  errList.append(float(temp[err_index]))
x = np.array(xList)
dat = np.array(datList)
err = np.array(errList)

# Require enough degrees of freedom
dof = len(x) - len(p_in)
if dof <= 0:
  print "ERROR: Not enough data points to fit to rational function"
  sys.exit(1)

# Return fit parameters and covariance matrix
out, pcov, infodict, errmsg, success = \
                    optimize.leastsq(errfunc, p_in[:], args=(x, dat, err), \
                                     full_output=1, maxfev=10000)
if success < 0 or success > 4:
  print "WARNING: Fit failed with the following error message:"
  print errmsg
#if abs(out[0]) > 1:          # Doesn't seem to be a problem
#  print "I would prefer p[0] ~ O(-0.1)..."

if fit_form == 11:
  print "(1 + x*%.6g) /" % (out[0]),
  print "(x*(%.6g + x*%.6g))" % (out[2], out[3]),
elif fit_form == 22:
  print "(1 + x*(%.6g + x*%.6g)) /" % (out[0], out[1]),
  print "(x*(%.6g + x*(%.6g + x*%.6g)))" % (out[2], out[3], out[4]),
elif fit_form == 23:
  print "(1 + x*(%.6g + x*%.6g)) /" % (out[0], out[1]),
  print "(x*(%.6g + x*(%.6g + x*(%.6g + x*%.6g))))" \
        % (out[2], out[3], out[4], out[5]),
elif fit_form == 33:
  print "(1 + x*(%.6g + x*(%.6g + x*%.6g))) /" % (out[0], out[1], out[2]),
  print "(x*(%.6g + x*(%.6g + x*(%.6g + x*%.6g))))" \
        % (out[3], out[4], out[5], out[6]),
elif fit_form == 4:
  print "%.6g + (%.6g + (%.6g + (%.6g + %.6g / x) / x) / x) / x" \
        % (out[0], out[1], out[2], out[3], out[4]),

# Compute chiSq and confidence level of fit
# The infodict returned by leastsq includes fvec = f(x) - y
chiSq = (infodict['fvec']**2).sum()
CL = 1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq)
print "# %.4g %d --> %.4g" % (chiSq, dof, CL)

#for i in range(len(x)):
##  print x[i], dat[i], err[i], func23(out, x[i])
#  print "%.4g %.4g" % (x[i], ((infodict['fvec'])[i])**2)
#print "Total: %.4g" % chiSq
# ------------------------------------------------------------------
