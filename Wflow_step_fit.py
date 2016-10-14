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
# or (motivated by 1503.01132)
#   1 / g_c^2 = (b / 12) * (c0 + c1/b + c2/b^2 + c3/b^3 + c4/b^4)
#   --> g_c^2 = 1 / (c0*b + c1 + c2/b + c3/b^2 + c4/b^3)
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
if fit_form == 22:
  func = lambda p, x: (1.0 + x * (p[0] + x * p[1])) \
                    / (x * (p[2] + x * (p[3] + x * p[4])))
  p_in = np.array([0.01, 0.01, 0.01, 0.01, 0.01])
elif fit_form == 5:
  func = lambda p, x: 1.0 / (p[0] * x + p[1] + (p[2] + (p[3] + p[4] / x)/x)/x)
  p_in = np.array([0.01, 0.01, 0.01, 0.01, 0.01])
else:
  print "Error: only (2, 2) rational function and 5th-order polynomial set up,",
  print "while", str(fit_form), "was read in"
  sys.exit(1)
errfunc = lambda p, x, y, err: (func(p, x) - y) / err

# Define corresponding Jacobian matrix
def jac(p, x, y, err):
  J = np.empty((x.size, p.size), dtype = np.float)
  if fit_form == 22:
    num = 1.0 + p[0] * x + p[1] * x**2
    den = p[2] * x + p[3] * x**2 + p[4] * x**3
    J[:, 0] = x / den
    J[:, 1] = x**2 / den
    J[:, 2] = -x * num / den**2
    J[:, 3] = -x**2 * num / den**2
    J[:, 4] = -x**3 * num / den**2
    for i in range(p.size):
      J[:, i] /= err
  elif fit_form == 5:
    den = p[0] * x + p[1] + (p[2] + (p[3] + p[4] / x) / x) / x
    J[:, 0] = -x / den**2
    J[:, 1] = -1.0 / den**2
    J[:, 2] = -1.0 / (x * den**2)
    J[:, 3] = -1.0 / (x**2 * den**2)
    J[:, 4] = -1.0 / (x**3 * den**2)
    for i in range(p.size):
      J[:, i] /= err
  return J
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

# Extract fit parameters from output
all_out = optimize.least_squares(errfunc, p_in, #bounds=(-10.0, 10.0),
                                 jac=jac, method='lm', args=(x, dat, err))
out = all_out.x

if all_out.success < 0 or all_out.success > 4:
  print "WARNING: Fit failed with the following error message:"
  print all_out.message
if fit_form == 22 and abs(out[0]) > 1:
  print "WARNING: Fit may be unstable..."

if fit_form == 22:
  print "(1 + x*(%.6g + x*%.6g)) /" % (out[0], out[1]),
  print "(x*(%.6g + x*(%.6g + x*%.6g)))" % (out[2], out[3], out[4]),
elif fit_form == 5:
  print "1 / (%.6g*x + %.6g + (%.6g + (%.6g + %.6g / x) / x) / x)" \
        % (out[0], out[1], out[2], out[3], out[4]),

# Compute chiSq and confidence level of fit
# The infodict returned by leastsq includes fvec = f(x) - y
chiSq = ((errfunc(out, x, dat, err))**2).sum()
CL = 1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq)
print "# %.4g %d --> %.4g" % (chiSq, dof, CL)

#for i in range(len(x)):
##  print x[i], dat[i], err[i], func23(out, x[i])
#  print "%.4g %.4g" % (x[i], ((infodict['fvec'])[i])**2)
#print "Total: %.4g" % chiSq
# ------------------------------------------------------------------
