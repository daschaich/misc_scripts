#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from scipy import optimize
from scipy import special
# ------------------------------------------------------------------
# Compute s=1.5 step scaling function for given g^2 and given tau
# More directly fit Wilson flow coupling as function of beta_F ("b"),
# For now we have only (2, 2) and (2, 3) rational functions for interpolation
#   1 / g_c^2 = (b / 12) * (c2 + c3*b + c4*b^2) / (1 + c0*b + c1*b^2)
#   --> g_c^2 = (1 + c0*b + c1*b^2) / (c2*b + c3*b^2 + c4*b^3)
# or
#   --> g_c^2 = (1 + c0*b + c1*b^2) / (c2*b + c3*b^2 + c4*b^3 + c5*b^4)
# Data files contain c=0.2, 0.25, 0.3 and 0.35
# We print (a/L)^2 SSF err, (a/L)^2-->0 linear extrapolation, and its slope

# Parse arguments: the scale factor s, the Wilson flow parameter c,
# the t-shift parameter tau (which will tell us what files to use),
# and the form of the rational function
# Optional fifth  argument tells us to use the plaquette observable
# The data files are produced by average_Wflow.py, which already includes
# the perturbative finite-volume + zero-mode corrections
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<s> <c> <tau> <fit_form> <obs>"
  sys.exit(1)
s_tag = str(sys.argv[1])
c_tag = str(sys.argv[2]).rstrip('0')    # Strip trailing zeroes
tau = str(sys.argv[3])    # Need to save as string for file formatting...
fit_form = int(sys.argv[4])
dirpat = "Run_APBC12_"

# Choose which observable to use -- require 'plaq' as specific argument
if len(sys.argv) > 5:
  if str(sys.argv[5]) == 'plaq':
    filetag = '/results/Wplaq-all-tau'
  else:
    print "Warning: Implicitly using clover observable"
    filetag = '/results/Wflow-all-tau'
else:
  filetag = '/results/Wflow-all-tau'

# Set L based on s_tag read in
# !!! Note order of L: decreasing small then large with fixed s=1.5
if s_tag == "3/2"
  #L = np.array([24, 20, 16, 12, 36, 30, 24, 18], dtype = np.int)
  L = np.array([24, 20, 16, 36, 30, 24], dtype = np.int)
elif s_tag == "2"
  L = np.array([18, 16, 12, 36, 32, 24], dtype = np.int)
elif s_tag == "4/3"
  L = np.array([24, 18, 12, 32, 24, 16], dtype = np.int)
else:
  print "Error: only s=3/2, 2 and 4/3 set up,",
  print "while", c_tag, "was read in"
  sys.exit(1)

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

# errfunc will be minimized via least-squares optimization
# Small p_in seem to help the greedy algorithm find the right minimum
if fit_form == 22:
  func = lambda p, x: (1.0 + x * (p[0] + x * p[1])) \
                    / (x * (p[2] + x * (p[3] + x * p[4])))
  p_in = [0.1, 0.1, 0.1, 0.1, 0.1]
elif fit_form == 23:
  func = lambda p, x: (1.0 + x * (p[0] + x * p[1])) \
                    / (x * (p[2] + x * (p[3] + x * (p[4] + x * p[5]))))
  p_in = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
  if c_tag == "0.2" or c_tag == "0.25":
    p_in = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
else:
  print "Error: only (2, 2) and (2, 3) rational functions set up,",
  print "while", str(fit_form), "was read in"
  sys.exit(1)
errfunc = lambda p, x, y, err: (func(p, x) - y) / err
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Carry out fits and store results
u_min = -1
u_max = -1
all_beta = []   # Will be list of all beta on each volume
params = []     # Will be list of vectors
covars = []     # Will be list of matrices
for i in range(len(L)):
  tofit = dirpat + str(L[i]) + str(L[i]) + filetag + tau + '.dat'

  # Make sure averaged data exist to be analyzed
  if not os.path.isfile(tofit):
    print "ERROR:", tofit, "not found"
    sys.exit(1)

  xList = []
  datList = []
  errList = []
  for line in open(tofit):
    if line.startswith('#') or line.startswith('!'):
      continue
    temp = line.split()
    xList.append(float(temp[0]))
    datList.append(float(temp[c_index]))
    errList.append(float(temp[err_index]))

  x = np.array(xList)
  all_beta.append(x)
  dat = np.array(datList)
  err = np.array(errList)

  # Require enough degrees of freedom
  dof = len(x) - len(p_in)
  if dof <= 0:
    print "ERROR: Not enough data points to fit to rational function"
    sys.exit(1)

  # Record available range of input gc^2=u, from min(L_min) to max(L_max)
  if L[i] == min(L):
    u_max = dat.max()
  elif L[i] == max(L):
    u_min = dat.min()

  # Save fit parameters and covariance matrix
  out, pcov, infodict, errmsg, success = \
                    optimize.leastsq(errfunc, p_in[:], args=(x, dat, err), \
                                     full_output=1)

  if success < 0 or success > 4:
    print "ERROR: L=%d fit failed with the following error message:" % L[i]
    print errmsg
    sys.exit(1)
  params.append(out)
  covars.append(pcov)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Compute step scaling function and its uncertainty at this beta_F
# Linear extrapolation to (a / L)^2 --> 0
# Recall s=1.5 is hard-coded
div = np.log(9.0 / 4.0)
print "# Fitting for %.4g <= u <= %.4g" % (u_min, u_max)
for gSq in np.arange(0, u_max, 0.01):    # Preserve uniform spacing
  if gSq < u_min:
    continue
  xList = []
  datList = []
  wList = []
  for i in range(len(L) / 2):
    I = i + len(L) / 2    # Index of large volume L --> s*L
    # Second argument is initial guess for where the zero might be
    beta = optimize.fsolve(lambda x : func(params[i], x) - gSq, 6.0)

    # Make sure beta is within available range on large volume
    # When we have an IRFP, we might as well just end the file here,
    # so let's comment out the error message so it doesn't break anything
    if beta < all_beta[I].min() or beta > all_beta[I].max():
      print "# ERROR: beta = %.4g falls outside [%.4g, %.4g] for L=%d" \
            % (beta, all_beta[I].min(), all_beta[I].max(), L[I])
      sys.exit(1)

    # This is the step-scaling function / discrete beta function
    # func returns a list for some reason -- cast it to a float
    u = float(func(params[i], beta))
    Sigma = float(func(params[I], beta))
    SSF = (Sigma - u) / div

    # Now we need uncertainties
    # Error propagation involves the derivatives of func w.r.t. each p
    deriv = np.empty(len(params[i]), dtype = np.float)

    # First get error on small volume
    p = np.array(params[i])
    num = 1 + p[0] * beta + p[1] * beta**2
    if fit_form == 22:
      den = p[2] * beta + p[3] * beta**2 + p[4] * beta**3
    elif fit_form == 23:
      den = p[2] * beta + p[3] * beta**2 + p[4] * beta**3 + p[5] * beta**4
    deriv[0] = beta / den
    deriv[1] = beta**2 / den
    deriv[2] = -beta * num / den**2
    deriv[3] = -beta**2 * num / den**2
    deriv[4] = -beta**3 * num / den**2
    if fit_form == 23:
      deriv[5] = -beta**4 * num / den**2
    covar = np.array(covars[i], dtype = np.float)
    u_err = np.sqrt(np.dot(deriv, np.dot(covar, deriv)))
  #  print "L=%d: %.6g %.4g" % (L[i], u, u_err)

    # Now get error on large volume
    p = np.array(params[I])
    num = 1 + p[0] * beta + p[1] * beta**2
    if fit_form == 22:
      den = p[2] * beta + p[3] * beta**2 + p[4] * beta**3
    elif fit_form == 23:
      den = p[2] * beta + p[3] * beta**2 + p[4] * beta**3 + p[5] * beta**4
    deriv[0] = beta / den
    deriv[1] = beta**2 / den
    deriv[2] = -beta * num / den**2
    deriv[3] = -beta**2 * num / den**2
    deriv[4] = -beta**3 * num / den**2
    if fit_form == 23:
      deriv[5] = -beta**4 * num / den**2
    covar = np.array(covars[I], dtype = np.float)
    Sigma_err = np.sqrt(np.dot(deriv, np.dot(covar, deriv)))
  #  print "L=%d: %.6g %.4g" % (L[I], Sigma, Sigma_err)

    # And just add them in quadrature
    err = np.sqrt(Sigma_err**2 + u_err**2)
  #  print "%.4g %.6g %.4g" % (1.0 / float(L[i])**2, SSF, err)
    xList.append(1.0 / float(L[i])**2)
    datList.append(SSF)
    wList.append(np.power(err, -1))    # Squared in fit...

  # Now for linear (a / L)^2 --> 0 extrapolation
  x = np.array(xList)
  dat = np.array(datList)
  weight = np.array(wList)

  # Handle case of no degrees of freedom
  dof = len(x) - 2
  if dof == 0:
    # Central values
    out = np.polyfit(x, dat, 1, full=False, w=weight, cov=False)
    intercept = out[1]
    # It may also be useful to monitor the slope when optimizing tau
    print "Slope %.4g %.4g" % (gSq, out[0])

    # Shifted by uncertainties
    tofit = dat + 1.0 / weight
    out = np.polyfit(x, tofit, 1, full=False, w=weight, cov=False)
    err = out[1] - intercept
    print "%.4g %.6g %.4g # 0 dof" % (gSq, intercept, err)
    sys.exit(0)

  # Print intercept from linear extrapolation, the last polynomial coefficient
  out, cov = np.polyfit(x, dat, 1, full=False, w=weight, cov=True)
  fit = np.poly1d(out)
  chi = (dat - fit(x)) * weight
  chiSq = (chi**2).sum()
  chiSq_dof = chiSq / float(dof)
  cov /= chiSq_dof
  intercept = out[1]
  int_err = np.sqrt(cov[1][1])
  CL = 1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq)
  print "%.4g %.6g %.4g # %.4g %.4g" % (gSq, intercept, int_err, chiSq, CL)

  # Print data points themselves vs. (a / L)^2, in increasing order
  for i in range(len(x)):
    print "%.4g %.4g %.6g %.4g" % (x[i], gSq, dat[i], 1.0 / weight[i])

  # It may also be useful to monitor the slope when optimizing tau
  slope = out[0]
  err = np.sqrt(cov[0][0])
  print "Slope %.4g %.4g %.2g" % (gSq, slope, err)
# ------------------------------------------------------------------
