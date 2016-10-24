#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from scipy import optimize
from scipy import special
# ------------------------------------------------------------------
# Compute step scaling function for given scale change (s), c and tau
# More directly fit Wilson flow coupling as function of beta_F ("b"),
# Consider either (2, 2) rational function interpolation
#   1 / g_c^2 = (b / 12) * (c2 + c3*b + c4*b^2) / (1 + c0*b + c1*b^2)
#   --> g_c^2 = (1 + c0*b + c1*b^2) / (c2*b + c3*b^2 + c4*b^3)
# or (motivated by 1503.01132)
#   1 / g_c^2 = (b / 12) * (c0 + c1/b + c2/b^2 + c3/b^3 + c4/b^4)
#   --> g_c^2 = 1 / (c0*b + c1 + c2/b + c3/b^2 + c4/b^3)
# Data files contain c=0.2, 0.25, 0.3 and 0.35
# We print (a/L)^2 SSF err, (a/L)^2-->0 linear extrapolation, and its slope

# Parse arguments: the scale factor s, the Wilson flow parameter c,
# the t-shift parameter tau (which will tell us what files to use),
# and the functional form for interpolation
# Optional fifth argument tells us to use the plaquette observable
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
# !!! Note order of L: decreasing small then large s*L
if s_tag == "3/2":
  #L = np.array([24, 20, 16, 12, 36, 30, 24, 18], dtype = np.int)
  L = np.array([24, 20, 16, 36, 30, 24], dtype = np.int)
  #L = np.array([24, 20, 36, 30], dtype = np.int)
  div = np.log(9.0 / 4.0)
elif s_tag == "2":
  L = np.array([18, 16, 12, 36, 32, 24], dtype = np.int)
  #L = np.array([18, 16, 36, 32], dtype = np.int)
  div = np.log(4.0)
elif s_tag == "4/3":
  L = np.array([24, 18, 12, 32, 24, 16], dtype = np.int)
  #L = np.array([24, 18, 32, 24], dtype = np.int)
  div = np.log(16.0 / 9.0)
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

# Derivatives for error propagation are similar
def get_derivs(p):
  deriv = np.empty(p.size, dtype = np.float)
  if fit_form == 22:
    num = 1 + p[0] * beta + p[1] * beta**2
    den = p[2] * beta + p[3] * beta**2 + p[4] * beta**3
    deriv[0] = beta / den
    deriv[1] = beta**2 / den
    deriv[2] = -beta * num / den**2
    deriv[3] = -beta**2 * num / den**2
    deriv[4] = -beta**3 * num / den**2
  elif fit_form == 5:
    den = p[0] * beta + p[1] + (p[2] + (p[3] + p[4] / beta) / beta) / beta
    deriv[0] = -beta / den**2
    deriv[1] = -1.0 / den**2
    deriv[2] = -1.0 / (beta * den**2)
    deriv[3] = -1.0 / (beta**2 * den**2)
    deriv[4] = -1.0 / (beta**3 * den**2)
  return deriv
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Carry out fits and store results
minList = []
maxList = []
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

  # Need to find max(min) and min(max)
  # to determine available range of input gc^2=u
  minList.append(dat.min())
  maxList.append(dat.max())

  # Extract and save fit parameters
  # and covariance matrix (J^T J)^{-1} where J is jacobian matrix
  all_out = optimize.least_squares(errfunc, p_in, #bounds=(-10.0, 10.0),
                                   jac=jac, method='lm', args=(x, dat, err))
  out = all_out.x
  tj = all_out.jac
  pcov = np.linalg.inv(np.dot(np.transpose(tj), tj))

  if all_out.success < 0 or all_out.success > 4:
    print "ERROR: L=%d fit failed with the following error message:" % L[i]
    print errmsg
    sys.exit(1)
  if fit_form == 22 and abs(out[0]) > 1:
    print "WARNING: Fit may be unstable..."
  params.append(out)
  covars.append(pcov)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Compute step scaling function and its uncertainty at this beta_F
# Linear extrapolation to (a / L)^2 --> 0
# Scale change s set up above
# Pad lower bound (imprecision in fsolve?) but fine to hit upper bound
u_min = max(minList) + 0.01
u_max = min(maxList)
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
    # First get error on small volume
    deriv = get_derivs(np.array(params[i]))
    covar = np.array(covars[i], dtype = np.float)
    u_err = np.sqrt(np.dot(deriv, np.dot(covar, deriv)))
#    print "L=%d: %.6g %.4g" % (L[i], u, u_err)

    # Now get error on large volume
    deriv = get_derivs(np.array(params[I]))
    covar = np.array(covars[I], dtype = np.float)
    Sigma_err = np.sqrt(np.dot(deriv, np.dot(covar, deriv)))
#    print "L=%d: %.6g %.4g" % (L[I], Sigma, Sigma_err)

    # And just add them in quadrature
    err = np.sqrt(Sigma_err**2 + u_err**2)
#    print "%.4g %.6g %.4g" % (1.0 / float(L[i])**2, SSF, err)
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
    slope = out[0]
    intercept = out[1]

    # Shifted by uncertainties
    tofit = dat + 1.0 / weight
    out = np.polyfit(x, tofit, 1, full=False, w=weight, cov=False)
    err = out[1] - intercept
    print "%.4g %.6g %.4g # 0 dof" % (gSq, intercept, err)

    # Print data points themselves vs. (a / L)^2, in increasing order
    # Include omitted volume to match up lines
    for i in range(len(x)):
      print "%.4g %.4g %.6g %.4g" % (x[i], gSq, dat[i], 1.0 / weight[i])
    if s_tag == "3/2":
      print "0.003906 %.4g NaN NaN" % gSq   # 1 / 16^2
    else:
      print "0.006944 %.4g NaN NaN" % gSq   # 1 / 12^2

    # It may also be useful to monitor the slope when optimizing tau
    print "Slope %.4g %.4g" % (gSq, slope)
    continue

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
