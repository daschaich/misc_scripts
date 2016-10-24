#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from scipy import optimize
from scipy import special
# ------------------------------------------------------------------
# Scan over all beta_F, extract g_c^2(L) from rational function interpolations
#   1 / g_c^2 = (b / 12) * (c2 + c3*b + c4*b^2) / (1 + c0*b + c1*b^2)
#   --> g_c^2 = (1 + c0*b + c1*b^2) / (c2*b + c3*b^2 + c4*b^3)
# and fit to constant plus power in L

# Parse arguments: the minimum L, the Wilson flow parameter c,
# and the t-shift parameter tau (which will tell us what files to use)
# Fourth argument is g_*^2
# Optional fifth argument tells us to use the plaquette observable
# The data files are produced by average_Wflow.py, which already includes
# the perturbative finite-volume + zero-mode corrections
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<Lmin> <c> <tau> <gstar^2> <obs>"
  sys.exit(1)
Lmin = int(sys.argv[1])
c_tag = str(sys.argv[2]).rstrip('0')    # Strip trailing zeroes
tau = str(sys.argv[3])    # Need to save as string for file formatting...
gstar = float(sys.argv[4])
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

# Set L based on input L_min
Llist = []
for i in [12, 16, 18, 20, 24, 30, 32, 36]:
  if i < Lmin:
    continue
  Llist.append(i)
L = np.array(Llist, dtype = np.int)
if len(L) < 3:
  print "Error: d.o.f. <= 0..."
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

# Interpolation
# errfunc will be minimized via least-squares optimization
# Small p_in seem to help the greedy algorithm find the right minimum
func = lambda p, x: (1.0 + x * (p[0] + x * p[1])) \
                  / (x * (p[2] + x * (p[3] + x * p[4])))
p_in = np.array([0.01, 0.01, 0.01, 0.01, 0.01])
errfunc = lambda p, x, y, err: (func(p, x) - y) / err

# Define corresponding Jacobian matrix
def jac(p, x, y, err):
  J = np.empty((x.size, p.size), dtype = np.float)
  num = 1.0 + p[0] * x + p[1] * x**2
  den = p[2] * x + p[3] * x**2 + p[4] * x**3
  J[:, 0] = x / den
  J[:, 1] = x**2 / den
  J[:, 2] = -x * num / den**2
  J[:, 3] = -x**2 * num / den**2
  J[:, 4] = -x**3 * num / den**2
  for i in range(p.size):
    J[:, i] /= err
  return J

# Derivatives for error propagation are similar
def get_derivs(p):
  deriv = np.empty(p.size, dtype = np.float)
  num = 1 + p[0] * beta + p[1] * beta**2
  den = p[2] * beta + p[3] * beta**2 + p[4] * beta**3
  deriv[0] = beta / den
  deriv[1] = beta**2 / den
  deriv[2] = -beta * num / den**2
  deriv[3] = -beta**2 * num / den**2
  deriv[4] = -beta**3 * num / den**2
  return deriv

# Power-law fit
powfunc = lambda p, x: p[0] * np.power(x, p[1])
errpow = lambda p, x, y, err: (powfunc(p, x) - y) / err
pow_in = np.array([1.0, 1.0], dtype = np.float)

# Corresponding Jacobian matrix
def powjac(p, x, y, err):
  J = np.empty((x.size, p.size))
  J[:, 0] = np.power(x, p[1]) / err
  J[:, 1] = p[0] * np.log(x) * np.power(x, p[1]) / err
  return J
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
  if abs(out[0]) > 1:
    print "WARNING: Fit may be unstable..."
  params.append(out)
  covars.append(pcov)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Fit (g^2(L) - gstar^2) = A*L**gamma
# with g^2(L) from interpolations
for beta in np.arange(3.75, 6.0, 0.01):
  dat = np.zeros(len(L), dtype = np.float)
  err = np.zeros_like(dat)
  for i in range(len(L)):
    # Get g^2(L) - gstar^2 for each L and the propagated error
    # Error propagation involves the derivatives of func w.r.t. each p
    dat[i] = np.abs(float(func(params[i], beta)) - gstar)
    deriv = get_derivs(np.array(params[i]))
    covar = np.array(covars[i], dtype = np.float)
    err[i] = np.sqrt(np.dot(deriv, np.dot(covar, deriv)))
#    print "L=%d: %.6g %.4g" % (L[i], dat[i], err[i])

  # Now fit (g^2(L) - gstar^2) = A*L**gamma

  all_out = optimize.least_squares(errpow, pow_in,
                                   jac=powjac, method='lm', args=(L, dat, err))
  out = all_out.x
  tj = all_out.jac
  pcov = np.linalg.inv(np.dot(np.transpose(tj), tj))

  # Compute chiSq and print exponent -- the last component of the output
  chiSq = ((errpow(out, L, dat, err))**2).sum()
  CL = 1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq)
  print "%.4g %.6g %.4g # %.4g/%d = %.4g --> %.4g" \
        % (beta, out[-1], np.sqrt(pcov[-1][-1]), chiSq, dof, chiSq / dof, CL)
# ------------------------------------------------------------------
