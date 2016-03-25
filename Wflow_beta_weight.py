#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from scipy import special
# ------------------------------------------------------------------
# Compute s=1.5 step scaling function for given c and tau
# More directly fit beta function as function of gc^2 ("gSq"),
#   beta / gSq^2 = c0 + c1*gSq + c2*gSq^2 + ...
# Data files contain c=0.2, 0.25, 0.3 and 0.35
# We print (a/L)^2 SSF err, (a/L)^2-->0 linear extrapolation, and its slope

# This version scales the 12-->18 error bars by the chi^2 from the
# original continuum extrapolations, if chi^2>1

# Parse arguments: the number of smearing steps, the Wilson flow parameter c,
# the t-shift parameter tau (which will tell us what files to use),
# and the degree of the polynomial interpolation
# The data files are produced by average_Wflow.py, which already includes
# the tree-level perturbative finite-volume + zero-mode corrections
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<#HYP> <c> <tau> <degree>"
  sys.exit(1)
HYP_switch = int(sys.argv[1])
c_tag = str(sys.argv[2]).rstrip('0')    # Strip trailing zeroes
tau = str(sys.argv[3])    # Need to save as string for file formatting...
N = int(sys.argv[4])

# Swap between once- and twice-smeared runs
if HYP_switch == 1:
  dirpat = "Run_APBC8_"
elif HYP_switch == 2:
  dirpat = "Run_2HYPAPBC8_"
else:
  print "Error: don't have %d-times smeared data" % HYP_switch
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
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Carry out fits and store results
L = np.array([20, 16, 12], dtype = np.int)
sL = np.array([30, 24, 18], dtype = np.int)
u_min = -1
u_max = -1
params = []   # Will be list of vectors
covars = []   # Will be list of matrices
div = np.log(9.0 / 4.0)
for i in range(len(L)):
  tofitV1 = dirpat + str(L[i]) + str(L[i]) \
                   + '/results/Wflow-all-tau' + tau + '.dat'

  tofitV2 = dirpat + str(sL[i]) + str(sL[i]) \
                   + '/results/Wflow-all-tau' + tau + '.dat'

  # Make sure averaged data exist to be analyzed
  if not os.path.isfile(tofitV1):
    print "ERROR:", tofitV1, "not found"
    sys.exit(1)

  # Make sure averaged data exist to be analyzed
  if not os.path.isfile(tofitV2):
    print "ERROR:", tofitV2, "not found"
    sys.exit(1)

  xList = []        # gc^2(L)
  datListV1 = []    # gc^2(L)
  errListV1 = []    # delta gc^2(L)
  datListV2 = []    # gc^2(sL)
  errListV2 = []    # delta gc^2(sL)
  for line in open(tofitV1):
    if line.startswith('#') or line.startswith('!'):
      continue
    temp = line.split()
    xList.append(float(temp[c_index]))        # gc^2(L)
    datListV1.append(float(temp[c_index]))    # gc^2(L)
    errListV1.append(float(temp[err_index]))

  for line in open(tofitV2):
    if line.startswith('#') or line.startswith('!'):
      continue
    temp = line.split()
    datListV2.append(float(temp[c_index]))    # gc^2(sL)
    errListV2.append(float(temp[err_index]))

  x = np.array(xList)
  datV1 = np.array(datListV1)
  errV1 = np.array(errListV1)
  datV2 = np.array(datListV2)
  errV2 = np.array(errListV2)
  dat = (datV2 - datV1) / (x**2)    # L --> sL beta function
  dat /= div
  err = np.sqrt(errV2**2 + errV1**2) / (x**2)
  weight = div / err                # Will be squared in polyfit

  # Record available range of input gc^2=u, from min(L=30) to max(L=12)
  if L[i] == 12:
    u_max = datV1.max()
  if sL[i] == 30:
    u_min = datV2.min()

  # Fit and print out confidence level
  out, cov = np.polyfit(x, dat, N, full=False, w=weight, cov=True)
  fit = np.poly1d(out)
  chi = (dat - fit(x)) * weight
  chiSq = (chi**2).sum()
  dof = len(x) - N - 1
  chiSq_dof = chiSq / float(dof)
  cov /= chiSq_dof
  CL = 1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq)
  print "# L=%d %.4g %d %.4g %.4g" % (L[i], chiSq, dof, chiSq_dof, CL)

  # Save fit parameters and covariance matrix
  params.append(out)
  covars.append(cov)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Compute step scaling function and its uncertainty at this beta_F
# Linear extrapolation to (a / L)^2 --> 0
# Recall s=1.5 is hard-coded
print "# Fitting for %.4g <= u <= %.4g" % (u_min, u_max)
for gSq in np.arange(0, u_max, 0.01):    # Preserve uniform spacing
  if gSq < u_min:
    continue
  xList = []
  datList = []
  wList = []
  for i in range(len(L)):
    # This is the step-scaling function / discrete beta function
    xList.append(1.0 / L[i]**2)
    fit = np.poly1d(params[i]) #defining the fit function
    datList.append((gSq**2) * fit(gSq))

    # Now we need uncertainties
    # Error propagation involves the derivatives of func w.r.t. each p
    N = len(params[i]) - 1
    deriv = np.empty(N + 1, dtype = np.float)
    for j in range(N + 1):
      deriv[j] = np.power(gSq, N + 2 - j)
    covar = np.array(covars[i], dtype = np.float)
    u_err = np.sqrt(np.dot(deriv, np.dot(covar, deriv)))
    wList.append(1.0 / u_err)     # Will be squared in polyfit
  #  print "L=%d: %.6g %.4g" % (L[i], u, u_err)

  # Now for linear (a / L)^2 --> 0 extrapolation
  x = np.array(xList)
  dat = np.array(datList)
  weight = np.array(wList)

  dof = len(x) - 2
  if dof == 0:
    print "Error: Let's require dof>0 for weighted extrapolation"
    sys.exit(1)

  # Print intercept from linear extrapolation, the last polynomial coefficient
  out, cov = np.polyfit(x, dat, 1, full=False, w=weight, cov=True)
  fit = np.poly1d(out)
  chi = (dat - fit(x)) * weight
  chiSq = (chi**2).sum()

  # !!! This is where the magic happens
  save = weight[-1]
  if chiSq > 1:
    weight[-1] /= chiSq

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
  weight[-1] = save
  for i in range(len(x)):
    print "%.4g %.4g %.6g %.4g" % (x[i], gSq, dat[i], 1.0 / weight[i])

  # It may also be useful to monitor the slope when optimizing tau
  slope = out[0]
  err = np.sqrt(cov[0][0])
  print "Slope %.4g %.4g %.2g" % (gSq, slope, err)
# ------------------------------------------------------------------
