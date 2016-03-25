#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from scipy import special
# ------------------------------------------------------------------
# Fit beta function as function of gc^2 ("gSq"),
#   beta / gSq^2 = c0 + c1*gSq + c2*gSq^2 + ...
# Data files contain c=0.2, 0.25, 0.3 and 0.35

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
# Carry out fits
L = np.array([20, 16, 12], dtype = np.int)
sL = np.array([30, 24, 18], dtype = np.int)
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

  bList = []        # beta_F
  xList = []        # gc^2(L)
  datListV1 = []    # gc^2(L)
  errListV1 = []    # delta gc^2(L)
  datListV2 = []    # gc^2(sL)
  errListV2 = []    # delta gc^2(sL)
  for line in open(tofitV1):
    if line.startswith('#') or line.startswith('!'):
      continue
    temp = line.split()
    bList.append(float(temp[0]))
    xList.append(float(temp[c_index]))        # gc^2(L)
    datListV1.append(float(temp[c_index]))    # gc^2(L)
    errListV1.append(float(temp[err_index]))

  for line in open(tofitV2):
    if line.startswith('#') or line.startswith('!'):
      continue
    temp = line.split()
    datListV2.append(float(temp[c_index]))    # gc^2(sL)
    errListV2.append(float(temp[err_index]))

  beta = np.array(bList)
  x = np.array(xList)
  datV1 = np.array(datListV1)
  errV1 = np.array(errListV1)
  datV2 = np.array(datListV2)
  errV2 = np.array(errListV2)
  dat = (datV2 - datV1) / (x**2)    # L --> sL beta function
  dat /= div
  err = np.sqrt(errV2**2 + errV1**2) / (x**2)
  weight = div / err                # Will be squared in polyfit

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

  # Print fit parameters
  print "f%d(x) = %.6g * x**2" % (L[i], out[N]),
  for j in range(1, N):
    print "+ %.6g * x**%d" % (out[N - j], j + 2),
  print "+ %.6g * x**%d" % (out[0], N + 2)

  # Print error function sqrt(derivs.cov.derivs)
  # with derivs = {x^{N+2}, ..., x^2}
  terms = np.zeros(2 * N + 1)
  for j in range(N + 1):
    for k in range(N + 1):
      terms[j + k] += cov[j][k]
  print "err%d(x) = sqrt(%.8g * x**4" % (L[i], terms[2 * N]),
  for j in range(1, 2 * N):
    print "+ %.8g * x**%d" % (terms[2 * N - j], j + 4),
  print "+ %.8g * x**%d)" % (terms[0], 2 * N + 4)

  # Try printing numerical values for fit and errors
  for gSq in np.arange(0, np.amax(x), 0.01):    # Preserve uniform spacing
    if gSq < np.amin(x):
      continue
    deriv = np.empty(N + 1, dtype = np.float)
    for j in range(N + 1):
      deriv[j] = np.power(gSq, N + 2 - j)
    covar = np.array(cov, dtype = np.float)
    err = np.sqrt(np.dot(deriv, np.dot(covar, deriv)))
#    print "FIT%d %.4g %.6g %.4g" % (L[i], gSq, (gSq**2) * fit(gSq), err)

  # Print data
  temp = len(x)
  for k in range(temp):
    j = temp - 1 - k
    print "%.4g %.6g %.4g" % (x[j], dat[j] * x[j]**2, x[j]**2 / weight[j])

  # Comment out the line below to print out contributions to chi^2
  continue
  print "Chi^2:"
  for i in range(len(x)):
    print "%.4g (%.4g) %.4g" % (beta[i], x[i], chi[i]**2)
  print "Total: %.4g %d %.4g %.4g" % (chiSq, dof, chiSq_dof, CL)
# ------------------------------------------------------------------
