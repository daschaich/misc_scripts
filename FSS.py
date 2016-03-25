#!/usr/bin/python
import glob
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
rcParams['figure.figsize'] = 8, 10    # 8x10-inch figure
# ------------------------------------------------------------------
# This script quantifies the quality of curve collapse for data
# in the given file (for fixed observable and coupling)

# Calculate Houdayer--Hartmann quality parameter for given c0, omega and ym
def calcD(ym, omega, c0):
  D = 0
  Nover = 0

  # Cycle over each point j in each data set i
  for i_str in all_L:
    i = all_L.index(i_str)
    L = float(i_str)
    for j in range(len(mf[i])):
#        print "L[%d]=%d, m[%d][%d]=%.2g" % (i, L, i, j, mf[i][j])
      xij = L * np.power(mf[i][j], 1 / ym)
      scale = 1 + c0 * np.power(mf[i][j], omega)
      yij = L * MH[i][j] / scale
      yijErr = L * err[i][j] / scale

      # Find and save bracketing points in all other data sets
      # First clear previous round
      count = 0;          cntP1 = 1
      xL[:] = 0.;         yL[:] = 0.;           wL[:] = 0.
      for ip_str in all_L:
        ip = all_L.index(ip_str)
        Lp = float(ip_str)
        if ip == i:   # Looking only at different L
          continue
        for jp in range(len(mf[ip]) - 1):
          thisX = Lp * np.power(mf[ip][jp], 1 / ym)
          nextX = Lp * np.power(mf[ip][jp + 1], 1 / ym)
          if thisX <= xij and xij <= nextX:
#              print "  L[%d]=%d, m[%d][%d]=%.2g" % (ip, Lp, ip, jp, mf[ip][jp])
#              print "  %.2g <= %.2g <= %.2g" % (thisX, xij, nextX)
            scale0 = 1 + c0 * np.power(mf[ip][jp], omega)
            scale1 = 1 + c0 * np.power(mf[ip][jp + 1], omega)
            xL[count] = thisX
            xL[cntP1] = nextX
            yL[count] = Lp * MH[ip][jp] / scale0
            yL[cntP1] = Lp * MH[ip][jp + 1] / scale1
            wL[count] = np.power(Lp * err[ip][jp] / scale0, -2)
            wL[cntP1] = np.power(Lp * err[ip][jp + 1] / scale1, -2)
            count += 2
            cntP1 += 2

            # Avoid rare case that nextX = xij and would be counted twice
            break

      # Use bracketing points to construct fiducial curve F
      # Convert lists to numpy arrays for easy multiplication
      if count > 0:
        Nover += 1     # We have overlap
        Ktot = np.sum(wL)     # Total weights
        Kx = np.sum(wL * xL)
        Ky = np.sum(wL * yL)
        Kxx = np.sum(wL * xL * xL)
        Kxy = np.sum(wL * xL * yL)
        Delta = Ktot * Kxx - Kx * Kx

        # Fiducial curve and uncertainty
        F = Kxx * Ky - Kx * Kxy + xij * (Ktot * Kxy - Kx * Ky)
        F /= Delta
        dFSq = (Kxx - 2 * xij * Kx + xij * xij * Ktot) / Delta

        # Deviation weighted by uncertainties
        D += np.power(yij - F, 2) / (dFSq + np.power(yijErr, 2))

  # Scale by number of overlapping points
  # Set large dummy deviation if there are no overlapping points
  # Hopefully that won't be a common problem
  if Nover == 0:
    print "Warning: no overlap for y_m=%.4g" % ym
    D = 1e6
  else:
    D /= Nover
#    print "D=%.4g at y_m=%.4g" % (D, ym)
  return D, Nover
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Scan over fixed range of ym and c0, check where D < min_D + 1
# and make sure that this should-be ellipse is contained within the range
def find_bounds(rect):
  ymin, ymax, c0min, c0max = rect
  all_ym = np.linspace(ymin, ymax, Nym)
  all_c0 = np.linspace(c0min, c0max, Nc0)
  D = np.zeros((Nym, Nc0), dtype = np.float)
  Nover = np.zeros((Nym, Nc0), dtype = np.int)
  for i in range(Nym):
    ym = all_ym[i]
    omega = -y0 / ym
    for j in range(Nc0):
      c0 = all_c0[j]
      D[i, j], Nover[i, j] = calcD(ym, omega, c0)

  bestD = D.min()
  coords = np.unravel_index(D.argmin(), D.shape)
  besty = all_ym[coords[0]]
  lo_y = besty
  hi_y = besty
  bestc = all_c0[coords[1]]
  lo_c = bestc
  hi_c = bestc
  plus_one = bestD + 1. / float(Nover[coords[0], coords[1]])
  for i in range(Nym):
    for j in range(Nc0):
      if D[i, j] <= plus_one:
        ym = all_ym[i]
        if ym < lo_y: lo_y = ym
        if ym > hi_y: hi_y = ym
        c0 = all_c0[j]
        if c0 < lo_c: lo_c = c0
        if c0 > hi_c: hi_c = c0

  best = [bestD, besty, bestc, Nover[coords[0], coords[1]]]
  bounds = [lo_y, hi_y, lo_c, hi_c]

  # Allow us to run with fixed c0
  c0step = (c0max - c0min) / float(Nc0)
  if c0step == 0:
    return D, rect, best, bounds

  # The new D_min may be smaller than before,
  # potentially leading to the ellipse overflowing the range
  # Try expanding the range by 5 steps in the offending direction
  # simultaneously shrinking it as much as possible on the other sides
  # (But be careful -- shrinking too much can produce an infinite loop!)
  # Hard bounds introduce another potential infinite loop,
  # dealing with which is annoying:
  lo_cTest = True
  if lo_c != c0min:                 lo_cTest = False
  if hard_c0 and lo_c == c0min_in:  lo_cTest = False

  hi_cTest = True
  if hi_c != c0max:                 hi_cTest = False
  if hard_c0 and hi_c == c0max_in:  hi_cTest = False
  if lo_y == ymin or hi_y == ymax or lo_cTest or hi_cTest:
    ystep = (ymax - ymin) / float(Nym)

    print "Warning: ellipse overflows scanned range"
    print "Re-running with the following bounds"
    print "  ymin:  %.4g -->" % ymin,
    shift_factor = 0.1 / min(ystep, c0step)
    if lo_y == ymin:  ymin -= shift_factor * ystep
    else:             ymin = max(ymin, lo_y - shift_factor * ystep)
    print "%.4g" % ymin

    print "  ymax:  %.4g -->" % ymax,
    if hi_y == ymax:  ymax += shift_factor * ystep
    else:             ymax = min(ymax, hi_y + shift_factor * ystep)
    print "%.4g" % ymax

    print "  c0min: %.4g -->" % c0min,
    if lo_c == c0min:  c0min -= shift_factor * c0step
    else:              c0min = max(c0min, lo_c - shift_factor * c0step)
    if hard_c0:        c0min = max(c0min, c0min_in)   # Force input bound
    print "%.4g" % c0min

    print "  c0max: %.4g -->" % c0max,
    if hi_c == c0max:  c0max += shift_factor * c0step
    else:              c0max = min(c0max, hi_c + shift_factor * c0step)
    if hard_c0:        c0max = min(c0max, c0max_in)   # Force input bound
    print "%.4g" % c0max

    rect = [ymin, ymax, c0min, c0max]
    D, rect, best, bounds = find_bounds(rect)

  return D, rect, best, bounds
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments: first is file to analyze,
# second is optional cut on M*L,
# third and fourth are optional hard bounds on c0 range
# For now hard-code ymin, ymax, Nym and Nc0
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file> <cut (default=0)>"
  print "       <hard bounds on c0 (min, max)>"
  sys.exit(1)
filename = str(sys.argv[1])
if len(sys.argv) > 2:
  cut = float(sys.argv[2])
  print "Cut on M*L set to %.4g" % cut
else:
  cut = 0.

if len(sys.argv) > 4:
  c0min_in = float(sys.argv[3])
  c0min = c0min_in
  c0max_in = float(sys.argv[4])
  c0max = c0max_in
  print "Bounds on c0 fixed at [%.4g, %.4g]" % (c0min, c0max)
  hard_c0 = True
else:
  c0min = -1
  c0max = 0
  hard_c0 = False

# !!! Hard-coded definitions
ymin = 1
ymax = 2
Nym = 10
if c0min == c0max:    Nc0 = 1
else:                 Nc0 = 10
y0 = -0.1     # Scaling dimension of gauge coupling
refinements = 2
refine_factor = 2

runtime = -time.time()
if not os.path.isfile(filename):
  print "ERROR:", filename, "does not exist"
  sys.exit(1)

# Only allow cut when c0 is fixed
if cut != 0. and c0min != c0max:
  print "ERROR: can't both apply cut and scan over c0"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read and parse data
# Format: mf MH err beta L Nt obs
beta = -1       # Check
obs = 'obs'     # Check
all_L = []      # System sizes (assuming constant aspect ratio)
mf = []         # Fermion mass -- list of lists
MH = []         # Hadron mass (or other observable) -- list of lists
err = []        # Error on observable -- list of lists
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  if beta == -1:
    beta = float(temp[3])
    obs = temp[6]

  # Check that line seems to be well-formed
  elif beta != float(temp[3]):
    print "ERROR: the coupling has changed from %.2g to %.2g" \
          % (beta, float(temp[3]))
    sys.exit(1)
  elif obs != temp[6]:
    print "ERROR: the observable has changed from %s to %s" \
          % (obs, temp[6])
    sys.exit(1)

  # Only load data that passes cut on M * L
  if cut > 0 and float(temp[1]) * float(temp[4]) < cut:
    continue

  # Now we should be good to read the line
  # First consider the case that this L has already been seen
  if temp[4] in all_L:
    toAdd = all_L.index(temp[4])
    (mf[toAdd]).append(float(temp[0]))
    (MH[toAdd]).append(float(temp[1]))
    (err[toAdd]).append(float(temp[2]))
  else:                         # Add new lists to mf, MH and err
    all_L.append(temp[4])       # Save as string for comparison above
    mf.append([])
    (mf[-1]).append(float(temp[0]))
    MH.append([])
    (MH[-1]).append(float(temp[1]))
    err.append([])
    (err[-1]).append(float(temp[2]))

# Check that all mf are ordered
# Die if not (need to be fancier to sort mf, MH and err jointly)
for i in range(len(all_L)):
  for j in range(1, len(mf[i])):
    prev = mf[i][j - 1]
    if mf[i][j] <= prev:
      print "ERROR: data for L=%s are not sorted:" % (all_L[i]),
      print "%.2g <= %.2g" % (mf[i][j], prev)
      sys.exit(1)

# Maximum possible number of overlapping points
maxPts = 2 * (len(all_L)  - 1)
xL = np.zeros(maxPts, dtype = np.float)
yL = np.zeros(maxPts, dtype = np.float)
wL = np.zeros(maxPts, dtype = np.float)   # Weights from errors on yL
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we have all the data, so let's analyze it
ystep = (ymax - ymin) / float(Nym)
c0step = (c0max - c0min) / float(Nc0)
print "Initial scan: %.4g <= ym <= %.4g by %.4g" % (ymin, ymax, ystep)
print "              %.4g <= c0 <= %.4g by %.4g" % (c0min, c0max, c0step)
print "              (%dx%d = %d points)" % (Nym, Nc0, Nym * Nc0)

# Allow find_bounds() to modify these variables
rect = [ymin, ymax, c0min, c0max]
D, rect, best, bounds = find_bounds(rect)
ymin, ymax, c0min, c0max = rect
bestD, besty, bestc, Nover = best
lo_y, hi_y, lo_c, hi_c = bounds

print "Minimum D = %.4g (N = %d) at ym = %.4g, c0 = %.4g" \
      % (bestD, Nover, besty, bestc)
print "Bounds: %.4g <= ym <= %.4g, %.4g <= c0 <= %.4g" \
      % (lo_y - ystep, hi_y + ystep, lo_c - c0step, hi_c + c0step)

for r in range(refinements):
  ymin = lo_y - ystep
  ymax = hi_y + ystep
  c0min = lo_c - c0step
  c0max = hi_c + c0step
  Nym *= refine_factor
  if c0step > 0:
    Nc0 *= refine_factor
  ystep = (ymax - ymin) / float(Nym)
  c0step = (c0max - c0min) / float(Nc0)
  print "\nRefinement %d: %.4g <= ym <= %.4g by %.4g" % (r + 1, ymin, ymax, ystep)
  print "              %.4g <= c0 <= %.4g by %.4g" % (c0min, c0max, c0step)
  print "              (%dx%d = %d points)" % (Nym, Nc0, Nym * Nc0)
  rect = [ymin, ymax, c0min, c0max]
  D, rect, best, bounds = find_bounds(rect)
  ymin, ymax, c0min, c0max = rect
  bestD, besty, bestc, Nover = best
  lo_y, hi_y, lo_c, hi_c = bounds
  print "Minimum D = %.4g (N = %d) at ym = %.4g, c0 = %.4g" \
        % (bestD, Nover, besty, bestc)
  print "Bounds: %.4g <= ym <= %.4g, %.4g <= c0 <= %.4g" \
        % (lo_y - ystep, hi_y + ystep, lo_c - c0step, hi_c + c0step)

# Allow us to run with fixed c0
if c0step == 0:
  runtime += time.time()
  print "Runtime: %.2g seconds" % runtime
  # To copy into gnuplot
  print "%s %.4g %.4g %.4g" \
        % (filename[-3:], besty, lo_y - ystep, hi_y + ystep)
  sys.exit(0)

# Plot final ellipse
#fig, (ellipse, collapse) = plt.subplots(2, 1, sharex=False, sharey=False)
myc = ['red', 'blue', 'green', 'black', 'purple']   # Colors
mym = ['x', 'o', 's', '.', 'D']                     # Markers
fig = plt.figure()
ellipse = fig.add_subplot(211, xlim=(ymin - ystep, ymax + ystep), \
                               ylim=(c0min - c0step, c0max + c0step))
title = 'where quality parameter $D \leq \ %.2f + 1 / %d$' % (bestD, Nover)
ellipse.set_title(r'"Ellipse" ' + title)
ellipse.set_xlabel(r'$y_m$')    # r for "raw" -- needed to parse tex
ellipse.set_ylabel(r'$c_0$')

all_ym = np.linspace(ymin, ymax, Nym)
all_c0 = np.linspace(c0min, c0max, Nc0)
ym = []
c0 = []
plus_one = bestD + 1. / float(Nover)
for i in range(Nym):
  for j in range(Nc0):
    if D[i, j] <= plus_one:
      ym.append(all_ym[i])
      c0.append(all_c0[j])
ellipse.scatter(ym, c0, c=myc[1], marker=mym[1])

# Plot final curve collapse
collapse = fig.add_subplot(212)
best = ", $c_0=$" + str("%.4g" % bestc) + ", $y_m=$" + str("%.4g" % besty)
collapse.set_title(r'Curve collapse at $y_0=$' + str(y0) + best)
collapse.set_xlabel(r'$L m^{1 / y_m}$')
collapse.set_ylabel(r'$M L$')

# Find plot range by cycling over each point j in each data set i
omega = -y0 / besty
maxX = 0
minY = 999
maxY = 0
for i in range(len(all_L)):
  L = float(all_L[i])
  tag = "$L=" + all_L[i] + "$"
  temp = [L * np.power(m, 1 / besty) for m in mf[i]]
  x = np.array(temp)
  temp = [L * MH[i][j] / (1 + bestc * np.power(mf[i][j], omega)) for j in range(len(mf[i]))]
  y = np.array(temp)
  temp = [L * err[i][j] / (1 + bestc * np.power(mf[i][j], omega)) for j in range(len(mf[i]))]
  yerr = np.array(temp)

  if np.amax(x) > maxX: maxX = np.amax(x)
  if np.amin(y) < minY: minY = np.amin(y)
  if np.amax(y) > maxY: maxY = np.amax(y)
  collapse.errorbar(x, y, yerr=yerr, c=myc[i], marker=mym[i],
                    linestyle='None', label=tag)

  # Print data for gnuplotting
#  print int(L)
#  for j in range(len(mf[i])):
#    print "%.4g %.4g %.4g" % (x[j], y[j], yerr[j])

collapse.legend(loc='best')
collapse.set_xlim(0, maxX * 1.1)
if cut == 0:
  collapse.set_ylim(minY / 1.1, maxY * 1.1)
else:
  collapse.set_ylim(cut, maxY * 1.1)    # May need refinement
plt.subplots_adjust(hspace = 0.4)
saveName = filename + ".pdf"
#plt.savefig(saveName, bbox_inches=0)    # Save without bounding-box whitespace

# Stop timing before plt.show()
# or else it will count until the plot is closed by the user
runtime += time.time()
print "Runtime: %.2g seconds" % runtime

# To copy into gnuplot
print "%s %.4g %.4g %.4g" % (filename[-3:], besty, lo_y - ystep, hi_y + ystep),
print "%.4g %.4g %.4g" % (bestc, lo_c - c0step, hi_c + c0step)
#plt.show()
# ------------------------------------------------------------------

