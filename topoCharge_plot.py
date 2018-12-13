#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
#from scipy.optimize import curve_fit
# ------------------------------------------------------------------
# Plot topological charge time series and histogram
# Fit histogram to gaussian and plot that as well

# Parse arguments: first is file to analyze, second is thermalization cut,
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<file> <cut>"
  sys.exit(1)
filename = str(sys.argv[1])
cut = int(sys.argv[2])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Extract filename, cutting ".dat" from end
tag = filename[:-4]
print tag

# Read the data
# Format: MDTU   c=0.2   c=0.3   c=0.4   c=0.5
# We're interested in c=0.5
if not os.path.isfile(filename):
  print "ERROR:", filename, "does not exist"
  sys.exit(1)
x, junk, junk, junk, y = np.loadtxt(filename, unpack=True)

# Translate thermalization cut into limit on data in file
start = -1
for i in range(0, len(y), 1):
  if x[i] >= cut:
    start = i
    break
if start == -1:     # Check that we picked up thermalized data
  print "ERROR: no data found after thermlization cut", cut
  sys.exit(1)

# Find the maximum and minimum after the thermalization cut
maxQ = int(round(max(y[start:])))
minQ = int(round(min(y[start:])))
print("maxQ = %d," % maxQ),
print("minQ = %d" % minQ)
if maxQ > -minQ:
  max_abs = maxQ
else:
  max_abs = -minQ

nbins = 2 * max_abs + 1

# Simple calculation of the susceptibility
count = 0
Qavg = 0
Q2 = 0
for i in range(start, len(y), 1):
#  print x[i], y[i]   # Sanity check
#  sys.exit(0)
  count += 1
  Q2 += y[i] * y[i]
  Qavg += y[i]

Q2 /= count
Qavg /= count
chiQ = Q2 - Qavg*Qavg

print "From direct calculation: <Q^2> = %.2f, <Q> = %.2f, Vchi_t = %.2f" \
      % (Q2, Qavg, chiQ)
#output = open(topdir+'chiQ.dat','w')
#output.write("%f %f\n" % (float(m), chiQ))
#output.close()

# Fit the data to a normal distribution -- V\chi_t is just sigma^2
(mu, sigma) = norm.fit(y[start:])
chi_t = sigma * sigma

print 'From Gaussian fit to histogram: <Q> = %.2f, Vchi_t = %.2f' % (mu, chi_t)

# Create an instance of plt and position the subplots
# (MNi) denotes the ith plot in an MxN array
fig = plt.figure()
timeseries = fig.add_subplot(211)
histogram = fig.add_subplot(223, xlim=(-max_abs - 1, max_abs + 1))
histogram_norm = fig.add_subplot(224, xlim=(-max_abs - 1, max_abs + 1))

# Make the histogram plots!
histogram.hist(y[start:], nbins, normed=False, align='mid',
               facecolor='blue', alpha=0.9)
histogram.set_xlabel('Q')
histogram.set_ylabel('#')
histogram.set_xticks(range(-max_abs - 1, max_abs + 1, 4))

n, bins, patches = histogram_norm.hist(y[start:], nbins, normed=True,
                                       align='mid', facecolor='green',
                                       alpha=0.9)
histogram_norm.set_xlabel('Q')
histogram_norm.set_ylabel('%')
#histogram_norm.set_ylabel('Percentage')
histogram_norm.set_xticks(range(-max_abs - 1, max_abs + 1, 4))  # !!! TODO: Figure out separation from max_abs; 4 is sub-optimal in some cases

# Plot the Gaussian fit
#print bins
#bins.append(16)
gauss = mlab.normpdf(bins, mu, sigma)
histogram_norm.plot(bins, gauss, 'r-', linewidth=2)

# Time series plot
timeseries.plot(x, y, 'r-')
timeseries.set_title(tag)
timeseries.set_xlabel('Trajectory')
timeseries.set_ylabel('Q')
plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

# Save a pdf
temp = (filename.split())[0]
outfile = tag + '.pdf'
plt.savefig(outfile)
# ------------------------------------------------------------------
