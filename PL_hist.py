#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of (c=0.5 Wilson-flowed) Polyakov loop magnitudes
# Save resulting plot as ./WL_hist_$tag.pdf
# Extract $tag from path rather than input argument

# Parse arguments: First is thermalization cut
# Second is how many bins to have in the histogram
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<cut> <bins>"
  sys.exit(1)
cut = int(sys.argv[1])
nbins = int(sys.argv[2])

# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

# Extract tag from path as everything after the last '/'
cwd = os.getcwd()
tag = (cwd.split('/'))[-1]

# Extract c=0.5 modulus as last (fourth) datum on each line of data file
# (Format: MDTU,c=0.2,0.3,0.4,0.5)
dat = []
for line in open('data/WPoly_mod.csv'):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= cut:
    continue
  dat.append(float(temp[3]))

# Print number of measurements to allow offline checks
print "%d measurements for %s" % (len(dat), tag)

# Create histogram
# TODO: The '4.0' below should really be Nc...
plt.figure(figsize=(6.40, 3.84))    # Gnuplot default
plt.hist(dat, bins=np.arange(0.0, 4.0, 0.05),
         log=False, normed=False, align='mid',
         edgecolor='blue', label=tag, histtype='step', hatch='//')

plt.xlim([0.0, 1.0])
plt.grid(False)

plt.xlabel('|PL_W|')
plt.ylabel('Count')
plt.legend()

# Save a pdf
outfile = 'WPoly_hist_' + tag + '.pdf'
plt.savefig(outfile, bbox_inches='tight')   # Reduce surrounding whitespace
# ------------------------------------------------------------------
