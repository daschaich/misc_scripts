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

# TODO: Should be able to determine Nc automatically from path...
Nc = 4.0
bin_width = Nc / float(sys.argv[2])

# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

# Extract tag from path as everything after the last '/'
cwd = os.getcwd()
tag = (cwd.split('/'))[-1]

# Extract c=0.5 modulus as last (fourth) datum on each line of data file
# Format: MDTU,c=0.2,0.3,0.4,0.5
dat = []
for line in open('data/Wpoly_mod.csv'):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= cut:
    continue
  dat.append(float(temp[4]))

# Print number of measurements to allow offline checks
print "%d measurements for %s" % (len(dat), tag)

# Create histogram
plt.figure(figsize=(6.40, 3.84))    # Gnuplot default
plt.hist(dat, bins=np.arange(0.0, Nc, bin_width),
         log=False, density=False, align='mid',
         edgecolor='blue', label=tag, histtype='step', hatch='//')

plt.xlim([0.0, Nc])
plt.grid(False)

plt.xlabel('|PL_W|')
plt.ylabel('Count')
plt.legend()

# Save a pdf
outfile = 'Wpoly_hist_' + tag + '.pdf'
plt.savefig(outfile, bbox_inches='tight')   # Reduce surrounding whitespace
# ------------------------------------------------------------------
