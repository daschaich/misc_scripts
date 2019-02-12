#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of Polyakov loop magnitudes from dygraph data files
# Show only the ensemble from which the script is called (for animations)
# Save resulting plot as ./Wpoly_hist_$tag.pdf
# Extract $tag from path rather than input argument

# Parse arguments: First whether we're considering poly or Wpoly
# Then thermalization cuta and how many bins to have in the histogram
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<file> <cut> <bins>"
  sys.exit(1)

# This should be either poly_mod or Wpoly_mod
poly_file = str(sys.argv[1])

cut = int(sys.argv[2])
xmax = float(sys.argv[3])
nbins = int(sys.argv[4])

# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

# Extract tag from path as everything after the last '/'
cwd = os.getcwd()
tag = (cwd.split('/'))[-1]

# Prepare data to plot, using c=0.5 if Wflowed
# Normalize both poly and Wpoly to be out of Nc rather than 1
# poly format: MDTU,poly_mod              (normalized to 1)
# Wpoly format: MDTU,c=0.2,0.3,0.4,0.5    (normalized to Nc)
dat = []
toOpen = 'data/' + poly_file + '.csv'
for line in open():
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= cut:
    continue
  if poly_file.startswith('W'):
    dat.append(float(temp[4]))
  else:
    dat.append(Nc * float(temp[1]))
  dat.append(float(temp[4]))

# Print number of measurements to allow offline checks
print "%d measurements for %s" % (len(dat), tag)

# Create histogram
plt.figure(figsize=(6.40, 3.84))    # Gnuplot default
plt.hist(dat, bins=np.arange(0.0, xmax, xmax / nbins),
         log=False, density=False, align='mid',
         edgecolor='blue', label=tag, histtype='step', hatch='//')

plt.xlim([0.0, xmax])
plt.grid(False)

if poly_file.startswith('W'):
  plt.title('Wflowed Polyakov loop magnitude, ' + title)
  plt.xlabel('|PL_W|')
else:
  plt.title('Polyakov loop magnitude, ' + title)
  plt.xlabel('|PL|')
plt.ylabel('Count')
plt.legend()

# Save a pdf
if poly_file.startswith('W'):
  outfile = 'Wpoly_hist_' + tag + '.pdf'
else:
  outfile = 'poly_hist_' + tag + '.pdf'
plt.savefig(outfile, bbox_inches='tight')   # Reduce surrounding whitespace
# ------------------------------------------------------------------
