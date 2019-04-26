#!/usr/bin/python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of Polyakov loop magnitudes from dygraph data files
# Show only the aspect ratio from which the script is called (for animations)
# Still combine both high and low starts for that aspect ratio
# Save resulting plot as ./[W]poly_hist_$tag.pdf
# Extract volume and Nt from path rather than input argument

# Parse arguments: First whether we're considering poly or Wpoly
# Then thermalization cuta and how many bins to have in the histogram
if len(sys.argv) < 8:
  print "Usage:", str(sys.argv[0]), "<file> <Nt> <beta> <mass>",
  print "<high cut> <low cut> <xmax> <bins>",
  sys.exit(1)

# This should be either poly_mod or Wpoly_mod
poly_file = str(sys.argv[1])

Nt = int(sys.argv[2])
beta = str(sys.argv[3])
mass = str(sys.argv[4])
high_cut = int(sys.argv[5])
low_cut = int(sys.argv[6])

xmax = float(sys.argv[3])
nbins = int(sys.argv[4])

# First make sure we're calling this from the right place
# Then extract tag from path as everything after the last '/'
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)
cwd = os.getcwd()
if 'SU4' in cwd:
  Nc = 4.0
else:
  print "ERROR: Only expecting SU4 at the moment"
  sys.exit(1)
tag = (cwd.split('/'))[-1]

# Prepare data to plot, using c=0.5 if Wflowed
# Normalize both poly and Wpoly to be out of Nc rather than 1
# poly format: MDTU,poly_mod              (normalized to 1)
# Wpoly format: MDTU,c=0.2,0.3,0.4,0.5    (normalized to Nc)
dat = []
toOpen = 'data/' + poly_file + '.csv'
for line in open(toOpen):
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

# Print number of measurements to allow offline checks
print "%d measurements for %s" % (len(dat), tag)

# Create histogram
plt.figure(figsize=(6.40, 3.84))    # Gnuplot default
plt.hist(dat, nbins, log=False, normed=False, align='mid',
         edgecolor='blue', label=tag, histtype='step', hatch='//')

plt.xlim([0.0, xmax])
plt.grid(False)

if poly_file.startswith('W'):
  plt.xlabel('|PL_W|')
else:
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
