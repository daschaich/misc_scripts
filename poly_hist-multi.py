#!/usr/bin/python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------
# Plot histogram of Polyakov loop magnitudes from dygraph data files
# Superimpose multiple ensembles specified by input arguments
# Combine both high and low starts for aspect ratios alpha=2, 3 and 4
# Save resulting plot as ./[W]poly_hist.pdf

# Parse arguments: First whether we're considering poly or Wpoly
# Then Nt, beta, mass to specify the ensembles to consider,
# the six corresponding thermalization cuts, high then low for each L,
# following by the title for the plot
# the upper bound for the y-axis,
# and finally a tag for the title for the plot
if len(sys.argv) < 12:
  print "Usage:", str(sys.argv[0]), "<file> <Nt> <beta> <mass>",
  print "<a2 high cut> <a2 low cut>",
  print "<a3 high cut> <a3 low cut>",
  print "<a4 high cut> <a4 low cut>",
  print "<y-axis upper bound> <plot title tag>",
  sys.exit(1)

# This should be either poly_mod or Wpoly_mod
poly_file = str(sys.argv[1])

Nt = int(sys.argv[2])
beta = str(sys.argv[3])
mass = str(sys.argv[4])
high_cuts = [int(sys.argv[5]), int(sys.argv[7]), int(sys.argv[9])]
low_cuts = [int(sys.argv[6]), int(sys.argv[8]), int(sys.argv[10])]

ymax = float(sys.argv[11])
title = str(sys.argv[12])

# Set up directories --- for now hard-code Nt=6, 8, three aspect ratios
Ndirs = 3
if Nt == 6:
  L = ['L=12', 'L=18', 'L=24']      # Labels
  tag = beta + '_high_m' + mass
  high_dirs = ['12nt6/b' + tag, '18nt6/b' + tag, '24nt6/b' + tag]
  tag = beta + '_low_m' + mass
  low_dirs = ['12nt6/b' + tag, '18nt6/b' + tag, '24nt6/b' + tag]
elif Nt == 8:
  L = ['L=16', 'L=24', 'L=32']      # Labels
  tag = beta + '_high_m' + mass
  high_dirs = ['16nt8/b' + tag, '24nt8/b' + tag, '32nt8/b' + tag]
  tag = beta + '_low_m' + mass
  low_dirs = ['16nt8/b' + tag, '24nt8/b' + tag, '32nt8/b' + tag]
else:
  print "ERROR: Only Nt=6 and 8 set up at the moment"
  sys.exit(1)

# Check that all data files exist and Nc=4
for dirname in high_dirs + low_dirs:
  tocheck = dirname + '/data/' + poly_file + '.csv'
  if not os.path.isfile(tocheck):
    print "ERROR:", tocheck, "not found"
    sys.exit(1)

cwd = os.getcwd()
if 'SU4' in cwd:
  Nc = 4.0
else:
  print "ERROR: Only expecting SU4 at the moment"
  sys.exit(1)

# Prepare data to plot, combining high and low starts
# Use c=0.5 if Wflowed
# Normalize both poly and Wpoly to be out of Nc rather than 1
# Keep track of maximum to set horizontal axis size
# poly format: MDTU,poly_mod              (normalized to 1)
# Wpoly format: MDTU,c=0.2,0.3,0.4,0.5    (normalized to Nc)
xmax = 0.0
dat = [[] for x in range(Ndirs)]
for i in range(Ndirs):
  # High start
  toOpen = high_dirs[i] + '/data/' + poly_file + '.csv'
  for line in open(toOpen):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= high_cuts[i]:
      continue
    if poly_file.startswith('W'):
      dat[i].append(float(temp[4]))
    else:
      dat[i].append(Nc * float(temp[1]))

  # Low start
  toOpen = low_dirs[i] + '/data/' + poly_file + '.csv'
  for line in open(toOpen):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= low_cuts[i]:
      continue
    if poly_file.startswith('W'):
      dat[i].append(float(temp[4]))
    else:
      dat[i].append(Nc * float(temp[1]))

  # TODO: A hack to let the script run before ensembles exist
  if len(dat[i]) == 0:
    dat[i].append(-10.0)

  if max(dat[i]) > xmax:
    xmax = max(dat[i])

  # Print number of measurements to allow offline checks
  print "%d measurements for %s" % (len(dat[i]), L[i])

# Create histogram
nbins = 20
plt.figure(figsize=(6.40, 3.84))    # Gnuplot default
plt.hist(dat[0], nbins, log=False, density=True, align='mid',
         edgecolor='blue', label=L[0], histtype='step', hatch='//')
plt.hist(dat[1], nbins, log=False, density=True, align='mid',
         edgecolor='green', label=L[1], histtype='step', hatch='\\\\')
plt.hist(dat[2], nbins, log=False, density=True, align='mid',
         edgecolor='red', label=L[2], histtype='step', hatch='oo')

#xmax = 1.1 * xmax   # A bit of padding at the right edge of the plot
plt.axis([0, xmax, 0.0, ymax])
plt.grid(False)

if poly_file.startswith('W'):
  plt.title('Wflowed Ploop mag, ' + title)
  plt.xlabel('|PL_W|')
else:
  plt.title('Polyakov loop magnitude, ' + title)
  plt.xlabel('|PL|')
plt.ylabel('Relative frequency')
plt.legend()

# Save a pdf
outfile = 'poly_hist.pdf'
if poly_file.startswith('W'):
  outfile = 'W' + outfile
plt.savefig(outfile, bbox_inches='tight')   # Reduce surrounding whitespace
# ------------------------------------------------------------------
