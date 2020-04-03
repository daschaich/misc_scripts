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
# Then specify ensemble by beta and mass
# Next thermalization cuts, horizontal maximum (ignored if negative),
# and how many bins to have in the histogram
if len(sys.argv) < 7:
  print "Usage:", str(sys.argv[0]), "<file> <beta> <mass>",
  print "<high cut> <low cut> <xmax> <bins>",
  sys.exit(1)

# This should be either poly_mod or Wpoly_mod
poly_file = str(sys.argv[1])

beta = str(sys.argv[2])
mass = str(sys.argv[3])
high_cut = int(sys.argv[4])
low_cut = int(sys.argv[5])
xmax = float(sys.argv[6])
nbins = int(sys.argv[7])

# Set up directories, tag and plot title from path and input
# First make sure we're calling this from the right place
high_dir = 'b' + beta + '_high_m' + mass
low_dir = 'b' + beta + '_low_m' + mass
for dirname in [high_dir, low_dir]:
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

tag = 'beta_F=' + beta
temp = cwd.split('/')   # Path ends "#f/#nt#"
#title = temp[-2] + ', ' + temp[-1] + ', m=' + mass
title = 'Nf=4, ' + temp[-1] + ', m=' + mass
#title = 'pure gauge, ' + temp[-1]          # Uncomment for pure-gauge

# Prepare data to plot, combining high and low starts
# Use c=0.5 if Wflowed
# Normalize both poly and Wpoly to be out of Nc rather than 1
# poly format: MDTU,poly_mod              (normalized to 1)
# Wpoly format: MDTU,c=0.2,0.3,0.4,0.5    (normalized to Nc)
dat = []
toOpen = high_dir + '/data/' + poly_file + '.csv'   # High start
for line in open(toOpen):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= high_cut:
    continue
  if poly_file.startswith('W'):
    dat.append(float(temp[4]))
  else:
    dat.append(Nc * float(temp[1]))

toOpen = low_dir + '/data/' + poly_file + '.csv'    # Low start
for line in open(toOpen):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= low_cut:
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

# If xmax<=0 then autozoom, otherwise use given xmax
if xmax > 0.0:
  plt.xlim([0, xmax])
plt.grid(False)

if poly_file.startswith('W'):   # Record t corresponding to c=0.5
  if 'nt4' in cwd:
    tmax = '0.5'
  elif 'nt6' in cwd:
    tmax = '1.13'
  elif 'nt8' in cwd:
    tmax = '2'
  elif 'nt12' in cwd:
    tmax = '4.5'
  else:
    print "ERROR: Didn't recognize Nt of current working directory"
    sys.exit(1)

  plt.title('Wflowed Ploop mag, ' + title + ', t=' + tmax)
  plt.xlabel('|PL_W|')
else:
  plt.title('Polyakov loop magnitude, ' + title)
  plt.xlabel('|PL|')
plt.ylabel('Count')
plt.legend()

# Save a pdf
outfile = 'poly_hist_b' + beta + '_m' + mass + '.pdf'
if poly_file.startswith('W'):
  outfile = 'W' + outfile
plt.savefig(outfile, bbox_inches='tight')   # Reduce surrounding whitespace
# ------------------------------------------------------------------
