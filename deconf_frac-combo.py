#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse dygraph data files to compute SU(4) deconfinement fraction
# No blocking, no uncertainties for fixed theta

# Combine high- and low-start runs
# Save output in both runs' results directories

# Parse arguments: first specify ensemble by beta and mass
# Then give both thermalization cuts, and finally the angle cut theta
# Require 0 <= theta <= pi/4
if len(sys.argv) < 6:
  print "Usage:", str(sys.argv[0]), "<beta> <mass>",
  print "<high cut> <low cut> <theta>"
  sys.exit(1)
beta = str(sys.argv[1])
mass = str(sys.argv[2])
hi_cut = int(sys.argv[3])
lo_cut = int(sys.argv[4])
theta = float(sys.argv[5])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Set up directories from path and input
hi_dir = 'b' + beta + '_high_m' + mass
lo_dir = 'b' + beta + '_low_m' + mass

# First make sure we're calling this from the right place
toCheck = [hi_dir + '/data', lo_dir + '/data']
for i in toCheck:
  if not os.path.isdir(i):
    print "ERROR:", i, "does not exist, aborting"
    sys.exit(1)
  obsfile = i + '/poly_arg.csv'
  if not os.path.isfile(obsfile):
    print "ERROR:", i, "does not exist, aborting"
    sys.exit(1)

# For convenience
pi_ov_four = 0.25 * np.pi
pi_ov_two = 0.5 * np.pi
three_pi_ov_four = 0.75 * np.pi
five_pi_ov_four = 1.25 * np.pi
three_pi_ov_two = 1.5 * np.pi
seven_pi_ov_four = 1.75 * np.pi
two_pi = 2.0 * np.pi
shift = theta / pi_ov_four
norm = pi_ov_four / (pi_ov_four - theta)

# Quick sanity check on input angle cut theta
if theta < 0.0 or theta > pi_ov_four:
  print "ERROR: theta=%.4g out of range [0, pi/4], aborting" % theta
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For poly_arg, have only one datum per line
# Although counts are integer, make them floats for division below
N_in = 0.0
N_out = 0.0

# First high-start run
obsfile = hi_dir + '/data/poly_arg.csv'
for line in open(obsfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= hi_cut:
    continue

  # Exactly one Z4 vacuum should be within pi_ov_four radians
  # Convert range from [-pi, pi) to [0, 2pi)
  arg = float(temp[1]) + np.pi
  if arg >= two_pi:                 # Sanity check
    print "ERROR: arg out of range at MDTU=%d, aborting" % int(temp[0])
    sys.exit(1)

  found = 0
  if arg <= pi_ov_four:             # Closest to theta=0 (part I)
    phi = arg
    found += 1
  elif arg <= three_pi_ov_four:     # Closest to theta=pi/2
    phi = np.fabs(arg - pi_ov_two)
    found += 1
  elif arg <= five_pi_ov_four:      # Closest to theta=pi
    phi = np.fabs(arg - np.pi)
    found += 1
  elif arg <= seven_pi_ov_four:     # Closest to theta=3pi/2
    phi = np.fabs(arg - three_pi_ov_two)
    found += 1
  else:                             # Closest to theta=0 (part II)
    phi = two_pi - arg
    found += 1

  # Two sanity checks that should never be triggered
  if not found == 1:
    print "ERROR: Couldn't categorize arg=%.4g at MDTU=%d, aborting" \
          % (arg, int(temp[0]))
    sys.exit(1)
  if phi < 0.0 or phi > pi_ov_four:
    print "ERROR: phi=%.4g out of range [0, pi/4] at MDTU=%d, aborting" \
          % (phi, int(temp[0]))
    sys.exit(1)

  # Record whether phi is "in" or "out" of small range [0, theta]
  if phi <= theta:
    N_in += 1.0
  else:
    N_out += 1.0

# Add the low-start run
obsfile = lo_dir + '/data/poly_arg.csv'
for line in open(obsfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= lo_cut:
    continue

  # Exactly one Z4 vacuum should be within pi_ov_four radians
  # Convert range from [-pi, pi) to [0, 2pi)
  arg = float(temp[1]) + np.pi
  if arg >= two_pi:                 # Sanity check
    print "ERROR: arg out of range at MDTU=%d, aborting" % int(temp[0])
    sys.exit(1)

  found = 0
  if arg <= pi_ov_four:             # Closest to theta=0 (part I)
    phi = arg
    found += 1
  elif arg <= three_pi_ov_four:     # Closest to theta=pi/2
    phi = np.fabs(arg - pi_ov_two)
    found += 1
  elif arg <= five_pi_ov_four:      # Closest to theta=pi
    phi = np.fabs(arg - np.pi)
    found += 1
  elif arg <= seven_pi_ov_four:     # Closest to theta=3pi/2
    phi = np.fabs(arg - three_pi_ov_two)
    found += 1
  else:                             # Closest to theta=0 (part II)
    phi = two_pi - arg
    found += 1

  # Two sanity checks that should never be triggered
  if not found == 1:
    print "ERROR: Couldn't categorize arg=%.4g at MDTU=%d, aborting" \
          % (arg, int(temp[0]))
    sys.exit(1)
  if phi < 0.0 or phi > pi_ov_four:
    print "ERROR: phi=%.4g out of range [0, pi/4] at MDTU=%d, aborting" \
          % (phi, int(temp[0]))
    sys.exit(1)

  # Record whether phi is "in" or "out" of small range [0, theta]
  if phi <= theta:
    N_in += 1.0
  else:
    N_out += 1.0

# Now print deconfinement fraction normalized to 1
tot = N_in + N_out
frac = norm * (N_in / tot - shift)
outfilename = hi_dir + '/results/deconf_frac-combo-' + str(theta) + '.dat'
outfile_hi = open(outfilename, 'w')
print >> outfile_hi, "%.8g # %d" % (frac, tot)
outfile_hi.close()
outfilename = lo_dir + '/results/deconf_frac-combo-' + str(theta) + '.dat'
outfile_lo = open(outfilename, 'w')
print >> outfile_lo, "%.8g # %d" % (frac, tot)
outfile_lo.close()
# ------------------------------------------------------------------



# ----------------------------------------------------------------
# For the Wilson-flowed Polyakov loop argument
# let's print out all four of c=0.2, 0.3, 0.4 and 0.5
N_in = [0.0 for x in range(4)]
N_out = [0.0 for x in range(4)]

# First high-start run
flowfile = hi_dir + '/data/Wpoly_arg.csv'
for line in open(flowfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= hi_cut:
    continue

  # Exactly one Z4 vacuum should be within pi_ov_four radians
  # Convert range from [-pi, pi) to [0, 2pi)
  for i in range(4):
    arg = float(temp[i + 1]) + np.pi
    if arg >= two_pi:                 # Sanity check
      print "ERROR: arg out of range at MDTU=%d, aborting" % int(temp[0])
      sys.exit(1)

    found = 0
    if arg <= pi_ov_four:             # Closest to theta=0 (part I)
      phi = arg
      found += 1
    elif arg <= three_pi_ov_four:     # Closest to theta=pi/2
      phi = np.fabs(arg - pi_ov_two)
      found += 1
    elif arg <= five_pi_ov_four:      # Closest to theta=pi
      phi = np.fabs(arg - np.pi)
      found += 1
    elif arg <= seven_pi_ov_four:     # Closest to theta=3pi/2
      phi = np.fabs(arg - three_pi_ov_two)
      found += 1
    else:                             # Closest to theta=0 (part II)
      phi = two_pi - arg
      found += 1

    # Two sanity checks that should never be triggered
    if not found == 1:
      print "ERROR: Couldn't categorize arg=%.4g at MDTU=%d, aborting" \
            % (arg, int(temp[0]))
      sys.exit(1)
    if phi < 0.0 or phi > pi_ov_four:
      print "ERROR: phi=%.4g out of range [0, pi/4] at MDTU=%d, aborting" \
            % (phi, int(temp[0]))
      sys.exit(1)

    # Record whether phi is "in" or "out" of small range [0, theta]
    if phi <= theta:
      N_in[i] += 1.0
    else:
      N_out[i] += 1.0

# Add the low-start run
flowfile = lo_dir + '/data/Wpoly_arg.csv'
for line in open(flowfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= lo_cut:
    continue

  # Exactly one Z4 vacuum should be within pi_ov_four radians
  # Convert range from [-pi, pi) to [0, 2pi)
  for i in range(4):
    arg = float(temp[i + 1]) + np.pi
    if arg >= two_pi:                 # Sanity check
      print "ERROR: arg out of range at MDTU=%d, aborting" % int(temp[0])
      sys.exit(1)

    found = 0
    if arg <= pi_ov_four:             # Closest to theta=0 (part I)
      phi = arg
      found += 1
    elif arg <= three_pi_ov_four:     # Closest to theta=pi/2
      phi = np.fabs(arg - pi_ov_two)
      found += 1
    elif arg <= five_pi_ov_four:      # Closest to theta=pi
      phi = np.fabs(arg - np.pi)
      found += 1
    elif arg <= seven_pi_ov_four:     # Closest to theta=3pi/2
      phi = np.fabs(arg - three_pi_ov_two)
      found += 1
    else:                             # Closest to theta=0 (part II)
      phi = two_pi - arg
      found += 1

    # Two sanity checks that should never be triggered
    if not found == 1:
      print "ERROR: Couldn't categorize arg=%.4g at MDTU=%d, aborting" \
            % (arg, int(temp[0]))
      sys.exit(1)
    if phi < 0.0 or phi > pi_ov_four:
      print "ERROR: phi=%.4g out of range [0, pi/4] at MDTU=%d, aborting" \
            % (phi, int(temp[0]))
      sys.exit(1)

    # Record whether phi is "in" or "out" of small range [0, theta]
    if phi <= theta:
      N_in[i] += 1.0
    else:
      N_out[i] += 1.0

# Now print deconfinement fractions normalized to 1
outfilename = hi_dir + '/results/Wdeconf_frac-combo-' + str(theta) + '.dat'
outfile_hi = open(outfilename, 'w')
print >> outfile_hi, "# c=0.2 c=0.3 c=0.4 c=0.5 # Nmeas"
outfilename = lo_dir + '/results/Wdeconf_frac-combo-' + str(theta) + '.dat'
outfile_lo = open(outfilename, 'w')
print >> outfile_lo, "# c=0.2 c=0.3 c=0.4 c=0.5 # Nmeas"
for i in range(4):
  tot = N_in[i] + N_out[i]
  frac = norm * (N_in[i] / tot - shift)
  print >> outfile_hi, "%.8g" % frac,
  print >> outfile_lo, "%.8g" % frac,
print >> outfile_hi, "# %d" % tot
print >> outfile_lo, "# %d" % tot
outfile_hi.close()
outfile_lo.close()
# ------------------------------------------------------------------
