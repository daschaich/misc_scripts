#!/usr/bin/python
import glob
import os
import sys
# ------------------------------------------------------------------
# Redo Wilson flow interpolations in the given output file
# Print to terminal for redirection to new file

# Parse argument: the file to consider
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
toOpen = str(sys.argv[1])

# Make sure we're calling this from the right place
if not os.path.isfile(toOpen):
  print "ERROR:", toOpen, "does not exist"
  sys.exit(1)

# Fix dt=0.01, which also implies t>0
dt = 0.01

# Cycle over lines in file
check = -1
toInterp = -1
for line in open(filename):
  if line.startswith('RUNNING COMPLETED'):
    if check == 1:    # Check that we have one measurement per file
      print filename, "reports two measurements"
    check = 1

  # Print all the standard stuff
  if not line.startswith('WFLOW '):
    print line.rstrip()
  # Deal with the real stuff
  # We will linearly interpolate t^2<E>, not <E> itself
  else:   # Format: WFLOW t plaq E tSqE der_tSqE check topo
    if 'interp' in line:
      toInterp = 1        # We will need to interpolate in the future
    else:
      if toInterp < 0:    # We can just print and save values
        print line.rstrip()
        temp = line.split()
        t_old = float(temp[1])
        tSqE_old = float(temp[4])
        check_old = float(temp[6])
        topo_old = float(temp[7])
      else:               # We need to interpolate up to current line
        temp = line.split()
        t = float(temp[1])
        interval = t - t_old

        slope_tSqE = (float(temp[4]) - tSqE_old) / interval
        slope_check = (float(temp[6]) - check_old) / interval
        slope_topo = (float(temp[7]) - topo_old) / interval
        prev = tSqE_old
        for i in np.arange(t_old + dt, t, dt, endpoint=False):
          delta = i - t_old
          tSqE = tSqE_old + delta * slope_tSqE
          E = tSqE / (i * i)
          der = i * (tSqE - prev) / dt
          prev = tSqE

          check = check_old + delta * slope_tSqE
          plaq = check / (12.0 * i * i)

          topo = topo_old + delta * slope_topo
          print "WFLOW %g %g %g %g %g %g %g\n" \
                % (i, plaq, E, tSqE, der, check, topo);

        # Now we are caught up and can print the current line
        print line.rstrip()

        # Prepare for future interpolations
        toInterp = -1
        t_old = float(temp[1])
        tSqE_old = float(temp[4])
        check_old = float(temp[6])
        topo_old = float(temp[7])

if check == -1:
  print filename, "did not complete"
  continue
# ------------------------------------------------------------------
