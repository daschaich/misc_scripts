#!/usr/bin/python
import glob
import math
import os
import sys
# ------------------------------------------------------------------
# Print out summary of all partially-quenched measurements
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isfile('ensembles'):
  print "ERROR: ensembles list does not exist"
  sys.exit(1)
elif not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

tags = []
for ensemble in open('ensembles'):
  tags.append(ensemble.rstrip())

# Print header
print "|--------------------------------------------------------|"
print "| Ensemble                    | m_v     | range    | N   |"
print "|--------------------------------------------------------|"

# Now cycle through ensembles, seeing which have measurements
tot = 0
for tag in tags:
  files = 'Out/pbp_' + tag + '-*'
  allFiles = glob.glob(files)
  totFiles = len(allFiles)
  tot += totFiles
  if totFiles > 0:
    # Figure out what valence masses have been measured
    allmV = []
    for filename in allFiles:
      # mV is in between last '-' and last '.'
      temp = (filename.split('-'))[-1]    # Everything after last '-'
      mV = '0.' + (temp.split('.'))[1]    # Have two '.' around mV
      if not mV in allmV:
        allmV.append(mV)

    # Easy case -- only a single valence mass for this ensemble
    if len(allmV) == 1:
      # Figure out range
      allMeas = []
      for filename in allFiles:
        cfg = (filename.split('.'))[-1]   # Everything after the last '.'
        allMeas.append(int(cfg))

      print "| %s" % tag.ljust(27),
      print "| %s" % mV.ljust(7),
      print "| %3d--%s" % (min(allMeas), str(max(allMeas)).ljust(3)),
      if mV == (tag.split('_'))[-1]:
        print "| %3d | Duplicates unitary measurements" % totFiles
      else:
        print "| %3d |" % totFiles
      print "|--------------------------------------------------------|"

    # Awkward case -- multiple valence masses for this ensemble
    elif len(allmV) > 1:
      allmV.sort(key=float)
      iter = 1
      print "| %s" % tag.ljust(27),
      for mV in allmV:          # Go through each mV
        if iter > 1 and iter < len(allmV):
          print "|                            ",
        elif iter == len(allmV):
          print "| Subtotal: %3d              " % totFiles,
        print "| %s" % mV.ljust(7),

        # Figure out range for this mV
        Nfilenames = 'Out/pbp_' + tag + '-' + mV + '.*'
        Nfiles = glob.glob(Nfilenames)
        subtot = len(Nfiles)
        allMeas = []
        for filename in Nfiles:
          cfg = (filename.split('.'))[-1]   # Everything after the last '.'
          allMeas.append(int(cfg))
        print "| %3d--%s" % (min(allMeas), str(max(allMeas)).ljust(3)),
        if mV == (tag.split('_'))[-1]:
          print "| %3d | Duplicates unitary measurements" % subtot
        else:
          print "| %3d |" % subtot
        iter += 1
      print "|--------------------------------------------------------|"
print "| Total: %4d                                            |" % tot
print "|--------------------------------------------------------|"
# ------------------------------------------------------------------
