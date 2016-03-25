#!/usr/bin/python
import glob
import math
import os
import sys
# ------------------------------------------------------------------
# Print out summary of all eigenvalue measurements
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isfile('ensembles'):
  print "ERROR: ensembles list does not exist"
  sys.exit(1)
elif not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

tags = []
for ensemble in open('ensembles'):
  tags.append(ensemble.rstrip())

# Print header
print "|-----------------------------------------------------------|"
print "| Ensemble                    | Nvec | range          | N   |"
print "|-----------------------------------------------------------|"

# Now cycle through ensembles, seeing which have measurements
tot = 0
for tag in tags:
  files = 'Out/eig_*_' + tag + '.*'
  allFiles = glob.glob(files)
  totFiles = len(allFiles)
  tot += totFiles
  if totFiles > 0:
    # Figure out what Nvecs have been measured
    allNvecs = []
    for filename in allFiles:
      Nvecs = (filename.split('_'))[1]
      if not Nvecs in allNvecs:
        allNvecs.append(Nvecs)

    # Easy case -- only a single Nvecs for this ensemble
    if len(allNvecs) == 1:
      # Figure out range
      allMeas = []
      for filename in allFiles:
        cfg = (filename.split('.'))[-1]   # Everything after the last '.'
        allMeas.append(int(cfg))

      print "| %s" % tag.ljust(27),
      print "| %4d" % int(allNvecs[0]),
      print "| %3d--%s" % (min(allMeas), str(max(allMeas)).ljust(3)),
      possible = max(allMeas) - min(allMeas) + 1
      print "(%s)" % str(possible).rjust(3),
      if totFiles == possible:
        print "| %3d |" % totFiles
      else:
        diff = possible - totFiles
        print "| %3d | (Missing %d)" % (totFiles, diff)
      print "|-----------------------------------------------------------|"

    # Awkward case -- multiple Nvecs for this ensemble
    elif len(allNvecs) > 1:
      allNvecs.sort(key=int)
      iter = 1
      print "| %s" % tag.ljust(27),
      for Nvecs in allNvecs:          # Go through each Nvecs
        if iter > 1 and iter < len(allNvecs):
          print "|                            ",
        elif iter == len(allNvecs):
          print "| Subtotal: %3d              " % totFiles,
        print "| %4d" % int(Nvecs),

        # Figure out range for this Nvecs
        Nfilenames = 'Out/eig_' + Nvecs + '_' + tag + '.*'
        Nfiles = glob.glob(Nfilenames)
        subtot = len(Nfiles)
        allMeas = []
        for filename in Nfiles:
          cfg = (filename.split('.'))[-1]   # Everything after the last '.'
          allMeas.append(int(cfg))
        print "| %3d--%s" % (min(allMeas), str(max(allMeas)).ljust(3)),
        possible = max(allMeas) - min(allMeas) + 1
        print "(%s)" % str(possible).rjust(3),
        if subtot == possible:
          print "| %3d |" % subtot
        else:
          diff = possible - subtot
          print "| %3d | (Missing %d)" % (subtot, diff)
        iter += 1
      print "|-----------------------------------------------------------|"
print "| Total: %4d                                               |" % tot
print "|-----------------------------------------------------------|"
# ------------------------------------------------------------------
