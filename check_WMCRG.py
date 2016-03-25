#!/usr/bin/python
import glob
import math
import os
import sys
import re
# ------------------------------------------------------------------
# Print out summary of all Wilson-flowed MCRG measurements
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
print "|-------------------------------------------------------------|"
print "| Ensemble                    | scheme | range          | N   |"
print "|-------------------------------------------------------------|"

# Now cycle through ensembles, seeing which have measurements
tot = 0
for tag in tags:
  files = 'Out/WMCRG*_' + tag + '.*'
  allFiles = glob.glob(files)
  totFiles = len(allFiles)
  tot += totFiles
  if totFiles > 0:
    # Figure out what schemes have been measured
    allschemes = []
    for filename in allFiles:
      # Luckily, this seems to pick up the correct (first) '_'
      test = re.search('WMCRG(.+?)_', filename)
      if test:
        scheme = test.group(1)
      else:
        print "ERROR: couldn't extract scheme from", filename
        sys.exit(1)
#      print filename, "-->", scheme
      if not scheme in allschemes:
        allschemes.append(scheme)

    # Unusual case -- only a single scheme for this ensemble
    if len(allschemes) == 1:
      # Figure out range
      allMeas = []
      for filename in allFiles:
        cfg = (filename.split('.'))[-1]   # Everything after the last '.'
        allMeas.append(int(cfg))

      print "| %s" % tag.ljust(27),
      print "|     %2d" % int(allschemes[0]),
      print "| %3d--%s" % (min(allMeas), str(max(allMeas)).ljust(3)),
      possible = max(allMeas) - min(allMeas) + 1
      print "(%s)" % str(possible).rjust(3),
      if totFiles == possible:
        print "| %3d |" % totFiles
      else:
        diff = possible - totFiles
        print "| %3d | (Missing %d)" % (totFiles, diff)
      print "|-------------------------------------------------------------|"

    # Expected case -- multiple schemes for this ensemble
    elif len(allschemes) > 1:
      allschemes.sort(key=int)
      iter = 1
      print "| %s" % tag.ljust(27),
      for scheme in allschemes:          # Go through each scheme
        if iter > 1 and iter < len(allschemes):
          print "|                            ",
        elif iter == len(allschemes):
          print "| Subtotal: %3d              " % totFiles,
        print "|     %2d" % int(scheme),

        # Figure out range for this scheme
        Nfilenames = 'Out/WMCRG' + scheme + '_' + tag + '.*'
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
      print "|-------------------------------------------------------------|"
print "| Total: %4d                                                 |" % tot
print "|-------------------------------------------------------------|"
# ------------------------------------------------------------------
