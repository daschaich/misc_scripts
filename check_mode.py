#!/usr/bin/python
import glob
import math
import os
import sys
# ------------------------------------------------------------------
# Print out summary of all stochastic mode number measurements
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
print "|---------------------------------------------------------------------------|"
print "| Ensemble                    | label                | range          | N   |"
print "|---------------------------------------------------------------------------|"

# Now cycle through ensembles, seeing which have measurements
tot = 0
for tag in tags:
  files = 'Out/mode_*' + tag + '*'
  allFiles = glob.glob(files)
  totFiles = len(allFiles)
  tot += totFiles
  if totFiles > 0:
    # Figure out what labels have been used
    labels = []
    for filename in allFiles:
      temp = (filename.split(tag))[1]   # Everything after the tag
      label = (temp.split('.'))[:-1]
      # There really should be a better way to do this
      if len(label) > 1:    # This happens if the label contains '.'
        temp = ''
        for i in range(len(label)):
          temp += label[i]
          if i + 1 < len(label):
            temp += '.'
        label = temp
      else:
        temp = label[0]
        label = temp
      if not label in labels:
        labels.append(label)

    # Easy case -- only a single label for this ensemble
    if len(labels) == 1:
      # Figure out range
      allMeas = []
      for filename in allFiles:
        cfg = (filename.split('.'))[-1]   # Everything after the last '.'
        allMeas.append(int(cfg))

      print "| %s" % tag.ljust(27),
      print "| %s" % labels[0].ljust(20),
      print "| %3d--%s" % (min(allMeas), str(max(allMeas)).ljust(3)),
      possible = max(allMeas) - min(allMeas) + 1
      print "(%s)" % str(possible).rjust(3),
      print "| %3d |" % totFiles
      print "|---------------------------------------------------------------------------|"

    # Awkward case -- multiple labels for this ensemble
    elif len(labels) > 1:
      labels.sort(key=str)
      iter = 1
      print "| %s" % tag.ljust(27),
      for label in labels:          # Go through each label
        if iter > 1 and iter < len(label):
          print "|                            ",
        elif iter == len(labels):
          print "| Subtotal: %3d              " % totFiles,
        print "| %s" % label.ljust(20),

        # Figure out range for this label
        if label is '':
          Nfilenames = 'Out/mode_*' + tag + '.*'
        else:
          Nfilenames = 'Out/mode_*' + tag + '*' + label + '.*'
        Nfiles = glob.glob(Nfilenames)
        subtot = len(Nfiles)
        allMeas = []
        for filename in Nfiles:
          cfg = (filename.split('.'))[-1]   # Everything after the last '.'
          allMeas.append(int(cfg))
        print "| %3d--%s" % (min(allMeas), str(max(allMeas)).ljust(3)),
        possible = max(allMeas) - min(allMeas) + 1
        print "(%s)" % str(possible).rjust(3),
        print "| %3d |" % subtot
        iter += 1
      print "|---------------------------------------------------------------------------|"
print "| Total: %4d                                                               |" % tot
print "|---------------------------------------------------------------------------|"
# ------------------------------------------------------------------
