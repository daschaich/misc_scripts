#!/usr/bin/python
import glob
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
import sys
# ------------------------------------------------------------------
# This script plots the Wilson flow t^2<E> as a function of flow time,
# averaging over a sliding window in MD time to illustrate thermalization

# Animation based on BSD-licensed example by Jake Vanderplas,
# http://jakevdp.github.com, vanderplas@astro.washington.edu
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments: first is directory, second is sliding window width,
# third is file name, up to the # at the end
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<dir> <bin [MDTU]> <tag [out.Wflow]>"
  sys.exit(1)
dirname = str(sys.argv[1])
windowSize = int(sys.argv[2])
filetag = str(sys.argv[3])

if not os.path.isdir(dirname):
  print "ERROR:", dirname, "does not exist"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make and sort a list of all the dirname/filetag.* files
cfgs = []
temp = dirname + '/' + filetag + '.*'
for filename in glob.glob(temp):
  cfgs.append(int((filename.split('.'))[-1]))  # Number after last .
cfgs.sort()

# Sanity check: monitor t to make sure we combine the right data
# Also grab (all!) tSqE to estimate y-axis range
t = []
tSqE = []
for i in cfgs:
  filename = filetag + '.' + str(i)
  toOpen = os.path.join(dirname, filename)
  for line in open(toOpen):
    if line.startswith('WFLOW '):
      temp = line.split()
      t.append(float(temp[1]))
      tSqE.append(float(temp[4]))
xmax = t[-1] + 1
ymax = 1.1 * max(tSqE)
#print ymax

# Reset to get right t
t = []
tSqE = []
filename = filetag + '.' + str(cfgs[-1])
toOpen = os.path.join(dirname, filename)
for line in open(toOpen):
  if line.startswith('WFLOW '):
    temp = line.split()
    t.append(float(temp[1]))
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# See how many windows we have, each gets a frame
iter = 0
start = cfgs[iter]
while start + windowSize < cfgs[-1]:
  iter += 1
  start = cfgs[iter]
Nframes = iter
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Function to animate plot constructs sliding window
def slide(i):
  # Set up window
  window = [cfgs.pop(0)]    # Remove smallest cfg number from list
  iter = 0
  end = cfgs[iter]
  while end < window[0] + windowSize:
    window.append(end)
    iter += 1
    end = cfgs[iter]
  label = str(window[0]) + "--" + str(window[0])

  # Parse the files in the window
  tSqE = np.zeros(len(t), dtype = np.float)
  tSqE_err = np.zeros(len(t), dtype = np.float)
  for cfg in window:
    filename = filetag + '.'
    filename += str(cfg)
    toOpen = os.path.join(dirname, filename)
    iter = 0
    for line in open(toOpen):
      if line.startswith('WFLOW '):
        temp = line.split()   # Convert line into list

        # First a sanity check
        toCheck = float(temp[1])
        if toCheck != t[iter]:
          print "ERROR: FLOW TIME MISMATCH %.2g %.2g" % (toCheck, t[iter])

        # Now we are good
        toAdd = float(temp[4])
        tSqE[iter] += toAdd
        tSqE_err[iter] += toAdd * toAdd
        iter += 1

  # Now print averages and standard deviations
  for iter in range(0, len(t)):
    tSqE[iter] /= len(window)
    tSqE_err[iter] /= len(window)
    tSqE_err[iter] -= tSqE[iter] * tSqE[iter]
    tSqE_err[iter] = np.sqrt(tSqE_err[iter] / (len(window) - 1))
#    print('{0:.4g} {1:.5g} {2:.5g}'.format(t[iter], tSqE[iter], \
#                                           tSqE_err[iter]))
  toplot.set_data(t, tSqE)
  toplot.set_label(label)
  return toplot,
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Animation stuff
# First set up the figure, the axis, and the plot element to be animated
fig = plt.figure()
#ax = plt.axes(xlim=(30, 32), ylim=(0.71, 0.74))
ax = plt.axes(xlim=(0, xmax), ylim=(0, ymax))
plt.xlabel('t')
plt.ylabel('t^2<E>')
toplot, = ax.plot([], [], ls='None', marker='.', label="label")
#plt.legend()

# Initialization function: plot the background of each frame
def init():
  toplot.set_data([], [])
  return toplot,# errors

# Go.  blit=True means only re-draw the parts that have changed
# Interval is the time delay between frames, in ms
anim = animation.FuncAnimation(fig, slide, init_func=init, repeat=False,
                               frames=Nframes, interval=10, blit=True)

# Save as mp4
#anim.save('therm.mp4', writer='mencoder', fps=20)

plt.show()
# ------------------------------------------------------------------
