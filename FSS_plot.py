#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# ------------------------------------------------------------------
# This script produces a curve collapse animation
# for the given directory and observable

# Animation based on BSD-licensed example by Jake Vanderplas,
# http://jakevdp.github.com, vanderplas@astro.washington.edu
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse argument: file to analyze
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<file> <y_m> <c0>"
  sys.exit(1)
filename = str(sys.argv[1])
ym = float(sys.argv[2])
c0 = float(sys.argv[3])
runtime = -time.time()

if not os.path.isfile(filename):
  print "ERROR:", filename, "does not exist"
  sys.exit(1)

# !!! Important definitions
y0 = -0.3         # Scaling dimension of gauge coupling
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read and parse data
# Format: mf MH err beta L Nt obs
beta = -1       # Check
obs = 'obs'     # Check
all_L = []      # System sizes (assuming constant aspect ratio)
mf = []         # Fermion mass -- list of lists
MH = []         # Hadron mass (or other observable) -- list of lists
err = []        # Error on observable -- list of lists
for line in open(filename):
  if line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  if beta == -1:
    beta = float(temp[3])
    obs = temp[6]
  # Check that data seem to be in order
  elif beta != float(temp[3]):
    print "ERROR: the coupling has changed from %.2g to %.2g" \
          % (beta, float(temp[3]))
    sys.exit(1)
  elif obs != temp[6]:
    print "ERROR: the observable has changed from %s to %s" \
          % (obs, temp[6])
    sys.exit(1)

  # Now we should be good to read the line
  # First consider the case that this L has already been seen
  # Allow for the data to be out of order
  if temp[4] in all_L:
    toAdd = all_L.index(temp[4])
    (mf[toAdd]).append(float(temp[0]))
    (MH[toAdd]).append(float(temp[1]))
    (err[toAdd]).append(float(temp[2]))
  else:                         # Add new lists to mf, MH and err
    all_L.append(temp[4])       # Save as string for comparison above
    mf.append([])
    (mf[-1]).append(float(temp[0]))
    MH.append([])
    (MH[-1]).append(float(temp[1]))
    err.append([])
    (err[-1]).append(float(temp[2]))
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# For now, just display optimal and exit
myc = ['red', 'blue', 'green', 'black', 'purple']   # Colors
mym = ['x', 'o', 's', '.', 'D']                     # Markers
fig, collapse = plt.subplots(1)
title = filename + " curve collapse at $y_0=$" + str(y0)
readin = ", $c_0=$" + str(c0) + ", $y_m=$" + str(ym)
collapse.set_title(r'' + title + readin)
collapse.set_xlabel(r'$L m^{1 / y_m}$')
collapse.set_ylabel(r'$M L$')

# Find plot range by cycling over each point j in each data set i
omega = -y0 / ym
maxX = 0
minY = 999
maxY = 0
for i in range(len(all_L)):
  L = float(all_L[i])
  tag = "$L=" + all_L[i] + "$"
  temp = [L * np.power(m, 1 / ym) for m in mf[i]]
  x = np.array(temp)
  temp = [L * MH[i][j] / (1 + c0 * np.power(mf[i][j], omega)) for j in range(len(mf[i]))]
  y = np.array(temp)
  temp = [L * err[i][j] / (1 + c0 * np.power(mf[i][j], omega)) for j in range(len(mf[i]))]
  yerr = np.array(temp)

  if np.amax(x) > maxX: maxX = np.amax(x)
  if np.amin(y) < minY: minY = np.amin(y)
  if np.amax(y) > maxY: maxY = np.amax(y)
  collapse.errorbar(x, y, yerr=yerr, c=myc[i], marker=mym[i], linestyle='None', label=tag)

collapse.legend(loc='best')
collapse.set_xlim(0, maxX * 1.1)
collapse.set_ylim(minY / 1.1, maxY * 1.1)
runtime += time.time()
print "Runtime: %.2g seconds" % runtime
plt.show()
sys.exit(0)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Function to animate plot depends on ym
def curve(i):
  ym = 1 + float(i) / granularity

  # !!! Need multiple data sets?
  toplot.set_data(t, tSqE)
  toplot.set_label(label)
  return toplot,
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Animation stuff
ymax = 3
ymin = 1
granularity = 0.01
Npts = int((ymax - ymin) / granularity) + 1

# Set up the figure, the axis, and the plot element to be animated
fig = plt.figure()
ax = plt.axes(xlim=(0, 5), ylim=(0, 5))   # !!! Make this fancier
plt.xlabel('Lm^{1/y}')
plt.ylabel('M*L')
toplot, = ax.plot([], [], ls='None', marker='.', label="label")
#plt.legend()

# Initialization function: plot the background of each frame
def init():
  toplot.set_data([], [])
  return toplot,

# Go.  blit=True means only re-draw the parts that have changed
# Interval is the time delay between frames, in ms
anim = animation.FuncAnimation(fig, curve, init_func=init, repeat=False,
                               frames=201, interval=1, blit=True)

# Save as mp4
#anim.save('therm.mp4', writer='mencoder', fps=20)

plt.show()
# ------------------------------------------------------------------
