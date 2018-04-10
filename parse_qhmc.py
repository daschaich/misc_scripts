#!/usr/bin/python
import os
import sys
import glob
import math
# ------------------------------------------------------------------
# Parse QHMC output files for a single ensemble,
# shuffling the extracted data into dedicated files for plotting

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Try to deal with the forces
# First Nf/4 (counting from 01) are for outer lattice fields
# Next Nf/4 are for Hasenbusch lattice fields
temp = os.getcwd()
if '4f' in temp:
  Nc=4.0
  Ftag = 'FF[02]'
elif '8f' in temp:
  Nc=3.0
  Ftag = 'FF[03]'
else:
  print "ERROR: Force extraction can only handle Nf=4 or Nf=8"
  sys.exit(1)

ERRFILE = open('ERRORS', 'w')
MISSINGFILES = open('MISSING', 'w')

# Physical observables
PLAQ = open('data/plaq.csv', 'w')
print >> PLAQ, "MDTU,plaq_ss,plaq_st"
PBP = open('data/pbp.csv', 'w')
print >> PBP, "MDTU,Re(pbp),Im(pbp)"
POLY = open('data/poly.csv', 'w')
print >> POLY, "ReTr(L),ImTr(L)"
XPOLY = open('data/xpoly.csv', 'w')
print >> XPOLY, "ReTr(W_0),ImTr(W_0)"
POLY_R = open('data/poly_r.csv', 'w')
print >> POLY_R, "MDTU,ReTr(L)"
POLY_MOD = open('data/poly_mod.csv', 'w')
print >> POLY_MOD, "MDTU,|Tr(L)|"
XPOLY_R = open('data/xpoly_r.csv', 'w')
print >> XPOLY_R, "MDTU,ReTr(W_0)"
XPOLY_MOD = open('data/xpoly_mod.csv', 'w')
print >> XPOLY_MOD, "MDTU,|Tr(W_0)|"
POLY_ARG = open('data/poly_arg.csv', 'w')
print >> POLY_ARG, "MDTU,arg(Tr(L))"
XPOLY_ARG = open('data/xpoly_arg.csv', 'w')
print >> XPOLY_ARG, "MDTU,arg(Tr(W_0))"
PLAQ_DIFF = open('data/plaq_diff.csv', 'w')
print >> PLAQ_DIFF, "MDTU,t,x,y,z,Norm"
LINK_DIFF = open('data/link_diff.csv', 'w')
print >> LINK_DIFF, "MDTU,t,x,y,z,Norm"
EIG = open('data/eig.csv', 'w')
print >> EIG, "MDTU,1,2,3,4,5,6,7,8,9,10,11,12"
WFLOW = open('data/Wflow.csv', 'w')
print >> WFLOW, "MDTU,c=0.2,c=0.25,c=0.3,c=0.35"
gprop = 128.0 * 3.14159**2 / (3.0 * 8.0)  # Only compute this once
TOPO = open('data/topo.csv', 'w')
print >> TOPO, "MDTU,c=0.2,c=0.3,c=0.4,c=0.5"
WPOLY = open('data/Wpoly.csv', 'w')
print >> WPOLY, "MDTU,c=0.2,c=0.3,c=0.4,c=0.5"

# Blocked observables
PLAQB = open('data/plaqB.csv', 'w')
print >> PLAQB, "MDTU,bl0,bl1,bl2,bl3,bl4"
POLYB = open('data/polyB.csv', 'w')
print >> POLYB, "ReTr(L_b),bl0,bl1,bl2,bl3,bl4"
XPOLYB = open('data/xpolyB.csv', 'w')
print >> XPOLYB, "Re(W_0),bl0,bl1,bl2,bl3,bl4"
POLY_RB = open('data/poly_rB.csv', 'w')
print >> POLY_RB, "MDTU,bl0,bl1,bl2,bl3,bl4"
POLY_MODB = open('data/poly_modB.csv', 'w')
print >> POLY_MODB, "MDTU,bl0,bl1,bl2,bl3,bl4"
XPOLY_RB = open('data/xpoly_rB.csv', 'w')
print >> XPOLY_RB, "MDTU,bl0,bl1,bl2,bl3,bl4"
XPOLY_MODB = open('data/xpoly_modB.csv', 'w')
print >> XPOLY_MODB, "MDTU,bl0,bl1,bl2,bl3,bl4"
POLY_ARGB = open('data/poly_argB.csv', 'w')
print >> POLY_ARGB, "MDTU,bl0,bl1,bl2,bl3,bl4"
XPOLY_ARGB = open('data/xpoly_argB.csv', 'w')
print >> XPOLY_ARGB, "MDTU,bl0,bl1,bl2,bl3,bl4"

# Evolution observables
ACCP = open('data/accP.csv', 'w')
print >> ACCP, "t,accP"
EXP_DS = open('data/exp_dS.csv', 'w')
print >> EXP_DS, "t,e^(-dS)"
DELTAS = open('data/deltaS.csv', 'w')
print >> DELTAS, "t,deltaS"
ABS_DS = open('data/abs_dS.csv', 'w')
print >> ABS_DS, "t,|deltaS|"
FORCE = open('data/force.csv', 'w')
print >> FORCE, "t,F0,F1,Fgauge"
CG_ITERS = open('data/cg_iters.csv', 'w')
print >> CG_ITERS, "t,cg_iters"
WALLTIME = open('data/walltime.csv', 'w')
print >> WALLTIME, "t,walltime"
WALLTU = open('data/wallTU.csv', 'w')
print >> WALLTU, "t,cost"

# Run parameters
NSTEP = open('data/Nstep.csv', 'w')
print >> NSTEP, "t,N0,N1,Ngauge"
STEPSIZE = open('data/stepsize.csv', 'w')
print >> STEPSIZE, "t,eps0,eps1,eps_gauge"
MH = open('data/MH.csv', 'w')
print >> MH, "t,MH"
TLENGTH = open('data/tlength.csv', 'w')
print >> TLENGTH, "t,L"
KEY = open('data/key.csv', 'w')
print >> KEY, "t,file"
TU = open('data/TU.csv', 'w')
print >> TU, "t,MDTU"

# Status checks and running sums for the ensemble as a whole
oldcfg = 0
oldstamp = "start"
traj = 0;
MDTU = 0;
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Cycle through files, first "out" and then "Sout", "Wflow" and "eig"
# list.txt is a list of $load-$save
for temp_tag in open('list.txt'):
  tag = temp_tag.rstrip()
  load, cfg = tag.split('-')
  # Initialize running sums and set dummy walltime, and nsteps
  # If the walltime isn't overwritten, then the run died
  # or its output file is corrupted
  pbp_r = 0;
  pbp_iter = 0;
  walltime = -1;
  stamp = "start";

  # Open file
  # If not found, move on to next file instead of killing whole program,
  # but print error message so I know there is a problem
  infile = 'Out/out.' + tag
  if not os.path.isfile(infile):
    print "Problem opening", infile
    print >> ERRFILE, "Problem opening", infile
    continue    # Skip to next file
#  print infile  # Monitor running status

  # If not starting from first file in this ensemble,
  # or if we seem to have skipped a file,
  # guess approximate starting trajectory
  traj_per_file = -1;
  L = -1;
  Nt = -1;
  vol = -1;

  for line in open(infile):
    if line.startswith('latsize'):
      Nt = int((line.split())[-1])
      L = int((line.split())[-2])
      vol = L**3 * Nt
      # Set L to minimum of nx and nt for Wilson flow
      if Nt < L:
        L = Nt
      break       # Don't go through whole file yet
    elif line.startswith('ntraj'):
      traj_per_file = int((line.split())[1])
      endtraj = traj + traj_per_file

  if traj_per_file < 0:
    print infile, "never defines number of trajectories"
    print >> ERRFILE, infile, "never defines number of trajectories"
    continue    # Skip to next file
  elif L < 0 or vol < 0:
    print infile, "never defines lattice volume"
    print >> ERRFILE, infile, "never defines lattice volume"
    continue    # Skip to next file

  if (traj == 0 and int(load) > 0) or (int(load) != oldcfg):
    print "Guessing approximate starting trajectory for", infile
    traj = int(load)
    endtraj = traj + traj_per_file
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Cycle through lines in the "out" file
  start_traj = 0      # In case we started from a unit configuration
  oldcfg = int(cfg)
  for line in open(infile):
    if line.startswith('traj_length'):
      tlength = float((line.split())[1])
      print >> TLENGTH, "%d,%g" % (endtraj, tlength)
    elif line.startswith('nsquares'):
      cpus = 1
      for i in range(-4, 0):    # range doesn't include last
        cpus *= int((line.split())[i])
    elif line.startswith('Hasenbusch_mass'):
      print >> MH, "%d,%g" % (endtraj, float((line.split())[1]))

    # Check that the file loaded the appropriate configuration
    elif line.startswith('plaq'):
      if stamp == "start":    # Loading configuration
        stamp = line.rstrip()
        if stamp != oldstamp and oldstamp != "start":
          print infile, "plaq doesn't match final", oldstamp
          print >> ERRFILE, infile, "plaq doesn't match final", oldstamp
      else:                   # Saving configuration
        oldstamp = line.rstrip()    # Prepare for next file
        stamp = "start"
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Now extract evolution observables
    # Acceptance comes first
    elif line.startswith('ACCEPT') or line.startswith('REJECT'):
      traj += 1
      MDTU += tlength
      print >> TU, "%d,%g" % (traj, MDTU)
      print >> KEY, "%g,%s" % (traj, tag)
      print >> DELTAS, "%d,%g" % (traj, float((line.split())[3]))
      print >> ABS_DS, "%d,%g" % (traj, abs(float((line.split())[3])))
      print >> EXP_DS, "%d,%g" % (traj, float((line.split())[6]))
      if line.startswith('ACCEPT'):
        print >> ACCP, "%d,1" % traj
      else:
        print >> ACCP, "%d,0" % traj

    elif line.startswith('loading lattice:'):
      start_traj = int((line.split('.'))[-2])

    elif line.startswith('traj '):
      print >> WALLTIME, "%d,%g" % (traj, float((line.split())[3]))
      temp = start_traj + int((line.split())[1])
      if traj != temp:
        print "Trajectory count mismatch in", infile,
        print "%d vs %d" % (traj, temp)
        print >> ERRFILE, "Trajectory count mismatch in", infile,
        print >> ERRFILE, "%d vs %d" % (traj, temp)

    # Forces and nsteps -- monitor maxima rather than rms...
    # Annoying dependence on Nf...
    elif line.startswith('GF'):
      Nstep_gauge = int((line.split())[5])
      stepsize_gauge = tlength / float(Nstep_gauge)
      force_gauge = float((line.split())[11])
    elif line.startswith('FF[01]'):   # Robust to Nf
      Nstep = int((line.split())[5])
      stepsize = tlength / float(Nstep)
      force0 = float((line.split())[11])
    elif line.startswith('FFtot'):    # Robust to Nf
      Nstep1 = int((line.split())[5]) - Nstep
      stepsize1 = tlength / float(Nstep1)
      print >> FORCE, "%d,%g,%g,%g" % (traj, force0, force1, force_gauge)
      print >> NSTEP, "%d,%d,%d,%d" % (traj, Nstep, Nstep1, Nstep_gauge)
      print >> STEPSIZE, "%d,%g,%g,%g" \
                         % (traj, stepsize, stepsize1, stepsize_gauge)
    elif line.startswith(Ftag):       # Nf-dependent Ftag set above
      force1 = float((line.split())[11])

    # Could split CG iterations into inner and outer steps
    # but just grab total for now
    elif line.startswith('CGtot'):
      iters = int((line.split())[5]) * int((line.split())[9])
      print >> CG_ITERS, "%d,%d" % (traj, iters)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Now extract physical observables and wrap up "out" file
    # Use consistent normalization for plaquette and pbp
    elif line.startswith('MEASplaq'):
      plaq_ss = Nc * float((line.split())[2])
      plaq_st = Nc * float((line.split())[4])
      print >> PLAQ, "%g,%g,%g" % (MDTU, plaq_ss, plaq_st)
    elif line.startswith('MEASploop'):
      poly_r = float((line.split())[5])
      poly_i = float((line.split())[6])
      print >> POLY, "%g,%g" % (poly_r, poly_i)
      print >> POLY_R, "%g,%g" % (MDTU, poly_r)
      poly_mod = math.sqrt(poly_r**2 + poly_i**2)
      print >> POLY_MOD, "%g,%g" % (MDTU, poly_mod)
      poly_arg = math.atan2(poly_i, poly_r)
      print >> POLY_ARG, "%g,%g" % (MDTU, poly_arg)

      # Spatial Wilson line
      poly_r = float((line.split())[2])
      poly_i = float((line.split())[3])
      print >> XPOLY, "%g,%g" % (poly_r, poly_i)
      print >> XPOLY_R, "%g,%g" % (MDTU, poly_r)
      poly_mod = math.sqrt(poly_r**2 + poly_i**2)
      print >> XPOLY_MOD, "%g,%g" % (MDTU, poly_mod)
      poly_arg = math.atan2(poly_i, poly_r)
      print >> XPOLY_ARG, "%g,%g" % (MDTU, poly_arg)

    # QHMC measures pbp after every trajectory
    # Just use a single stochastic source
    elif line.startswith('MEASpbp'):
      print >> PBP, "%g,%g" % (MDTU, float((line.split())[4]) / 2.0)

    elif line.startswith('total time:'):
      walltime = float((line.split())[2])
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Check to see if run seems to have finished properly
  if walltime < 0:
    print infile, "didn't print final timing"
    print >> ERRFILE, infile, "didn't print final timing"
  else:   # We are good to go
    TUtime = walltime * float(cpus) / (60.0 * tlength * traj_per_file)
    print >> WALLTU, "%d,%g" % (traj, TUtime)
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Now deal with the corresponding "Sout" file, if it is present
  infile = 'Out/Sout.' + cfg
  if os.path.isfile(infile):
    pd = []   # Plaquette differences
    ld = []   # Link differences
    check = -1
    for line in open(infile):
      # Try to check plaquette...
      if line.startswith('CHECK PLAQ: '):
        order_check = float((line.split())[2]) + float((line.split())[3])
        out_check = 2.0 * Nc * float((oldstamp.split())[-1])
        diff = abs(order_check - out_check)
        if diff > 1e-5:
          print infile, "starting plaquette doesn't match saved:",
          print "|%g - %g| = %g" % (order_check, out_check, diff)
          print >> ERRFILE, infile, "starting plaquette doesn't match saved:",
          print >> ERRFILE, "|%g - %g| = %g" % (order_check, out_check, diff)

      elif line.startswith('StaggPlaq t '):
        pd.append(float((line.split())[2]) - float((line.split())[3]))
      elif line.startswith('StaggPlaq x '):
        pd.append(float((line.split())[2]) - float((line.split())[3]))
      elif line.startswith('StaggPlaq y '):
        pd.append(float((line.split())[2]) - float((line.split())[3]))
      elif line.startswith('StaggPlaq z '):
        pd.append(float((line.split())[2]) - float((line.split())[3]))
      elif line.startswith('pbpt: '):
        ld.append(float((line.split())[3]) - float((line.split())[4]))
      elif line.startswith('pbpx: '):
        ld.append(float((line.split())[3]) - float((line.split())[4]))
      elif line.startswith('pbpy: '):
        ld.append(float((line.split())[3]) - float((line.split())[4]))
      elif line.startswith('pbpz: '):
        ld.append(float((line.split())[3]) - float((line.split())[4]))

      elif line.startswith('RUNNING COMPLETED'):
        check = 1

    # Post-processing and printing
    if check < 0:
      print infile, "did not complete"
      print >> ERRFILE, infile, "did not complete"

    if len(pd) != 4 or len(ld) != 4:
      print "Measured %d and %d differences in %s" \
            % (len(pd), len(ld), infile)
    else:
      pd.append(0.0)
      ld.append(0.0)
      for i in range(4):
        pd[4] += pd[i]**2
        ld[4] += ld[i]**2
      pd[4] = math.sqrt(pd[4])
      ld[4] = math.sqrt(ld[4])
      print >> PLAQ_DIFF, "%g,%g,%g,%g,%g,%g" \
                          % (MDTU, pd[0], pd[1], pd[2], pd[3], pd[4])
      print >> LINK_DIFF, "%g,%g,%g,%g,%g,%g" \
                          % (MDTU, ld[0], ld[1], ld[2], ld[3], ld[4])
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Now deal with the corresponding "eig" file, if it is present
  infile = 'Out/eig.' + cfg
  if os.path.isfile(infile):
    print "ERROR: haven't yet set up processing for", infile
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Now deal with the corresponding "Wflow" file, if it is present
  # Print gSq for c=0.2, 0.25, 0.3 and 0.35
  # Print topological charge for c=0.2, 0.3, 0.4 and 0.5
  # These files now contain MCRG-blocked results: only consider t=0
  # Also ignore finite-volume correction in definition of gSq
  infile = 'Out/Wflow.' + cfg
  if not os.path.isfile(infile):
    print >> MISSINGFILES, infile
  else:
    FNAL = -1
    check = -1
    c = -1.0    # sqrt(8t) / L
    cOld = 1.0
    gSq = []
    topo = []   # Topological charge
    poly = []   # Wilson-flowed Polyakov loop

    # RG-blocked observables
    plaq_blocked = -1     # Check to see how much to print
    poly_blocked = -1
    plaq = ["null", "null", "null", "null", "null"]   # Initialize output
    p_r = ["null", "null", "null", "null", "null"]    # Polyakov loop
    p_i = ["null", "null", "null", "null", "null"]
    x_r = ["null", "null", "null", "null", "null"]    # Spatial loop
    x_i = ["null", "null", "null", "null", "null"]
    p_mod = ["null", "null", "null", "null", "null"]
    p_arg = ["null", "null", "null", "null", "null"]
    x_mod = ["null", "null", "null", "null", "null"]
    x_arg = ["null", "null", "null", "null", "null"]

    # Go
    for line in open(infile):
      # See where we ran this (annoying legacy normalization issue at Fermilab)
      if 'lqcdproj' in line:
        FNAL = 1
      # Try to check plaquette...
      elif line.startswith('CHECK PLAQ: '):
        order_check = float((line.split())[2]) + float((line.split())[3])
        out_check = 2.0 * 3.0 * float((oldstamp.split())[-1])
        diff = abs(order_check - out_check)
        if diff > 1e-5:
          print infile, "starting plaquette doesn't match saved:",
          print "|%g - %g| = %g" % (order_check, out_check, diff)
          print >> ERRFILE, infile, "starting plaquette doesn't match saved:",
          print >> ERRFILE, "|%g - %g| = %g" % (order_check, out_check, diff)

      elif line.startswith('WFLOW '):
        c = math.sqrt(8.0 * float((line.split())[1])) / L
        if cOld < 0.2 and c >= 0.2:
          gSq.append(gprop * float((line.split())[4]))
          topo.append(float((line.split())[7]))
        elif cOld < 0.25 and c >= 0.25:
          gSq.append(gprop * float((line.split())[4]))
        elif cOld < 0.3 and c >= 0.30:
          gSq.append(gprop * float((line.split())[4]))
          topo.append(float((line.split())[7]))
        elif cOld < 0.35 and c >= 0.35:
          gSq.append(gprop * float((line.split())[4]))
        elif cOld < 0.4 and c >= 0.40:
          topo.append(float((line.split())[7]))
        elif cOld < 0.5 and c >= 0.50:
          topo.append(float((line.split())[7]))
        cOld = c

      # Wilson-flowed Polyakov loop
      elif line.startswith('POLYA ORIG '):
        temp = line.split()
        cp = math.sqrt(8.0 * float((line.split())[2])) / L
        if cp == 0.0:
          p_r[0] = float(temp[3])
          p_i[0] = float(temp[4])
          x_r[0] = float(temp[5])
          x_i[0] = float(temp[6])
        # Should pick up cp=0.2, 0.3, 0.4 and 0.5 in that order
        elif 0.19 < cp:
          poly.append(float(temp[3]))

      # RG-blocked observables (with no Wilson flow)
      elif line.startswith('LOOPS 0 '):
        temp = line.split()
        # Should pick up unblocked then each blocking level in order
        if int(temp[2]) == 0:
          # Check that the blocking level is correct
          plaq_blocked = int(temp[4])
          plaq[plaq_blocked] = float(temp[6])
      elif line.startswith('POLYA NHYP 0 '):
        temp = line.split()
        # Should pick up each blocking level in order
        if temp[4] == "0.6":
          poly_blocked = int(temp[3])
          p_r[poly_blocked] = float(temp[5])
          p_i[poly_blocked] = float(temp[6])
          x_r[poly_blocked] = float(temp[7])
          x_i[poly_blocked] = float(temp[8])

      elif line.startswith('RUNNING COMPLETED'):
        check = 1

    # Post-processing and printing
    if check < 0:
      print infile, "did not complete"
      print >> ERRFILE, infile, "did not complete"

    # Annoying legacy normalization issue in Fermilab measurements
    if FNAL > 0:
      for i in range(len(topo)):
        topo[i] *= vol * 0.02533029591058444286**2    # 1/4pi^2

    print >> WFLOW, "%g,%g,%g,%g,%g" % (MDTU, gSq[0], gSq[1], gSq[2], gSq[3])
    print >> TOPO, "%g,null,null,null,%g" % (MDTU, topo[3])
    if len(poly) == 4:
      print >> WPOLY, "%g,%g,%g,%g,%g" \
                      % (MDTU, poly[0], poly[1], poly[2], poly[3])

    # Lots of RG-blocked stuff to print
    if plaq_blocked > 0:
      toprint = str(MDTU) + ','
      for i in range(4):
        toprint += str(plaq[i]) + ','
      toprint += str(plaq[4])
      print >> PLAQB, toprint

    if poly_blocked > 0:
      for i in range(poly_blocked + 1):
        p_mod[i] = math.sqrt(p_r[i]**2 + p_i[i]**2)
        p_arg[i] = math.atan2(p_i[i], p_r[i])
        x_mod[i] = math.sqrt(x_r[i]**2 + x_i[i]**2)
        x_arg[i] = math.atan2(x_i[i], x_r[i])
      print >> POLYB, "%g,%g,null,null,null,null" % (p_r[0], p_i[0])
      print >> POLYB, "%g,null,%g,null,null,null" % (p_r[1], p_i[1])
      print >> XPOLYB, "%g,%g,null,null,null,null" % (x_r[0], x_i[0])
      print >> XPOLYB, "%g,null,%g,null,null,null" % (x_r[1], x_i[1])
      if poly_blocked >= 2:
        print >> POLYB, "%g,null,null,%g,null,null" % (p_r[2], p_i[2])
        print >> XPOLYB, "%g,null,null,%g,null,null" % (x_r[2], x_i[2])
      if poly_blocked >= 3:
        print >> POLYB, "%g,null,null,null,%g,null" % (p_r[3], p_i[3])
        print >> XPOLYB, "%g,null,null,null,%g,null" % (x_r[3], x_i[3])
      if poly_blocked >= 4:
        print >> POLYB, "%g,null,null,null,null,%g" % (p_r[4], p_i[4])
        print >> XPOLYB, "%g,null,null,null,null,%g" % (x_r[4], x_i[4])

      p_r_out = str(MDTU) + ','
      p_mod_out = str(MDTU) + ','
      p_arg_out = str(MDTU) + ','
      x_r_out = str(MDTU) + ','
      x_mod_out = str(MDTU) + ','
      x_arg_out = str(MDTU) + ','
      for i in range(4):
        p_r_out += str(p_r[i]) + ','
        p_mod_out += str(p_mod[i]) + ','
        p_arg_out += str(p_arg[i]) + ','
        x_r_out += str(x_r[i]) + ','
        x_mod_out += str(x_mod[i]) + ','
        x_arg_out += str(x_arg[i]) + ','
      p_r_out += str(p_r[4])
      p_mod_out += str(p_mod[4])
      p_arg_out += str(p_arg[4])
      x_r_out += str(x_r[4])
      x_mod_out += str(x_mod[4])
      x_arg_out += str(x_arg[4])
      print >> POLY_RB, p_r_out
      print >> POLY_MODB, p_mod_out
      print >> POLY_ARGB, p_arg_out
      print >> XPOLY_RB, x_r_out
      print >> XPOLY_MODB, x_mod_out
      print >> XPOLY_ARGB, x_arg_out
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Clean up and close down
ERRFILE.close()
MISSINGFILES.close()
PLAQ.close()
PBP.close()
POLY.close()
XPOLY.close()
POLY_R.close()
POLY_MOD.close()
XPOLY_R.close()
XPOLY_MOD.close()
POLY_ARG.close()
XPOLY_ARG.close()
PLAQ_DIFF.close()
LINK_DIFF.close()
EIG.close()
WFLOW.close()
TOPO.close()
WPOLY.close()
PLAQB.close()
POLYB.close()
XPOLYB.close()
POLY_RB.close()
POLY_MODB.close()
XPOLY_RB.close()
XPOLY_MODB.close()
POLY_ARGB.close()
XPOLY_ARGB.close()
ACCP.close()
EXP_DS.close()
DELTAS.close()
ABS_DS.close()
FORCE.close()
CG_ITERS.close()
WALLTIME.close()
WALLTU.close()
NSTEP.close()
STEPSIZE.close()
MH.close()
TLENGTH.close()
KEY.close()
TU.close()
# ------------------------------------------------------------------
