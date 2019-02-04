#!/usr/bin/env python
import os
import sys
import numpy as np
# ------------------------------------------------------------------
# Compute latent heat and related quantities from dygraph data files
# Combine both high and low starts for the given parameters
# Will check aspect ratios alpha=2, 3 and 4, so do so all at once

# Parse arguments: First whether we're considering poly or Wpoly
# Then Nt, beta, mass to specify the ensembles to consider,
# and the six corresponding thermalization cuts, high then low for each L,
# and finally the separation between deconfined vs. confined
if len(sys.argv) < 11:
  print "Usage:", str(sys.argv[0]), "<file> <Nt> <beta> <mass>",
  print "<a2 high cut> <a2 low cut>",
  print "<a3 high cut> <a3 low cut>",
  print "<a4 high cut> <a4 low cut> <boundary>",
  sys.exit(1)

# poly_file should be either poly_mod or Wpoly_mod
# It will determine which of the plaquette vs. clover energy we loo at...
poly_file = str(sys.argv[1])
if poly_file.startswith('W'):
  print "ERROR: Haven't yet set up Wilson-flowed analysis..."
  sys.exit(1)
else:
  plaq_file = 'plaq'

Nt = int(sys.argv[2])
beta = str(sys.argv[3])
mass = str(sys.argv[4])
high_cuts = [int(sys.argv[5]), int(sys.argv[7]), int(sys.argv[9])]
low_cuts = [int(sys.argv[6]), int(sys.argv[8]), int(sys.argv[10])]

div = float(sys.argv[11])

# Set up directories --- for now hard-code Nt=6, 8, three aspect ratios
Ndirs = 3
if Nt == 6:
  L = ['L=12', 'L=18', 'L=24']      # Labels
  tag = beta + '_high_m' + mass
  high_dirs = ['12nt6/b' + tag, '18nt6/b' + tag, '24nt6/b' + tag]
  tag = beta + '_low_m' + mass
  low_dirs = ['12nt6/b' + tag, '18nt6/b' + tag, '24nt6/b' + tag]
elif Nt == 8:
  L = ['L=16', 'L=24', 'L=32']      # Labels
  tag = beta + '_high_m' + mass
  high_dirs = ['16nt8/b' + tag, '24nt8/b' + tag, '32nt8/b' + tag]
  tag = beta + '_low_m' + mass
  low_dirs = ['16nt8/b' + tag, '24nt8/b' + tag, '32nt8/b' + tag]
else:
  print "ERROR: Only Nt=6 and 8 set up at the moment"
  sys.exit(1)

# Check that all data files exist and Nc=4
for dirname in high_dirs + low_dirs:
  tocheck = dirname + '/data/' + poly_file + '.csv'
  if not os.path.isfile(tocheck):
    print "ERROR:", tocheck, "not found"
    sys.exit(1)

cwd = os.getcwd()
if 'SU4' in cwd:
  Nc = 4.0
  eps_SB = np.pi**2
else:
  print "ERROR: Only expecting SU4 at the moment"
  sys.exit(1)

# Constant normalization factors from tree-level perturbation theory
# Omit factors of Nt^4 until final combination
# Plugging in beta_func=-11/(6pi^2) and g^{-2}=beta_F/(4N)
norm_minus = 11.0 * Nc * np.power(np.pi, -2)
Karsch = 1.0 - 0.236626 * 4.0 * Nc / float(beta)
norm_plus = 2.0 * float(beta) * Karsch

# Quick check
#print norm_minus, norm_plus

# Load Polyakov loop and plaquette data
# Divide latter into deconfined vs. confined based on div
# Normalize poly to be out of Nc, plaq to be out of 1
# poly format: MDTU,poly_mod              (normalized to 1)
# Wpoly format: MDTU,c=0.2,0.3,0.4,0.5    (normalized to Nc)
# plaq format: MDTU,plaq_ss,plaq_st       (normalized to Nc)
for i in range(Ndirs):
  poly_dat = []
  sum_dec = []    # Deconfined plaq_st + plaq_ss
  sum_con = []    # Confined   plaq_st + plaq_ss
  dif_dec = []    # Deconfined plaq_st - plaq_ss
  dif_con = []    # Confined   plaq_st - plaq_ss

  # High start
  toOpen = high_dirs[i] + '/data/' + poly_file + '.csv'
  for line in open(toOpen):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= high_cuts[i]:
      continue
    if poly_file.startswith('W'):
      poly_dat.append(float(temp[4]))
    else:
      poly_dat.append(Nc * float(temp[1]))

  count = 0               # Sync plaq with poly...
  toOpen = high_dirs[i] + '/data/' + plaq_file + '.csv'
  for line in open(toOpen):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= high_cuts[i]:
      continue
    if poly_file.startswith('W'):
      print "ERROR: Haven't yet set up Wilson-flowed analysis..."
      sys.exit(1)
    else:
      plaq_ss = float(temp[1]) / Nc
      plaq_st = float(temp[2]) / Nc
      if poly_dat[count] > div:
        sum_dec.append(plaq_st + plaq_ss)
        dif_dec.append(plaq_st - plaq_ss)
      else:
        sum_con.append(plaq_st + plaq_ss)
        dif_con.append(plaq_st - plaq_ss)
      count += 1

  # Low start
  toOpen = low_dirs[i] + '/data/' + poly_file + '.csv'
  for line in open(toOpen):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= low_cuts[i]:
      continue
    if poly_file.startswith('W'):
      poly_dat.append(float(temp[4]))
    else:
      poly_dat.append(Nc * float(temp[1]))

  toOpen = low_dirs[i] + '/data/' + plaq_file + '.csv'
  for line in open(toOpen):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= low_cuts[i]:
      continue
    if poly_file.startswith('W'):
      print "ERROR: Haven't yet set up Wilson-flowed analysis..."
      sys.exit(1)
    else:
      plaq_ss = float(temp[1]) / Nc
      plaq_st = float(temp[2]) / Nc
      if poly_dat[count] > div:
        sum_dec.append(plaq_st + plaq_ss)
        dif_dec.append(plaq_st - plaq_ss)
      else:
        sum_con.append(plaq_st + plaq_ss)
        dif_con.append(plaq_st - plaq_ss)
      count += 1

  # Sanity check
  temp = len(sum_dec) + len(sum_con)
  if not len(poly_dat) == temp:
    print "ERROR: %d poly meas vs. %d plaq meas..." % (len(poly_dat), temp)
    sys.exit(1)

  # Print number of measurements to allow offline checks
  print "%s\n%d = %d + %d measurements" \
        % (L[i], len(poly_dat), len(sum_dec), len(sum_con))

  # Average plaquette data and extract quantities of interest
  # I think this should just be a container for pointers...
  lists = [sum_con, sum_dec, dif_con, dif_dec]
  aves = np.empty(len(lists), dtype = np.float)
  errs = np.empty_like(aves)
  for j in range(len(lists)):
    dat = np.array(lists[j])
    N = np.size(dat)
    if N == 0:
      aves[j] = np.nan
      errs[j] = np.nan
    else:
      aves[j] = np.mean(dat, dtype = np.float64)
      errs[j] = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
  sum_con_ave, sum_dec_ave, dif_con_ave, dif_dec_ave = aves
  sum_con_err, sum_dec_err, dif_con_err, dif_dec_err = errs

  # Print a^4(e-3p) and a^4(e+p) as in Wingate and Ohta Figs. 12 and 13
  em3p_dec = norm_minus * sum_dec_ave
  em3p_dec_err = norm_minus * sum_dec_err
  epp_dec = norm_plus * dif_dec_ave
  epp_dec_err = norm_plus * dif_dec_err
  print "Deconfined a^4(e - 3p): %.6g %.4g" % (em3p_dec, em3p_dec_err)
  print "Deconfined a^4(e +  p): %.6g %.4g" % (epp_dec, epp_dec_err)
  em3p_con = norm_minus * sum_con_ave
  em3p_con_err = norm_minus * sum_con_err
  epp_con = norm_plus * dif_con_ave
  epp_con_err = norm_plus * dif_con_err
  print "  Confined a^4(e - 3p): %.6g %.4g" % (em3p_con, em3p_con_err)
  print "  Confined a^4(e +  p): %.6g %.4g" % (epp_con, epp_con_err)

  # Deconfined - confined differences
  De_minus = em3p_dec - em3p_con
  De_minus_err = np.sqrt(em3p_dec_err**2 + em3p_con_err**2)
  De_plus = epp_dec - epp_con
  De_plus_err = np.sqrt(epp_dec_err**2 + epp_con_err**2)

  # Latent heat and pressure change, restoring Nt^4
  De_eps = 0.25 * Nt**4 * (3.0 * De_plus + De_minus)
  De_p = 0.25 * Nt**4 * (De_plus - De_minus)
  De_eps_err = 0.25 * Nt**4 * np.sqrt(9.0 * De_plus_err**2 + De_minus_err**2)
  De_p_err = 0.25 * Nt**4 * np.sqrt(De_plus_err**2 + De_minus_err**2)

  print "De_eps %.6g %.4g" % (De_eps, De_eps_err)
  print "SBnorm %.6g %.4g" % (De_eps / eps_SB, De_eps_err / eps_SB)
  print "De_p   %.6g %.4g" % (De_p, De_p_err)

  # TODO: Can probably still improve this workflow for gnuplotting
  print "%s %.6g %.4g %.6g %.4g" \
        % (beta, em3p_dec, em3p_dec_err, epp_dec, epp_dec_err)
  print "%s %.6g %.4g %.6g %.4g\n" \
        % (beta, em3p_con, em3p_con_err, epp_con, epp_con_err)
# ------------------------------------------------------------------
