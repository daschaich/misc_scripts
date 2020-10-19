#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
import h5py
# ------------------------------------------------------------------
# Package SU(4) data and results/attributes into HDF5 file

# Cycle over all streams and write to ~/LSD/SU4/SU4_data.hdf5
# Group paths will specify Nf, Nt, L, beta_F, (valence) mass and start
# Attributes for each stream:
#   Number of trajectories, acceptance rate, (Wpoly) autocorrelation time,
#   thermalization cut, block size, number of blocks,
#   and pbp autocorrelation time for reference
# Datasets for each stream, each with ave, err (suscept) as attributes:
#   plaquette, |PL|, arg(PL), real(PL), chiral condensate, exp(-Delta S),
#   list of Wflow measurements, Wflow_aniso, |PL_W|, arg(PL_W), real(PL_W)
# arg(PL) and arg(PL_W) not averaged, but determine deconfinement fraction
# Combined deconf. frac. and susceptibilities as attributes of combo group
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Annoyingly, need to use absolute Barkla path
path = "/users/schaich/LSD/SU4/"
if not os.path.isdir(path):
  print("ERROR: " + path + " not found")
  sys.exit(1)
os.chdir(path)
f = h5py.File("SU4_data.h5", 'w')

# Top-level groups for each Nf
for Nf in glob.glob('*f'):
  f.create_group(Nf)

  # Second- and third-level groups for each (Nt, L) volume
  os.chdir(path + Nf)
  for vol in glob.glob('*nt*'):
    temp = vol.split('nt')
    L = 'L' + temp[0]
    Nt = 'Nt' + temp[1]
    this_vol = Nf + '/' + Nt + '/' + L
    f.create_group(this_vol)

    # Fourth- and fifth-level groups for each (beta_F, mass) ensemble
    # For 0f ensembles this is the valence mass for pbp
    # Then sixth-level groups for each stream ('high', 'low', 'combo')
    os.chdir(path + Nf + '/' + vol)
    for stream in glob.glob('b[1-9]*'):   # Avoid 'backup' dirs
      os.chdir(path + Nf + '/' + vol + '/' + stream)
      toCheck = 'results/Wpoly_mod.autocorr'
      if not os.path.isfile(toCheck):     # Skip unfinished streams
        print("Skipping %s" % Nf + '/' + vol + '/' + stream)
        continue

      temp = stream.split('_')
      if vol == '24nt48':                 # Simpler format
        beta = temp[0]
        start = 'placeholder'
        mass = temp[1]
        this_str = this_vol + '/' + mass + '/' + beta
        this_grp = f.create_group(this_str)
      else:
        beta = temp[0]
        start = temp[1]
        mass = temp[2]
        this_ens = this_vol + '/' + mass + '/' + beta
        this_str = this_ens + '/' + start
        this_grp = f.create_group(this_str)

      # Record thermalization cut and block size
      check = -1
      therm = path + Nf + '/' + vol + '/therm.sh'
      for line in open(therm):
        if stream in line:
          temp = line.split()
          this_grp.attrs['thermalization cut'] = int(temp[2])
          this_grp.attrs['block size'] = int(temp[3])
          check = int(temp[2])
          break
      if check < 0:
        print("ERROR: Thermalization cut not found for %s: " % this_str)
        sys.exit(1)

      # ------------------------------------------------------------
      # Record acceptance, auto-correlation and Nblocks for stream
      for line in open('results/accP.dat'):
        temp = line.split()
        Nblocks = int(temp[-1])
        this_grp.attrs['Nblocks'] = Nblocks
        this_grp.attrs['acceptance'] = float(temp[0])
      for line in open('results/Wpoly_mod.autocorr'):
        temp = line.split()
        this_grp.attrs['autocorrelation time'] = float(temp[0])
      for line in open('results/pbp.autocorr'):
        if not line.startswith('W'):      # Ignore warning message
          temp = line.split()
          this_grp.attrs['autocorrelation time (pbp)'] = float(temp[0])

      # Now set up data and attributes/results for this stream
      # First do plaq, averaging ss & st and also setting Ntraj
      plaq_arr = []
      for line in open('data/plaq.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        dat = 0.5 * (float(temp[1]) + float(temp[2]))
        plaq_arr.append(dat)
      Ntraj = len(plaq_arr)
      this_grp.attrs['Ntraj'] = Ntraj

      plaq = np.array(plaq_arr)
      dset = this_grp.create_dataset('plaq', data=plaq)

      # Record results as attributes
      for line in open('results/plaq.dat'):
        temp = line.split()
        dset.attrs['ave'] = float(temp[0])
        dset.attrs['err'] = float(temp[1])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
          print("%s vs %d in plaq.dat" % (temp[-1], Nblocks))
          sys.exit(1)
      for line in open('results/plaq.suscept'):
        temp = line.split()
        dset.attrs[temp[0]] = float(temp[1])
        dset.attrs[temp[0] + '_err'] = float(temp[2])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
          print("%s vs %d in plaq.suscept" % (temp[-1], Nblocks))
          sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Next do simple observables measured every trajectory=MDTU
      # All these have a single datum per line
      # Skipping cg_iters, xpoly, wall_time
      for obs in ['poly_mod', 'poly_arg', 'poly_r', 'pbp', 'exp_dS']:
        obsfile = 'data/' + obs + '.csv'
        dset = this_grp.create_dataset(obs, (Ntraj,), dtype='f')

        traj = 0
        for line in open(obsfile):
          temp = line.split(',')
          if line.startswith('M') or line.startswith('t'):
            continue
          dset[traj] = float(temp[1])
          traj += 1

        if traj != Ntraj:
          print("ERROR: Ntraj mismatch in %s, " % this_str, end='')
          print("%d vs %d in %s" % (traj, Ntraj, obsfile))
          sys.exit(1)

        # Results as attributes, checking Nblocks
        resfile = 'results/' + obs + '.dat'
        if os.path.isfile(resfile):           # No ave for poly_arg
          for line in open(resfile):
            if line.startswith('#'):
              continue
            temp = line.split()
            dset.attrs['ave'] = float(temp[0])
            dset.attrs['err'] = float(temp[1])
            if not int(temp[-1]) == Nblocks:
              print("ERROR: Nblocks mismatch in %s, " % this_str, end='')
              print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
              sys.exit(1)

        # Same for susceptibility, skewness, kurtosis if present
        suscfile = 'results/' + obs + '.suscept'
        if os.path.isfile(suscfile):        # No suscept for some obs
          for line in open(suscfile):
            temp = line.split()
            dset.attrs[temp[0]] = float(temp[1])
            dset.attrs[temp[0] + '_err'] = float(temp[2])
            if not int(temp[-1]) == Nblocks:
              print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
              print("%s vs %d in %s.suscept" % (temp[-1], Nblocks, obs))
              sys.exit(1)

        # Deconfinement fractions as poly_arg attributes
        if obs == 'poly_arg':
          for theta in ['0.1', '0.15', '0.2', '0.25', '0.3']:
            fracfile = 'results/deconf_frac-' + theta + '.dat'
            if not os.path.isfile(fracfile):
              continue
            for line in open(fracfile):
              temp = line.split()
              dset.attrs['deconf frac theta=' + theta] = float(temp[0])
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Now for Wilson-flowed observables
      # First go through topo.dat to set Nmeas
      # Save trajectories at which Wilson flow measurements were run
      # Only have c=0.5 topological charge (for clean time-series plot)
      # We don't average the topological charge, so it has no attributes
      topo_arr = []
      Wmeas_arr = []
      for line in open('data/topo.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        Wmeas_arr.append(int(temp[0]))
        topo_arr.append(float(temp[4]))
      Nmeas = len(Wmeas_arr)
      Wmeas = np.array(Wmeas_arr)
      topo = np.array(topo_arr)
      this_grp.create_dataset('Wmeas', data=Wmeas)
      dset = this_grp.create_dataset('topo', data=topo)
      dset.attrs['columns'] = ['c=0.5']

      # Also treat Wpoly_r separately, because of naming issue
      dset = this_grp.create_dataset('Wpoly_r', (Nmeas,4,), dtype='f')
      dset.attrs['columns'] = ['c=0.2', 'c=0.3', 'c=0.4', 'c=0.5']
      meas = 0
      for line in open('data/Wpoly.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        dset[meas] = [float(temp[1]), float(temp[2]), \
                      float(temp[3]), float(temp[4])]
        meas += 1

      if meas != Nmeas:
        print("ERROR: Nmeas mismatch in %s, " % this_str, end='')
        print("%d vs %d in Wpoly_r" % (meas, Nmeas))
        sys.exit(1)

      # Results as attributes, checking Nblocks
      for line in open('results/Wpoly.dat'):
        if line.startswith('#'):
          continue
        temp = line.split()
        dset.attrs['c=0.2 ave'] = float(temp[0])
        dset.attrs['c=0.2 err'] = float(temp[1])
        dset.attrs['c=0.3 ave'] = float(temp[2])
        dset.attrs['c=0.3 err'] = float(temp[3])
        dset.attrs['c=0.4 ave'] = float(temp[4])
        dset.attrs['c=0.4 err'] = float(temp[5])
        dset.attrs['c=0.5 ave'] = float(temp[6])
        dset.attrs['c=0.5 err'] = float(temp[7])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s, " % this_str, end='')
          print("%s vs %d in Wpoly_r" % (temp[-1], Nblocks))
          sys.exit(1)

      # Same for susceptibility, skewness and kurtosis,
      # except that we only compute these for c=0.5
      for line in open('results/Wpoly.suscept'):
        temp = line.split()
        dset.attrs['c=0.5 ' + temp[0]] = float(temp[1])
        dset.attrs['c=0.5 ' + temp[0] + '_err'] = float(temp[2])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
          print("%s vs %d in Wpoly_r.suscept" % (temp[-1], Nblocks))
          sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Move on to other Wilson-flowed observables
      # All measured for c=0.2, 0.3, 0.4 and 0.5
      for obs in ['Wpoly_mod', 'Wpoly_arg', 'Wflow_aniso']:
        obsfile = 'data/' + obs + '.csv'
        if not os.path.isfile(obsfile):       # Only Nt=8 has Wpoly_arg
          continue

        dset = this_grp.create_dataset(obs, (Nmeas,4), dtype='f')
        dset.attrs['columns'] = ['c=0.2', 'c=0.3', 'c=0.4', 'c=0.5']
        meas = 0
        for line in open(obsfile):
          if line.startswith('M'):
            continue
          temp = line.split(',')
          dset[meas] = [float(temp[1]), float(temp[2]), \
                        float(temp[3]), float(temp[4])]
          meas += 1

        if meas != Nmeas:
          print("ERROR: Nmeas mismatch in %s, " % this_str, end='')
          print("%d vs %d in %s" % (meas, Nmeas, obsfile))
          sys.exit(1)

        # Results as attributes, checking Nblocks
        resfile = 'results/' + obs + '.dat'
        if os.path.isfile(resfile):           # No ave for Wpoly_arg
          for line in open(resfile):
            if line.startswith('#'):
              continue
            temp = line.split()
            dset.attrs['c=0.2 ave'] = float(temp[0])
            dset.attrs['c=0.2 err'] = float(temp[1])
            dset.attrs['c=0.3 ave'] = float(temp[2])
            dset.attrs['c=0.3 err'] = float(temp[3])
            dset.attrs['c=0.4 ave'] = float(temp[4])
            dset.attrs['c=0.4 err'] = float(temp[5])
            dset.attrs['c=0.5 ave'] = float(temp[6])
            dset.attrs['c=0.5 err'] = float(temp[7])
            if not int(temp[-1]) == Nblocks:
              print("ERROR: Nblocks mismatch in %s, " % this_str, end='')
              print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
              sys.exit(1)

        # Same for susceptibility, skewness, kurtosis if present,
        # except that we only compute these for c=0.5
        suscfile = 'results/' + obs + '.suscept'
        if os.path.isfile(suscfile):        # No suscept for some obs
          for line in open(suscfile):
            temp = line.split()
            dset.attrs['c=0.5 ' + temp[0]] = float(temp[1])
            dset.attrs['c=0.5 ' + temp[0] + '_err'] = float(temp[2])
            if not int(temp[-1]) == Nblocks:
              print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
              print("%s vs %d in %s.suscept" % (temp[-1], Nblocks, obs))
              sys.exit(1)

        # Deconfinement fractions as Wpoly_arg attributes
        if obs == 'Wpoly_arg':
          for theta in ['0.1', '0.15', '0.2', '0.25', '0.3']:
            fracfile = 'results/Wdeconf_frac-' + theta + '.dat'
            if not os.path.isfile(fracfile):
              continue
            tag = 'Wdeconf frac theta=' + theta
            for line in open(fracfile):
              if line.startswith('#'):
                continue
              temp = line.split()
              dset.attrs[tag + ' c=0.2'] = float(temp[0])
              dset.attrs[tag + ' c=0.3'] = float(temp[1])
              dset.attrs[tag + ' c=0.4'] = float(temp[2])
              dset.attrs[tag + ' c=0.5'] = float(temp[3])
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Finally record combined results in separate group if present
      if start == 'low':
        this_str = this_ens + '/combo'
        this_grp = f.create_group(this_str)

        # Treat Wpoly_mod and Wpoly_r separately, to record c=0.5
        tag = 'Wpoly_mod c=0.5 '
        for line in open('results/Wpoly_mod.suscept-combo'):
          temp = line.split()
          this_grp.attrs[tag + temp[0]] = float(temp[1])
          this_grp.attrs[tag + temp[0] + '_err'] = float(temp[2])
          Nblocks = int(temp[-1])
        this_grp.attrs['Nblocks'] = Nblocks

        tag = 'Wpoly_r c=0.5 '
        for line in open('results/Wpoly.suscept-combo'):
          temp = line.split()
          this_grp.attrs[tag + temp[0]] = float(temp[1])
          this_grp.attrs[tag + temp[0] + '_err'] = float(temp[2])
          if not int(temp[-1]) == Nblocks:
            print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
            print("%s vs %d in Wpoly.suscept" % (temp[-1], Nblocks))
            sys.exit(1)

        # Other combined susceptibility results
        for obs in ['plaq', 'pbp', 'poly_mod', 'poly_r']:
          suscfile = 'results/' + obs + '.suscept-combo'
          for line in open(suscfile):
            temp = line.split()
            this_grp.attrs[obs + ' ' + temp[0]] = float(temp[1])
            this_grp.attrs[obs + ' ' + temp[0] + '_err'] = float(temp[2])
            if not int(temp[-1]) == Nblocks:
              print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
              print("%s vs %d in %s.suscept" % (temp[-1], Nblocks, obs))
              sys.exit(1)

        # Combined deconfinement fractions if present
        for theta in ['0.1', '0.15', '0.2', '0.25', '0.3']:
          fracfile = 'results/deconf_frac-combo-' + theta + '.dat'
          if not os.path.isfile(fracfile):
            continue
          for line in open(fracfile):
            temp = line.split()
            this_grp.attrs['deconf frac theta=' + theta] = float(temp[0])

          fracfile = 'results/Wdeconf_frac-combo-' + theta + '.dat'
          if not os.path.isfile(fracfile):
            continue
          tag = 'Wdeconf frac theta=' + theta
          for line in open(fracfile):
            if line.startswith('#'):
              continue
            temp = line.split()
            this_grp.attrs[tag + ' c=0.2'] = float(temp[0])
            this_grp.attrs[tag + ' c=0.3'] = float(temp[1])
            this_grp.attrs[tag + ' c=0.4'] = float(temp[2])
            this_grp.attrs[tag + ' c=0.5'] = float(temp[3])
# ------------------------------------------------------------------
