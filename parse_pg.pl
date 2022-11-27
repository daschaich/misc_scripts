#!/usr/bin/perl -w
use strict;
use warnings;
use Switch;
use Cwd 'abs_path';
# ------------------------------------------------------------------
# This script parses pure-gauge output files for a single ensemble,
# shuffling the extracted data into dedicated files for plotting.

if (! -d "Out") {
  die "ERROR: Out/ does not exist"
}

my $i;
my $junk;
open ERRFILE, "> ERRORS" or die "Error opening ERRORS ($!)\n";
open MISSINGFILES, "> MISSING" or die "Error opening MISSING ($!)\n";

# Set up header lines in new data/ files
open KEY, "> data/key.csv" or die "Error opening data/key.csv ($!)\n";
open STEPSIZE, "> data/stepsize.csv" or die "Error opening data/stepsize.csv ($!)\n";
open NSTEP, "> data/Nstep.csv" or die "Error opening data/Nstep.csv ($!)\n";
open TLENGTH, "> data/tlength.csv" or die "Error opening data/tlength.csv ($!)\n";
open TU, "> data/TU.csv" or die "Error opening data/TU.csv ($!)\n";
open DELTAS, "> data/deltaS.csv" or die "Error opening data/deltaS.csv ($!)\n";
open EXP_DS, "> data/exp_dS.csv" or die "Error opening data/exp_dS.csv ($!)\n";
open ABS_DS, "> data/abs_dS.csv" or die "Error opening data/abs_dS.csv ($!)\n";
open ACCP, "> data/accP.csv" or die "Error opening data/accP.csv ($!)\n";
open FORCE, "> data/force.csv" or die "Error opening data/force.csv ($!)\n";
open PLAQ, "> data/plaq.csv" or die "Error opening data/plaq.csv ($!)\n";
open POLY, "> data/poly.csv" or die "Error opening data/poly.csv ($!)\n";
open POLY_R, "> data/poly_r.csv" or die "Error opening data/poly_r.csv ($!)\n";
open POLY_MOD, "> data/poly_mod.csv" or die "Error opening data/poly_mod.csv ($!)\n";
open POLY_ARG, "> data/poly_arg.csv" or die "Error opening data/poly_arg.csv ($!)\n";
open XPOLY, "> data/xpoly.csv" or die "Error opening data/xpoly.csv ($!)\n";
open XPOLY_R, "> data/xpoly_r.csv" or die "Error opening data/xpoly_r.csv ($!)\n";
open XPOLY_MOD, "> data/xpoly_mod.csv" or die "Error opening data/xpoly_mod.csv ($!)\n";
open XPOLY_ARG, "> data/xpoly_arg.csv" or die "Error opening data/xpoly_arg.csv ($!)\n";
open PBP, "> data/pbp.csv" or die "Error opening data/pbp.csv ($!)\n";
open WALLTIME, "> data/walltime.csv" or die "Error opening data/walltime.csv ($!)\n";
open WALLTU, "> data/wallTU.csv" or die "Error opening data/wallTU.csv ($!)\n";
open PLAQB, "> data/plaqB.csv" or die "Error opening data/plaqB.csv ($!)\n";
open POLYB, "> data/polyB.csv" or die "Error opening data/polyB.csv ($!)\n";
open POLY_RB, "> data/poly_rB.csv" or die "Error opening data/poly_rB.csv ($!)\n";
open POLY_MODB, "> data/poly_modB.csv" or die "Error opening data/poly_modB.csv ($!)\n";
open POLY_ARGB, "> data/poly_argB.csv" or die "Error opening data/poly_argB.csv ($!)\n";
open XPOLYB, "> data/xpolyB.csv" or die "Error opening data/xpolyB.csv ($!)\n";
open XPOLY_RB, "> data/xpoly_rB.csv" or die "Error opening data/xpoly_rB.csv ($!)\n";
open XPOLY_MODB, "> data/xpoly_modB.csv" or die "Error opening data/xpoly_modB.csv ($!)\n";
open XPOLY_ARGB, "> data/xpoly_argB.csv" or die "Error opening data/xpoly_argB.csv ($!)\n";
open WFLOW, "> data/Wflow.csv" or die "Error opening data/Wflow.csv ($!)\n";
open TOPO, "> data/topo.csv" or die "Error opening data/topo.csv ($!)\n";
open WPOLY, "> data/Wpoly.csv" or die "Error opening data/Wpoly.csv ($!)\n";
open WPOLY_MOD, "> data/Wpoly_mod.csv" or die "Error opening data/Wpoly_mod.csv ($!)\n";
open WFLOW_SS, "> data/Wflow_ss.csv" or die "Error opening data/Wflow_ss.csv ($!)\n";
open WFLOW_ST, "> data/Wflow_st.csv" or die "Error opening data/Wflow_st.csv ($!)\n";
open WFLOW_ANISO, "> data/Wflow_aniso.csv" or die "Error opening data/Wflow_aniso.csv ($!)\n";

print KEY "t,file\n";
print STEPSIZE "t,eps\n";
print NSTEP "t,N\n";
print TLENGTH "t,L\n";
print TU "t,MDTU\n";
print DELTAS "t,deltaS\n";
print EXP_DS "t,e^(-dS)\n";
print ABS_DS "t,|deltaS|\n";
print ACCP "t,accP\n";
print FORCE "t,F0,F1,Fgauge\n";
print PLAQ "MDTU,plaq_ss,plaq_st\n";
print POLY "ReTr(L),ImTr(L)\n";
print POLY_R "MDTU,ReTr(L)\n";
print POLY_MOD "MDTU,|Tr(L)|\n";
print POLY_ARG "MDTU,arg(Tr(L))\n";
print XPOLY "ReTr(W_0),ImTr(W_0)\n";
print XPOLY_R "MDTU,ReTr(W_0)\n";
print XPOLY_MOD "MDTU,|Tr(W_0)|\n";
print XPOLY_ARG "MDTU,arg(Tr(W_0))\n";
#print PBP "MDTU,Re(pbp),Im(pbp)\n";
print PBP "MDTU,Re(pbp)\n";
print WALLTIME "t,walltime\n";
print WALLTU "t,cost\n";
print PLAQB "MDTU,bl0,bl1,bl2,bl3,bl4\n";
print POLYB "ReTr(L_b),bl0,bl1,bl2,bl3,bl4\n";
print POLY_RB "MDTU,bl0,bl1,bl2,bl3,bl4\n";
print POLY_MODB "MDTU,bl0,bl1,bl2,bl3,bl4\n";
print POLY_ARGB "MDTU,bl0,bl1,bl2,bl3,bl4\n";
print XPOLYB "Re(W_0),bl0,bl1,bl2,bl3,bl4\n";
print XPOLY_RB "MDTU,bl0,bl1,bl2,bl3,bl4\n";
print XPOLY_MODB "MDTU,bl0,bl1,bl2,bl3,bl4\n";
print XPOLY_ARGB "MDTU,bl0,bl1,bl2,bl3,bl4\n";
print WFLOW "MDTU,c=0.2,c=0.25,c=0.3,c=0.35\n";
print TOPO "MDTU,c=0.2,c=0.3,c=0.4,c=0.5\n";
print WPOLY "MDTU,c=0.2,c=0.3,c=0.4,c=0.5\n";
print WPOLY_MOD "MDTU,c=0.2,c=0.3,c=0.4,c=0.5\n";
print WFLOW_SS "MDTU,c=0.2,c=0.3,c=0.4,c=0.5\n";
print WFLOW_ST "MDTU,c=0.2,c=0.3,c=0.4,c=0.5\n";
print WFLOW_ANISO "MDTU,c=0.2,c=0.3,c=0.4,c=0.5\n";
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Declare observables to be printed
my $cpus;
my $stepsize;
my $Nstep;
my $tlength;
my $acc;
my $force;
my $dS;
my $exp_dS;
my $abs_dS;
my $iters;
my $plaq_ss;
my $plaq_st;
my $ploop_r;
my $ploop_i;
my $Wploop_r;
my $Wploop_i;
my $p_mod;
my $p_arg;
my $xloop_r;
my $xloop_i;
my $x_mod;
my $x_arg;
my $r_e;
my $r_o;
my $i_e;
my $i_o;
my $pbp;
my $ave_time;
my $TUtime;
my $bl;
my $alpha;

# Only calculate this once
my $Nc = -1.0;
my $path = abs_path("Out");
if ($path =~ /SU4/) {
  $Nc = 4.0;
}
elsif ($path =~ /SU6/) {
  $Nc = 6.0;
}
elsif ($path =~ /SU8/) {
  $Nc = 8.0;
}
else {
  die "ERROR: Nc unrecognized in $path"
}
my $gprop = 128.0 * 3.14159**2 / (3.0 * ($Nc * $Nc - 1.0));

# These are not currently used
my $startS;
my $endS;
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# These fellows help track status and perform checks
my $infile = -1;
my $check = -1;           # Check that file is present or completed
my $temp;                 # To hold topological charge input
my $temp_i;               # To hold Wilson line input
my $walltime = -1;        # Check that file completed successfully
my $level = 0;            # Track inner vs. outer level for step size
my $load;                 # MDTU of loaded configuration
my $cfg;                  # MDTU of saved configuration
my $oldcfg = 0;           # Check if any files are missing
my $stamp = "start";
my $Wflow_stamp = "start";
my $oldstamp = "start";   # Check that correct configuration was used

# Running sums for the ensemble as a whole
my $traj = 0;
my $endtraj = 0;          # Reset from traj, so not really running sum
my $MDTU = 0;

# Running sums for each trajectory
my $pbp_r = 0;
my $pbp_i = 0;
my $pbp_iter = 0;
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Cycle through files, first "out" and then "Wflow"
# list.txt is a list of $load-$save
open FILES, "< list.txt" or die "Error opening list.txt ($!)\n";
my @files = <FILES>;
close FILES;
FILE: for my $file (@files) {
  chomp ($file);          # Get rid of linebreaks in filenames
  ($load, $cfg) = split /-/, $file;
  # Initialize running sums and set dummy walltime
  # If the latter isn't overwritten, then the run died
  # or its output file is corrupted
  $pbp_r = 0;
  $pbp_i = 0;
  $pbp_iter = 0;
  $walltime = -1;
  $stamp = "start";

  # Open file
  # If not found, move on to next file instead of killing whole program,
  # but print error message so I know there is a problem
  $infile = "Out/out.$file";
  $check = open IN, "< $infile";
  if (!$check) {
    print STDERR "Problem opening Out/out.$file: $!\n";
    print ERRFILE "Problem opening Out/out.$file: $!\n";
    next FILE;
  }
#  print STDOUT "$infile\n"; # Monitor running status
  my @in = <IN>;
  close IN;

  # If not starting from file 1 in this ensemble,
  # or if we seem to have skipped a file,
  # guess approximate starting trajectory
  my $traj_per_file = -1;
  my $L = -1;
  my $Nt = -1;
  my $vol = -1;
  LINE: for my $line (@in) {
    if ($line =~ /^PLACEHOLDER/) {
      # Placeholder file -- error has been addressed as well as possible,
      # but don't print nonsense wall clock time
      $walltime = -2;
    }
    elsif ($line =~ /^trajecs /) {
      ($junk, $traj_per_file) = split /\s+/, $line;
      $endtraj = $traj + $traj_per_file;
      last LINE;   # Don't go through whole file yet
    }
    elsif ($line =~ /^nx /) {
      ($junk, $L) = split /\s+/, $line;
    }
    elsif ($line =~ /^nt /) {
      ($junk, $Nt) = split /\s+/, $line;
      $vol = $L**3 * $Nt;
      # Now set $L to minimum of nx and nt for Wilson flow
      if ($Nt < $L) {
        $L = $Nt;
      }
    }
  }
  if ($traj_per_file == -1) {
    print STDERR "$infile: number of trajectories is never defined\n";
    print ERRFILE "$infile: number of trajectories is never defined\n";
    next FILE;     # Skip this file, which appears malformed
  }
  elsif ($L == -1 || $vol == -1) {
    print STDERR "$infile: Lattice volume is never defined\n";
    print ERRFILE "$infile: Lattice volume is never defined\n";
    next FILE;     # Skip this file, which appears malformed
  }

  if ($traj == 0 && $load > 0 || $load != $oldcfg) {
    print STDERR "$infile: guessing approximate starting trajectory\n";
    $traj = $load;
    $endtraj = $traj + $traj_per_file;
  }
  # --------------------------------------------------------------



  # --------------------------------------------------------------
  # Cycle through lines in the "out" file
  for my $line (@in) {
    # Check for unitarity problems
    if ($line =~ /^Unitarity problem /) {
      print STDERR "$infile: unitarity problem reported\n";
      print ERRFILE "$infile: unitarity problem reported\n";
      next FILE;     # Skip this file, which appears malformed
    }
    # Check for premature termination (e.g., layout problem)
    if ($line =~ /^termination/) {
      print STDERR "$infile: premature termination reported\n";
      print ERRFILE "$infile: premature termination reported\n";
      next FILE;     # Skip this file, which appears malformed
    }
  }

  # At this point, we should be able to begin
  $oldcfg = $cfg;
  for my $line (@in) {
    # Extract constant run parameters
    if ($line =~ /^hmc_steps /) {
      ($junk, $Nstep) = split /\s+/, $line;
      print NSTEP "$endtraj,$Nstep\n";
    }

    elsif ($line =~ /^traj_length/) {
      ($junk, $tlength) = split /\s+/, $line;
      print TLENGTH "$endtraj,$tlength\n";
      $stepsize = $tlength / $Nstep;
      print STEPSIZE "$endtraj,$stepsize\n";
    }

    elsif ($line =~ /^Machine = /) {
      ($junk, $junk, $junk, $junk, $junk, $cpus, $junk) = split /\s+/, $line;
    }
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Check that the file loaded the appropriate configuration
    elsif ($line =~ /^Restoring binary SciDAC file /) { # Loading configuration, but don't have Time stamp
      $stamp = "loaded";
    }
    elsif ($line =~ /^Time stamp /) {
      chomp ($line);              # Remove linebreak from end of line
      if ($stamp eq "start") {    # Loading configuration
        $stamp = join ' ', split ' ', $line;    # Replace multiple spaces with single spaces
        if ($stamp ne $oldstamp && $oldstamp ne "start") {
          print STDERR "$infile: loaded $stamp doesn't match saved $oldstamp\n";
          print ERRFILE "$infile: loaded $stamp doesn't match saved $oldstamp\n";
        }
      }
      else {                # Saving configuration
        $oldstamp = join ' ', split ' ', $line;  # Prepare for next file or MCRG file
        $stamp = "start";
      }
    }
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Now extract evolution observables and physical observables
    # Acceptance comes before measurements
    elsif ($line =~ / delta S = /) {
      $traj++;            # Remains across files in this ensemble!
      $MDTU += $tlength;  # Remains across files in this ensemble!
      # Round off digits in MDTU
      $junk = $MDTU;
      $MDTU = sprintf("%.3f", $junk);
      print TU "$traj,$MDTU\n";
      print KEY "$MDTU,$file\n";

      ($acc, $junk, $junk, $junk, $dS, $junk, $junk, $junk, $startS, $junk, $junk, $junk, $endS) = split /\s+/, $line;
      print DELTAS "$traj,$dS\n";
      $exp_dS = exp(-$dS);
      print EXP_DS "$traj,$exp_dS\n";

      # For RMS, don't have an easy way to average over many measurements
      # Instead just print out absolute value and consider its running average
      $abs_dS = abs($dS);
      print ABS_DS "$traj,$abs_dS\n";

      # Will be smeared out by running averages
      if ($acc =~ /ACCEPT/) {
        print ACCP "$traj,1\n";
      }
      else {
        print ACCP "$traj,0\n";
      }
    }

    # Force --- interested in max rather than average
    elsif ($line =~ /MONITOR_FORCE /) {
      ($junk, $junk, $force) = split /\s+/, $line;
      print FORCE "$traj,$force\n";
    }
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Gauge measurements come next: plaquette and Polyakov loop
    # Normalize Wilson lines to 1
    elsif ($line =~ /^GMES/) {
      ($junk, $temp, $temp_i, $plaq_ss, $plaq_st) = split /\s+/, $line;
      print PLAQ "$MDTU,$plaq_ss,$plaq_st\n";

      $ploop_r = $temp / $Nc;
      $ploop_i = $temp_i / $Nc;
      print POLY "$ploop_r,$ploop_i\n";
      print POLY_R "$MDTU,$ploop_r\n";
      $p_mod = sqrt($ploop_r**2 + $ploop_i**2);
      print POLY_MOD "$MDTU,$p_mod\n";
      $p_arg = atan2($ploop_i, $ploop_r);
      print POLY_ARG "$MDTU,$p_arg\n";
    }

    # Wilson line in the x direction
    elsif ($line =~ /^POLYA/) {
      ($junk, $junk, $junk, $temp, $temp_i) = split /\s+/, $line;
      $xloop_r = $temp / $Nc;
      $xloop_i = $temp_i / $Nc;
      print XPOLY "$xloop_r,$xloop_i\n";
      print XPOLY_R "$MDTU,$xloop_r\n";
      $x_mod = sqrt($xloop_r**2 + $xloop_i**2);
      print XPOLY_MOD "$MDTU,$x_mod\n";
      $x_arg = atan2($xloop_i, $xloop_r);
      print XPOLY_ARG "$MDTU,$x_arg\n";
    }
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Finally, the chiral condensate is only measured every once in a while
    # but is measured stochastically, so we need to average
    elsif ($line =~ /^PBP:/) {
#      ($junk, $junk, $junk, $r_e, $r_o, $i_e, $i_o, $junk, $i, $junk, $Nsrc, $junk, $iters) = split /\s+/, $line;
      ($junk, $junk, $junk, $r_e, $r_o, $i_e, $i_o, $junk) = split /\s+/, $line;
      $pbp_r += $r_e + $r_o;
      $pbp_i += $i_e + $i_o;

      # !!!Hard-code Nsrc=5 to avoid fatal errors from malformed input
      $pbp_iter++;
      my $Nsrc = 5;
      if ($pbp_iter == $Nsrc) {  # Average, print and reset
        $pbp_r /= (2 * $Nsrc);
        $pbp_i /= (2 * $Nsrc);
#        print PBP "$MDTU,$pbp_r,$pbp_i\n"
        print PBP "$MDTU,$pbp_r\n";
        $pbp_r = 0;
        $pbp_i = 0;
        $pbp_iter = 0;
      }
    }

    # Store total walltime to average at the end
    elsif ($line =~ /^Time = /) {
      ($junk, $junk, $walltime, $junk) = split /\s+/, $line;
    }
  } # Done cycling through output file
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Check to see if run seems to have finished properly
  if ($walltime == -1) {
    print STDERR "$infile: final timing never reported\n";
    print ERRFILE "$infile: final timing never reported\n";
  }
  elsif ($walltime == -2) {
    # Placeholder file -- error has been addressed as well as possible,
    # but don't print nonsense wall clock time
  }
  else {    # We are good to go
    # Average walltime over all trajectories
    $ave_time = $walltime / $traj_per_file;
    print WALLTIME "$traj,$ave_time\n";
    $TUtime = $ave_time * $cpus / (60 * abs($tlength));
    print WALLTU "$traj,$TUtime\n";
  }
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
WFLOW:
  # Check out Wilson flow, printing gSq for c=0.2, 0.25, 0.3 and 0.35
  # These files now contain MCRG-blocked results: only consider t=0
  # Also ignore finite-volume correction in definition of gSq
  $check = -1;
  $infile = "Out/Wflow.$cfg";
  $check = open WFLOW_IN, "< $infile";
  # File may not be present if configuration was not saved
  if (!$check) {
    print MISSINGFILES "$infile\n";
    next FILE;
  }
  my @Wflow_in = <WFLOW_IN>;
  close WFLOW_IN;

  # We have a file, so let's cycle over its lines
  $check = -1;              # Check whether file completed successfully
  my $Wflow_t = -1;         # Flow time
  my $E_ss = -1;            # Space--space t^2*E
  my $E_st = -1;            # Space--time t^2*E
  my $tSqE = -1;            # To check, related to coupling
  my $c = -1;               # sqrt(8t) / L
  my $cOld = 1;             # To see where we pass thresholds
  my $cp = -1;              # sqrt(8t) / L
  my @gSq = ("null", "null", "null", "null");
  my @topo = ("null", "null", "null", "null");  # Topological charge
  my @poly = ("null", "null", "null", "null");  # Wilson-flowed Polyakov loop real part
  my @mod = ("null", "null", "null", "null");   # Wilson-flowed Polyakov loop modulus
  my @Wflow_ss = ("null", "null", "null", "null");  # t^2 * E_ss
  my @Wflow_st = ("null", "null", "null", "null");  # t^2 * E_st
  my @aniso = ("null", "null", "null", "null");  # Anisotropy E_ss / E_st

  # RG-blocked observables
  my $loop;                   # Which loop is this (plaquette is zero)
  my $plaq_blocked = -1;      # Check to see whether these are here to print
  my $poly_blocked = -1;
  my @plaq = ("null", "null", "null", "null", "null");   # Initialize output
  my @ploop_r = ("null", "null", "null", "null", "null");
  my @ploop_i = ("null", "null", "null", "null", "null");
  my @xloop_r = ("null", "null", "null", "null", "null");
  my @xloop_i = ("null", "null", "null", "null", "null");
  my @p_mod = ("null", "null", "null", "null", "null");
  my @p_arg = ("null", "null", "null", "null", "null");
  my @x_mod = ("null", "null", "null", "null", "null");
  my @x_arg = ("null", "null", "null", "null", "null");

  # Go
  $Wflow_stamp = "start";   # Check loading/saving configuration
  for my $line (@Wflow_in) {
    if ($line =~ /^Time stamp /) {
      if ($Wflow_stamp eq "start") {  # Loading the configuration
        chomp ($line);                # Remove linebreak from end of line
        $Wflow_stamp = join ' ', split ' ', $line;    # Replace multiple spaces with single spaces
        if ($Wflow_stamp ne $oldstamp) {
          print STDERR "$infile: loaded $Wflow_stamp doesn't match saved $oldstamp\n";
          print ERRFILE "$infile: loaded $Wflow_stamp doesn't match saved $oldstamp\n";
        }
      }
      else {                          # Saving the configuration
        $Wflow_stamp = "start";
      }
    }
    elsif ($line =~ /^WFLOW /) { # Save tSqE in case we want to look at that
      ($junk, $Wflow_t, $junk, $junk, $tSqE, $junk, $junk, $temp, $E_ss, $E_st) = split /\s+/, $line;
      $c = sqrt(8*$Wflow_t) / $L;
      if ($cOld < 0.2 && $c >= 0.2) {
        $gSq[0] = $gprop * $tSqE;
        $topo[0] = $temp;
        $Wflow_ss[0] = $E_ss;
        $Wflow_st[0] = $E_st;
        $aniso[0] = $E_ss / $E_st;
      }
      elsif ($cOld < 0.25 && $c >= 0.25) {
        $gSq[1] = $gprop * $tSqE;
      }
      elsif ($cOld < 0.3 && $c >= 0.3) {
        $gSq[2] = $gprop * $tSqE;
        $topo[1] = $temp;
        $Wflow_ss[1] = $E_ss;
        $Wflow_st[1] = $E_st;
        $aniso[1] = $E_ss / $E_st;
      }
      elsif ($cOld < 0.35 && $c >= 0.35) {
        $gSq[3] = $gprop * $tSqE;
      }
      elsif ($cOld < 0.4 && $c >= 0.4) {
        $topo[2] = $temp;
        $Wflow_ss[2] = $E_ss;
        $Wflow_st[2] = $E_st;
        $aniso[2] = $E_ss / $E_st;
      }
      elsif ($cOld < 0.5 && $c >= 0.5) {
        $topo[3] = $temp;
        $Wflow_ss[3] = $E_ss;
        $Wflow_st[3] = $E_st;
        $aniso[3] = $E_ss / $E_st;
      }
      $cOld = $c;
    }
    # Wilson-flowed Polyakov loop
    # Normalized to Nc rather than 1
    elsif ($line =~ /^POLYA ORIG /) {
      ($junk, $junk, $Wflow_t, $Wploop_r, $Wploop_i, $xloop_r, $xloop_i) = split /\s+/, $line;
      $cp = sqrt(8.0 * $Wflow_t) / $L;
      if ($cp == 0) {
        $ploop_r[0] = $Wploop_r;
        $ploop_i[0] = $Wploop_i;
        $xloop_r[0] = $xloop_r;
        $xloop_i[0] = $xloop_i;
      }
      if (0.19 < $cp && $cp < 0.21) {
        $poly[0] = $Wploop_r;
        $mod[0] = sqrt($Wploop_r**2 + $Wploop_i**2);
      }
      elsif (0.29 < $cp && $cp < 0.31) {
        $poly[1] = $Wploop_r;
        $mod[1] = sqrt($Wploop_r**2 + $Wploop_i**2);
      }
      elsif (0.39 < $cp && $cp < 0.41) {
        $poly[2] = $Wploop_r;
        $mod[2] = sqrt($Wploop_r**2 + $Wploop_i**2);
      }
      elsif (0.49 < $cp && $cp < 0.51) {
        $poly[3] = $Wploop_r;
        $mod[3] = sqrt($Wploop_r**2 + $Wploop_i**2);
      }
    }

    # RG-blocked observables (with no Wilson flow)
    elsif ($line =~ /^LOOPS 0 /) {
      ($junk, $junk, $loop, $junk, $bl, $alpha, $junk) = split /\s+/, $line;
      if ($loop == 0 && $alpha == 0) {
        ($junk, $junk, $loop, $junk, $bl, $alpha, $plaq[0]) = split /\s+/, $line;
      }
      elsif ($loop == 0 && $alpha == 0.6 && 1 <= $bl && $bl <= 4) {
        $plaq_blocked = 1;
        ($junk, $junk, $loop, $junk, $bl, $alpha, $plaq[$bl]) = split /\s+/, $line;
      }
    }
    elsif ($line =~ /^POLYA NHYP 0 /) {
      ($junk, $junk, $junk, $bl, $alpha, $junk) = split /\s+/, $line;
      if ($alpha == 0.6) {
        ($junk, $junk, $junk, $bl, $alpha, $ploop_r[$bl], $ploop_i[$bl], $xloop_r[$bl], $xloop_i[$bl]) = split /\s+/, $line;
        if ($bl > $poly_blocked) {
          $poly_blocked = $bl;
        }
      }
    }
    elsif ($line =~ /RUNNING COMPLETED/) {
      $check = 1;
    }
  } # Done with Wilson flow file
  if ($check == -1) {
    print STDERR "$infile did not complete\n";
    print ERRFILE "$infile did not complete\n";
  }
  print WFLOW "$MDTU,$gSq[0],$gSq[1],$gSq[2],$gSq[3]\n";
  print TOPO "$MDTU,null,null,null,$topo[3]\n";
  print WPOLY "$MDTU,$poly[0],$poly[1],$poly[2],$poly[3]\n";
  print WPOLY_MOD "$MDTU,$mod[0],$mod[1],$mod[2],$mod[3]\n";
  print WFLOW_SS "$MDTU,$Wflow_ss[0],$Wflow_ss[1],$Wflow_ss[2],$Wflow_ss[3],\n";
  print WFLOW_ST "$MDTU,$Wflow_st[0],$Wflow_st[1],$Wflow_st[2],$Wflow_st[3],\n";
  print WFLOW_ANISO "$MDTU,$aniso[0],$aniso[1],$aniso[2],$aniso[3],\n";

  # Lots of RG-blocked stuff to print
  # Print plaquettes
  if ($plaq_blocked > 0) {
    print PLAQB "$MDTU,$plaq[0],$plaq[1],$plaq[2],$plaq[3],$plaq[4]\n";
  }

  # Print Wilson lines -- assume that all elements filled properly
  if ($poly_blocked > 0) {
    for ($i = 0; $i <= $poly_blocked; $i++) {
      $p_mod[$i] = sqrt($ploop_r[$i]**2 + $ploop_i[$i]**2);
      $p_arg[$i] = atan2($ploop_i[$i], $ploop_r[$i]);
      $x_mod[$i] = sqrt($xloop_r[$i]**2 + $xloop_i[$i]**2);
      $x_arg[$i] = atan2($xloop_i[$i], $xloop_r[$i]);
    }

    print POLYB "$ploop_r[0],$ploop_i[0],null,null,null,null\n";
    print POLYB "$ploop_r[1],null,$ploop_i[1],null,null,null\n";
    print XPOLYB "$xloop_r[0],$xloop_i[0],null,null,null,null\n";
    print XPOLYB "$xloop_r[1],null,$xloop_i[1],null,null,null\n";
    if ($poly_blocked >= 2) {
      print POLYB "$ploop_r[2],null,null,$ploop_i[2],null,null\n";
      print XPOLYB "$xloop_r[2],null,null,$xloop_i[2],null,null\n";
    }
    if ($poly_blocked >= 3) {
      print POLYB "$ploop_r[3],null,null,null,$ploop_i[3],null\n";
      print XPOLYB "$xloop_r[3],null,null,null,$xloop_i[3],null\n";
    }
    if ($poly_blocked == 4) {
      print POLYB "$ploop_r[4],null,null,null,null,$ploop_i[4]\n";
      print XPOLYB "$xloop_r[4],null,null,null,null,$xloop_i[4]\n";
    }
    print POLY_RB "$MDTU,$ploop_r[0],$ploop_r[1],$ploop_r[2],$ploop_r[3],$ploop_r[4]\n";
    print POLY_MODB "$MDTU,$p_mod[0],$p_mod[1],$p_mod[2],$p_mod[3],$p_mod[4]\n";
    print POLY_ARGB "$MDTU,$p_arg[0],$p_arg[1],$p_arg[2],$p_arg[3],$p_arg[4]\n";
    print XPOLY_RB "$MDTU,$xloop_r[0],$xloop_r[1],$xloop_r[2],$xloop_r[3],$xloop_r[4]\n";
    print XPOLY_MODB "$MDTU,$x_mod[0],$x_mod[1],$x_mod[2],$x_mod[3],$x_mod[4]\n";
    print XPOLY_ARGB "$MDTU,$x_arg[0],$x_arg[1],$x_arg[2],$x_arg[3],$x_arg[4]\n";
  } # Done printing MCRG data
} # Done cycling through files
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Clean up and close down
close ERRFILE;
close KEY;
close STEPSIZE;
close NSTEP;
close TLENGTH;
close TU;
close DELTAS;
close EXP_DS;
close ABS_DS;
close ACCP;
close FORCE;
close PLAQ;
close POLY;
close POLY_R;
close POLY_MOD;
close POLY_ARG;
close XPOLY;
close XPOLY_R;
close XPOLY_MOD;
close XPOLY_ARG;
close PBP;
close WALLTIME;
close WALLTU;
close PLAQB;
close POLYB;
close POLY_RB;
close POLY_MODB;
close POLY_ARGB;
close XPOLYB;
close XPOLY_RB;
close XPOLY_MODB;
close XPOLY_ARGB;
close WFLOW;
close TOPO;
close WPOLY;
close WPOLY_MOD;
close WFLOW_SS;
close WFLOW_ST;
close WFLOW_ANISO;
exit(0);
# ------------------------------------------------------------------
