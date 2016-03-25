#!/usr/bin/perl -w
use strict;
use warnings;
use Switch;
# ------------------------------------------------------------------
# This script parses MILC output files for a single ensemble,
# shuffling the extracted data into dedicated files for plotting.

die "Usage: $0 <path>\n"
  if (@ARGV != 1);

my $path = shift;
my $i;
my $junk;
open ERRFILE, "> ERRORS" or die "Error opening ERRORS ($!)\n";
open MISSINGFILES, "> MISSING" or die "Error opening MISSING ($!)\n";

# Set up header lines in new data/ files
open KEY, "> $path/data/key.csv" or die "Error opening $path/data/key.csv ($!)\n";
open STEPSIZE, "> $path/data/stepsize.csv" or die "Error opening $path/data/stepsize.csv ($!)\n";
open NSTEP, "> $path/data/Nstep.csv" or die "Error opening $path/data/Nstep.csv ($!)\n";
open MH, "> $path/data/MH.csv" or die "Error opening $path/data/MH.csv ($!)\n";
open TLENGTH, "> $path/data/tlength.csv" or die "Error opening $path/data/tlength.csv ($!)\n";
open TU, "> $path/data/TU.csv" or die "Error opening $path/data/TU.csv ($!)\n";
open DELTAS, "> $path/data/deltaS.csv" or die "Error opening $path/data/deltaS.csv ($!)\n";
open EXP_DS, "> $path/data/exp_dS.csv" or die "Error opening $path/data/exp_dS.csv ($!)\n";
open ABS_DS, "> $path/data/abs_dS.csv" or die "Error opening $path/data/abs_dS.csv ($!)\n";
open ACCP, "> $path/data/accP.csv" or die "Error opening $path/data/accP.csv ($!)\n";
open FORCE, "> $path/data/force.csv" or die "Error opening $path/data/force.csv ($!)\n";
open PLAQ, "> $path/data/plaq.csv" or die "Error opening $path/data/plaq.csv ($!)\n";
open POLY, "> $path/data/poly.csv" or die "Error opening $path/data/poly.csv ($!)\n";
open POLY_R, "> $path/data/poly_r.csv" or die "Error opening $path/data/poly_r.csv ($!)\n";
open POLY_MOD, "> $path/data/poly_mod.csv" or die "Error opening $path/data/poly_mod.csv ($!)\n";
open POLY_ARG, "> $path/data/poly_arg.csv" or die "Error opening $path/data/poly_arg.csv ($!)\n";
open XPOLY, "> $path/data/xpoly.csv" or die "Error opening $path/data/xpoly.csv ($!)\n";
open XPOLY_R, "> $path/data/xpoly_r.csv" or die "Error opening $path/data/xpoly_r.csv ($!)\n";
open XPOLY_MOD, "> $path/data/xpoly_mod.csv" or die "Error opening $path/data/xpoly_mod.csv ($!)\n";
open XPOLY_ARG, "> $path/data/xpoly_arg.csv" or die "Error opening $path/data/xpoly_arg.csv ($!)\n";
open CG_ITERS, "> $path/data/cg_iters.csv" or die "Error opening $path/data/cg_iters.csv ($!)\n";
open PBP, "> $path/data/pbp.csv" or die "Error opening $path/data/pbp.csv ($!)\n";
open WALLTIME, "> $path/data/walltime.csv" or die "Error opening $path/data/walltime.csv ($!)\n";
open WALLTU, "> $path/data/wallTU.csv" or die "Error opening $path/data/wallTU.csv ($!)\n";
open PLAQB, "> $path/data/plaqB.csv" or die "Error opening $path/data/plaqB.csv ($!)\n";
open POLYB, "> $path/data/polyB.csv" or die "Error opening $path/data/polyB.csv ($!)\n";
open POLY_RB, "> $path/data/poly_rB.csv" or die "Error opening $path/data/poly_rB.csv ($!)\n";
open POLY_MODB, "> $path/data/poly_modB.csv" or die "Error opening $path/data/poly_modB.csv ($!)\n";
open POLY_ARGB, "> $path/data/poly_argB.csv" or die "Error opening $path/data/poly_argB.csv ($!)\n";
open XPOLYB, "> $path/data/xpolyB.csv" or die "Error opening $path/data/xpolyB.csv ($!)\n";
open XPOLY_RB, "> $path/data/xpoly_rB.csv" or die "Error opening $path/data/xpoly_rB.csv ($!)\n";
open XPOLY_MODB, "> $path/data/xpoly_modB.csv" or die "Error opening $path/data/xpoly_modB.csv ($!)\n";
open XPOLY_ARGB, "> $path/data/xpoly_argB.csv" or die "Error opening $path/data/xpoly_argB.csv ($!)\n";
open PLAQ_DIFF, "> $path/data/plaq_diff.csv" or die "Error opening $path/data/plaq_diff.csv ($!)\n";
open LINK_DIFF, "> $path/data/link_diff.csv" or die "Error opening $path/data/link_diff.csv ($!)\n";
open EIG, "> $path/data/eig.csv" or die "Error opening $path/data/eig.csv ($!)\n";
open WFLOW, "> $path/data/Wflow.csv" or die "Error opening $path/data/Wflow.csv ($!)\n";
open TOPO, "> $path/data/topo.csv" or die "Error opening $path/data/topo.csv ($!)\n";

print KEY "t,file\n";
print STEPSIZE "t,eps0,eps1,eps_gauge\n";
print NSTEP "t,N0,N1,Ngauge\n";
print MH "t,MH\n";
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
print CG_ITERS "t,cg_iters\n";
print PBP "MDTU,Re(pbp),Im(pbp)\n";
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
print PLAQ_DIFF "MDTU,t,x,y,z,Norm\n";
print LINK_DIFF "MDTU,t,x,y,z,Norm\n";
print EIG "MDTU,1,2,3,4,5,6,7,8,9,10,11,12\n";
print WFLOW "MDTU,c=0.2,c=0.25,c=0.3,c=0.35\n";
print TOPO "MDTU,c=0.2,c=0.3,c=0.4,c=0.5\n";
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Declare observables to be printed
my $cpus;
my $stepsize;
my $stepsize1;
my $stepsize_gauge;
my $Nstep;
my $Nstep1;
my $Nstep_gauge;
my $MH;
my $tlength;
my $acc;
my $force0;
my $force1;
my $force_gauge;
my $dS;
my $exp_dS;
my $abs_dS;
my $iters;
my $plaq_ss;
my $plaq_st;
my $ploop_r;
my $ploop_i;
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
my $gprop = 128 * 3.14159**2 / (3 * 8);

# These fellows are not currently used
my $startS;
my $endS;
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# These fellows help track status and perform checks
my $infile = -1;
my $check = -1;           # Check that file is present or completed
my $temp;                 # To hold topological charge input
my $walltime = -1;        # Check that file completed successfully
my $level = 0;            # Track inner vs. outer level for step size
my $force_counter = 0;    # Track how many forces we have seen
my $load;                 # MDTU of loaded configuration
my $cfg;                  # MDTU of saved configuration
my $oldcfg = 0;           # Check if any files are missing
my $stamp = "start";
my $mcrg_stamp = "";
my $order_stamp = "";
my $eig_stamp = "";
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
# Cycle through files, first "out_" and then "mcrg_", "Sout_" and "eig_"
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

  # Open file, which may be in one of two directories
  # If not found, move on to next file instead of killing whole program,
  # but print error message so I know there is a problem
  $infile = "$path/Out/out.$file";
  $check = open IN, "< $infile";
  if (!$check) {
    print STDERR "Problem opening $path/Out/out.$file: $!\n";
    print ERRFILE "Problem opening $path/Out/out.$file: $!\n";
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
    if ($line =~/^PLACEHOLDER/) {
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
  # Cycle through lines in the "out_" file
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
    if ($line =~ /^traj_length/) {
      ($junk, $tlength) = split /\s+/, $line;
      print TLENGTH "$endtraj,$tlength\n";
    }

    elsif ($line =~ /^nstep /) {
      if ($level == 0) {
        ($junk, $Nstep) = split /\s+/, $line;
        $stepsize = $tlength / $Nstep;
        $level = 1;
      }
      else {
        ($junk, $Nstep1) = split /\s+/, $line;
        $stepsize1 = $stepsize / (2 * $Nstep1);
        $Nstep1 *= 2 * $Nstep;
        $level = 0;
      }
    }
    elsif ($line =~ /^nstep_gauge /) {
      ($junk, $Nstep_gauge) = split /\s+/, $line;
      $stepsize_gauge = $stepsize1 / (2 * $Nstep_gauge);
      $Nstep_gauge *= 2 * $Nstep1;
      print NSTEP "$endtraj,$Nstep,$Nstep1,$Nstep_gauge\n";
      print STEPSIZE "$endtraj,$stepsize,$stepsize1,$stepsize_gauge\n";
    }

    elsif ($line =~ /^Machine = /) {
      ($junk, $junk, $junk, $junk, $junk, $cpus, $junk) = split /\s+/, $line;
    }

    elsif ($line =~ /^Hasenbusch_mass /) {
      ($junk, $MH) = split /\s+/, $line;
      print MH "$endtraj,$MH\n";
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

   # Forces -- order rearranged during code development
   elsif ($line =~ /MONITOR_FORCE_FERMION0 /) {
     ($junk, $force0) = split /\s+/, $line;
     $force_counter++;
     if ($force_counter == 3) {
       print FORCE "$traj,$force0,$force1,$force_gauge\n";
       $force_counter = 0;
     }
   }
   elsif ($line =~ /MONITOR_FORCE_FERMION1 /) {
     ($junk, $force1) = split /\s+/, $line;
     $force_counter++;
     if ($force_counter == 3) {
       print FORCE "$traj,$force0,$force1,$force_gauge\n";
       $force_counter = 0;
     }
   }
   elsif ($line =~ /MONITOR_FORCE_GAUGE /) {
     ($junk, $force_gauge) = split /\s+/, $line;
     $force_counter++;
     if ($force_counter == 3) {
       print FORCE "$traj,$force0,$force1,$force_gauge\n";
       $force_counter = 0;
     }
   }
   # ------------------------------------------------------------

   # ------------------------------------------------------------
   # Gauge measurements come next:
   # plaquette and Polyakov loop, and CG iterations for some reason
   elsif ($line =~ /^GMES/) {
     ($junk, $ploop_r, $ploop_i, $iters, $plaq_ss, $plaq_st) = split /\s+/, $line;
     print PLAQ "$MDTU,$plaq_ss,$plaq_st\n";
     print CG_ITERS "$traj,$iters\n";

     print POLY "$ploop_r,$ploop_i\n";
     print POLY_R "$MDTU,$ploop_r\n";
     $p_mod = sqrt($ploop_r**2 + $ploop_i**2);
     print POLY_MOD "$MDTU,$p_mod\n";
     $p_arg = atan2($ploop_i, $ploop_r);
     print POLY_ARG "$MDTU,$p_arg\n";
   }

    # Wilson line in the x direction
    elsif ($line =~ /^POLYA/) {
      ($junk, $ploop_r, $ploop_i, $xloop_r, $xloop_i) = split /\s+/, $line;
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
MCRG:
  # Now deal with the corresponding MCRG file, if it is there
  # Need to extract the configuration number from $file="$load-$save"
  $check = -1;
  $infile = "$path/Out/mcrg.$cfg";
  $check = open MCRG, "< $infile";
  # File may not be present if configuration was not saved
  if (!$check) {
#    print MISSINGFILES "$infile\n";
    goto ORDER;
  }
  my @mcrg = <MCRG>;
  close MCRG;

  # We have a file, so let's cycle over its lines
  # !!! Smearing parameter alpha=0.6 hard-coded!!!
  $check = -1;                # Check whether file completed successfully
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
  for my $line (@mcrg) {
    # Check that the file loaded the same configuration as above
    if ($line =~ /^Time stamp /) {
      chomp ($line);              # Remove linebreak from end of line
      $mcrg_stamp = join ' ', split ' ', $line;    # Replace multiple spaces with single spaces
      if ($mcrg_stamp ne $oldstamp) {
        print STDERR "$infile: loaded $mcrg_stamp doesn't match saved $oldstamp\n";
        print ERRFILE "$infile: loaded $mcrg_stamp doesn't match saved $oldstamp\n";
      }
    }
    elsif ($line =~ /^LOOPS /) {
      ($junk, $loop, $junk, $bl, $alpha, $junk) = split /\s+/, $line;
      if ($loop == 0 && $alpha == 0) {
        ($junk, $loop, $junk, $bl, $alpha, $plaq[0]) = split /\s+/, $line;
      }
      elsif ($loop == 0 && $alpha == 0.6 && 1 <= $bl && $bl <= 4) {
        $plaq_blocked = 1;
        ($junk, $loop, $junk, $bl, $alpha, $plaq[$bl]) = split /\s+/, $line;
      }
    }

    elsif ($line =~ /^POLYA NHYP /) {
      ($junk, $junk, $bl, $alpha, $junk) = split /\s+/, $line;
      if ($alpha == 0.6) {
        ($junk, $junk, $bl, $alpha, $ploop_r[$bl], $ploop_i[$bl], $xloop_r[$bl], $xloop_i[$bl]) = split /\s+/, $line;
        if ($bl > $poly_blocked) {
          $poly_blocked = $bl;
        }
      }
    }
    elsif ($line =~ /RUNNING COMPLETED/) {
      $check = 1;
    }
  } # Done cycling through lines in MCRG file
  if ($check == -1) {
    print STDERR "$infile did not complete\n";
    print ERRFILE "$infile did not complete\n";
  }

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
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
ORDER:
  # Now deal with the corresponding Sout file, if it is there
  $check = -1;
  $infile = "$path/Out/Sout.$cfg";
  $check = open ORDER, "< $infile";
  # File may not be present if configuration was not saved
  if (!$check) {
#    print MISSINGFILES "$infile\n";
    goto EIGEN;
  }
  my @order = <ORDER>;
  close ORDER;

  # We have a file, so let's cycle over its lines
  $check = -1;               # Check whether file completed successfully
  my @plaq_diff = ("null", "null", "null", "null", "null");
  my @link_diff = ("null", "null", "null", "null", "null");
  for my $line (@order) {
    if ($line =~ /^Time stamp /) {
      chomp ($line);              # Remove linebreak from end of line
      $order_stamp = join ' ', split ' ', $line;    # Replace multiple spaces with single spaces
      if ($order_stamp ne $oldstamp) {
        print STDERR "$infile: loaded $order_stamp doesn't match saved $oldstamp\n";
        print ERRFILE "$infile: loaded $order_stamp doesn't match saved $oldstamp\n";
      }
    }
    elsif ($line =~ /^StaggPlaq t /) {
      ($junk, $junk, $r_e, $r_o) = split /\s+/, $line;
      $plaq_diff[0] = ($r_e - $r_o);
    }
    elsif ($line =~ /^StaggPlaq x /) {
      ($junk, $junk, $r_e, $r_o) = split /\s+/, $line;
      $plaq_diff[1] = ($r_e - $r_o);
    }
    elsif ($line =~ /^StaggPlaq y /) {
      ($junk, $junk, $r_e, $r_o) = split /\s+/, $line;
      $plaq_diff[2] = ($r_e - $r_o);
    }
    elsif ($line =~ /^StaggPlaq z /) {
      ($junk, $junk, $r_e, $r_o) = split /\s+/, $line;
      $plaq_diff[3] = ($r_e - $r_o);
    }
    elsif ($line =~ /^pbpt: /) {
      ($junk, $junk, $junk, $r_e, $r_o, $i_e, $i_o, $junk) = split /\s+/, $line;
      $link_diff[0] = ($r_e - $r_o);
    }
    elsif ($line =~ /^pbpx: /) {
      ($junk, $junk, $junk, $r_e, $r_o, $i_e, $i_o, $junk) = split /\s+/, $line;
      $link_diff[1] = ($r_e - $r_o);
    }
    elsif ($line =~ /^pbpy: /) {
      ($junk, $junk, $junk, $r_e, $r_o, $i_e, $i_o, $junk) = split /\s+/, $line;
      $link_diff[2] = ($r_e - $r_o);
    }
    elsif ($line =~ /^pbpz: /) {
      ($junk, $junk, $junk, $r_e, $r_o, $i_e, $i_o, $junk) = split /\s+/, $line;
      $link_diff[3] = ($r_e - $r_o);
    }
    elsif ($line =~ /RUNNING COMPLETED/) {
      $check = 1;
    }
  } # Done with order parameter file
  if ($plaq_diff[0] ne "null") {   # Calculate norms
    $plaq_diff[4] = 0;
    for ($i = 0; $i < 4; $i++) {
      $plaq_diff[4] += $plaq_diff[$i]**2;
    }
    $plaq_diff[4] = sqrt($plaq_diff[4]);
  }
  # Some older measurements didn't include the link difference
  if ($link_diff[0] ne "null") {
    $link_diff[4] = 0;
    for ($i = 0; $i < 4; $i++) {
      $link_diff[4] += $link_diff[$i]**2;
    }
    $link_diff[4] = sqrt($link_diff[4]);
  }
  if ($check == -1) {
    print STDERR "$infile did not complete\n";
    print ERRFILE "$infile did not complete\n";
  }
  print PLAQ_DIFF "$MDTU,$plaq_diff[0],$plaq_diff[1],$plaq_diff[2],$plaq_diff[3],$plaq_diff[4]\n";
  print LINK_DIFF "$MDTU,$link_diff[0],$link_diff[1],$link_diff[2],$link_diff[3],$link_diff[4]\n";
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
EIGEN:
  # Now deal with the corresponding eigenvalues file, if it is there
  $check = -1;
  $infile = "$path/Out/eig.$cfg";
  $check = open EIG_IN, "< $infile";
  # File may not be present if configuration was not saved
  if (!$check) {
#    print MISSINGFILES "$infile\n";
    goto WFLOW;
  }
  my @eig_in = <EIG_IN>;
  close EIG_IN;

  # We have a file, so let's cycle over its lines
  # Can't get all the eigenvalues, so let's check out the lowest six
  $check = -1;               # Check whether file completed successfully
  my @eig = ("null", "null", "null", "null", "null", "null");
  for my $line (@eig_in) {
    if ($line =~ /^Time stamp /) {
      chomp ($line);              # Remove linebreak from end of line
      $eig_stamp = join ' ', split ' ', $line;    # Replace multiple spaces with single spaces
      if ($eig_stamp ne $oldstamp) {
        print STDERR "$infile: loaded $eig_stamp doesn't match saved $oldstamp\n";
        print ERRFILE "$infile: loaded $eig_stamp doesn't match saved $oldstamp\n";
      }
    }
    elsif ($line =~ /^EIGENVALUE 0 /) {
      ($junk, $junk, $eig[0], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /^EIGENVALUE 1 /) {
      ($junk, $junk, $eig[1], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /^EIGENVALUE 2 /) {
      ($junk, $junk, $eig[2], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /^EIGENVALUE 3 /) {
      ($junk, $junk, $eig[3], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /^EIGENVALUE 4 /) {
      ($junk, $junk, $eig[4], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /^EIGENVALUE 5 /) {
      ($junk, $junk, $eig[5], $junk) = split /\s+/, $line;
    }
    elsif ($line =~ /WARNING/) {
      print STDERR "$infile saturated eigenvalue iterations\n";
      print ERRFILE "$infile saturated eigenvalue iterations\n";
    }
    elsif ($line =~ /RUNNING COMPLETED/) {
      $check = 1;
    }
  } # Done with eigenvalues file
  if ($check == -1) {
    print STDERR "$infile did not complete\n";
    print ERRFILE "$infile did not complete\n";
  }
  print EIG "$MDTU,$eig[0],$eig[1],$eig[2],$eig[3],$eig[4],$eig[5]\n";
  # ------------------------------------------------------------------



  # ----------------------------------------------------------------
WFLOW:
  # Check out Wilson flow, printing gSq for c=0.2, 0.25, 0.3 and 0.35
  # These files should not contain MCRG blocking; ignore it for now
  # Also ignore finite-volume correction in definition of gSq
  $check = -1;
  $infile = "$path/Out/Wflow.$cfg";
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
  my $norm = 1;             # Annoying normalization issue on different sites...
  my $Wflow_t = -1;         # Flow time
  my $tSqE = -1;            # To check, related to coupling
  my $c = -1;               # sqrt(8t) / L
  my $cOld = 1;             # To see where we pass thresholds
  my @gSq = ("null", "null", "null", "null");
  my @topo = ("null", "null", "null", "null");  # Topological charge
  $Wflow_stamp = "start";   # Check loading/saving configuration
  for my $line (@Wflow_in) {
    if ($line =~ /fnal.gov/) {      # Annoying normalization issue on different sites...
      $norm = $vol * 0.02533029591058444286**2;  # 1/4pi^2
    }
    elsif ($line =~ /^Time stamp /) {
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
      ($junk, $Wflow_t, $junk, $junk, $tSqE, $junk, $junk, $temp) = split /\s+/, $line;
      $c = sqrt(8*$Wflow_t) / $L;
      if ($cOld < 0.2 && $c >= 0.2) {
        $gSq[0] = $gprop * $tSqE;
        $topo[0] = $temp * $norm;
      }
      elsif ($cOld < 0.25 && $c >= 0.25) {
        $gSq[1] = $gprop * $tSqE;
      }
      elsif ($cOld < 0.3 && $c >= 0.3) {
        $gSq[2] = $gprop * $tSqE;
        $topo[1] = $temp * $norm;
      }
      elsif ($cOld < 0.35 && $c >= 0.35) {
        $gSq[3] = $gprop * $tSqE;
      }
      elsif ($cOld < 0.4 && $c >= 0.4) {
        $topo[2] = $temp * $norm;
      }
      elsif ($cOld < 0.5 && $c >= 0.5) {
        $topo[3] = $temp * $norm;
      }
      $cOld = $c;
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
} # Done cycling through files
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Clean up and close down
close ERRFILE;
close MISSINGFILES;
close KEY;
close STEPSIZE;
close NSTEP;
close MH;
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
close CG_ITERS;
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
close PLAQ_DIFF;
close LINK_DIFF;
close EIG;
close WFLOW;
close TOPO;
exit(0);
# ------------------------------------------------------------------
