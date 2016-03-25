#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw[min];
# ------------------------------------------------------------------
# Run jackknifed linear fits of log(nu) vs. log(la_ave)
die "Usage: $0 <Nvec> <min points> <name> <files>
       Shift hard-coded to range / 5.
       Initial fit range=0.0125 adjusted to accommodate min points.\n"
  if (@ARGV < 6);         # Require multiple data files

my $Nvec = shift;         # Number of eigenvalues to fit
my $FITMIN = shift;       # Require that at least this many points be fit
my $name = shift;         # Tag for output flies
my $range = 0.0125;       # Initial fit range, to be expanded if necessary
my $shift = $range / 5;   # Amount to shift fit range

my $meas, my $N, my $fit;                       # For loop indices
my $nu, my $dat, my $err, my $temp, my $junk;   # To hold split data

my $runtime = -time();
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First define and initialize accumulator arrays
my @tot_la, my @tot_laSq, my @la_ave, my @la_err;
for ($N = 0; $N < $Nvec; $N++) {
  push @tot_la, 0;
  push @tot_laSq, 0;
  push @la_ave, 0;
  push @la_err, 0;
}

# Now push all data into arrays, also accumulating totals
my $Nmeas = 0;      # Count number of files (jackknife samples)
my $min = 0;        # $range less than minimum value of $la_ave[$FITMIN - 1]
my $minmax = 99;    # Minimum value of last eigenvalue
my @files = (), my @eig, my @all_eig = ();
while (@ARGV > 0) {
  my $file = shift;
  push @files, $file;
  open IN, "< $file" or die "Error opening $file ($!)\n";
  my @in = <IN>;
  close IN;

  @eig = ();        # Reset array for each file
  for my $line (@in) {
    if ($line =~ /^EIGENVALUE /) {
      ($junk, $nu, $temp, $junk) = split /\s+/, $line;
      if ($nu >= $Nvec) {
        next;     # Don't bother pushing data we won't use later
      }
      $dat = sqrt($temp);
      push @eig, $dat;
      $tot_la[$nu] += $dat;
      $tot_laSq[$nu] += $temp;
    }
  }
  push @all_eig, [@eig];   # All data from all files in 2d array
  $Nmeas++;
}

# Find min and minmax
# Make min the smaller of the first average
# or $range less than the $FITMIN - 1 average
my @toMin = (), my @toMinMax = ();
for ($meas = 0; $meas < $Nmeas; $meas++) {
  $temp = $tot_la[0] - $all_eig[$meas][0];
  $temp /= ($Nmeas - 1);
  push @toMin, $temp;
  $temp = $tot_la[$FITMIN - 1] - $all_eig[$meas][$FITMIN - 1];
  $temp /= ($Nmeas - 1);
  push @toMin, ($temp - $range);
  push @toMinMax, $all_eig[$meas][$Nvec - 1];
}
$min = min(@toMin);
$minmax = min(@toMinMax);
#printf "%.4g, %.4g\n", $min, $minmax;
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct jackknife samples through single elimination
my $la_start = $min, my $la_end = $min + $range;  # Fit ranges in lambda
my $lo = 0, my $hi = 0;                           # Fit ranges in log(lambda)
my @fit = (), my $fit_pts, my $fitChi;            # More output from fitter
RESTART:
my $tag = "${name}_F${FITMIN}_D${range}_N$Nvec"; # Output tag
my @JKgamma = (), my @JKNpts = (), my @JKchiSq = ();
my @gamma = (), my @Npts = (), my @chiSq = ();
for ($meas = 0; $meas < $Nmeas; $meas++) {
  for ($N = 0; $N < $Nvec; $N++) {
    $la_ave[$N] = $tot_la[$N] - $all_eig[$meas][$N];
    $la_err[$N] = $tot_laSq[$N] - $all_eig[$meas][$N]**2;
    $la_ave[$N] /= ($Nmeas - 1);
    $la_err[$N] /= ($Nmeas - 1);
    $la_err[$N] -= $la_ave[$N]**2;
    $la_err[$N] = sqrt($la_err[$N]);
#    printf "%.4g\n", $la_ave[$N]";                        # CHECK
  }

  # Print file to be fit (error on log(la) is delta_lambda / lambda
  open TOFIT, "> TOFIT_$tag" or die "Error opening TOFIT_$tag ($!)\n";
  for ($N = 0; $N < $Nvec; $N++) {
    printf TOFIT "%.6g %.6g %.6g\n",
    log($N + 1), log($la_ave[$N]), ($la_err[$N] / $la_ave[$N]);
  }
  close TOFIT;

  # Reset output arrays for each jackknife sample
  @JKgamma = ();

  # Go through fit ranges in lambda
  $la_start = $min;
  $la_end = $min + $range;
  while ($la_end < $minmax) {
#    printf "[%.4g, %.4g] --> ", $la_start, $la_end;       # CHECK

    # Determine corresponding fit range in log(nu)
    $lo = -1;
    $hi = -1;
    for ($N = 0; $N < $Nvec; $N++) {  # May be more elegant way to do this
      if ($lo == -1 && $la_start <= $la_ave[$N]) {
        $lo = log($N + 1) - 0.00001;      # First data point after $la_start
      }
      if ($hi == -1 && $la_end < $la_ave[$N] && $N > 0) {
        $hi = log($N) + 0.00001;          # Last data point before $la_end
      }
    }
#    printf "[%.4g, %.4g] --> %.4g\n", $lo, $hi, $hi-$lo;  # CHECK
    # It seems this can happen even after shifting $min up to $la_ave[0]...
    # Different jackknife samples can have different numbers of skipped fits
    # Don't try to deal with it, just declare fit range too small and die
    if ($hi - $lo < 0.00003) {
      print "WARNING: NULL FIT.  ";
      printf "INCREASING RANGE %.4g --> ", $range;
      $range += $shift;
      printf "%.4g\n", $range;
      $junk = `rm -f TOFIT_$tag`;
      goto RESTART;
    }

    # Save fit results to file so we can get more than one line
    @fit = split /\n/, `fit2 $lo $hi < TOFIT_$tag | tail -n 3`;
    ($junk, $junk, $dat, $junk, $err) = split /\s+/, $fit[0];
    ($junk, $junk, $junk, $junk, $junk, $fit_pts) = split /\s+/, $fit[1];
    chomp($fit_pts);    # Strip newline character from end
    ($junk, $fitChi) = split /\s+/, $fit[2];

    # Since we still have problems after dealing with the above,
    # Require that our fit has at least $FITMIN data points
#    printf "%.4g %.4g %.4g %d %.4g\n", $la_start, $dat, $err, $fit_pts, $fitChi;
    if ($fit_pts < $FITMIN) {
      print "WARNING: $fit_pts < $FITMIN POINTS IN FIT.  ";
      printf "INCREASING RANGE %.4g --> ", $range;
      $range += $shift;
      printf "%.4g\n", $range;
      $junk = `rm -f TOFIT_$tag`;
      goto RESTART;
    }
    push @JKgamma, (4 * $dat - 1);
    push @JKNpts, $fit_pts;
    push @JKchiSq, $fitChi;

    # Shift to next fit range
    $la_start += $shift;
    $la_end = $la_start + $range;
  }

  # Store results from all jackknife samples in 2d arrays
  push @gamma, [@JKgamma];
  push @Npts, [@JKNpts];
  push @chiSq, [@JKchiSq];
}
$junk = `rm -f TOFIT_$tag`;
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
print "Writing output file jack_lam_gamma_$tag\n";
open OUT, "> jack_lam_gamma_$tag" or die "Error opening jack_lam_gamma_$tag ($!)\n";
print OUT "# Jackknifed results from $Nmeas measurements";
print OUT " (listed at bottom)\n";
print OUT "# Delta_lambda = $range\n";
print OUT "# lambda gamma delta_gamma Npts chiSq/dof\n";

# This should reproduce the correct number of fit ranges
my $Nfits = int(($minmax - $range - $min) / $shift);
$la_start = $min;
$la_end = $min + $range;
my $ga_ave = 0, my $ga_Sqd = 0, my $ga_err = 0;
my $chi_ave = 0; my $N_ave = 0;
for ($fit = 0; $fit < $Nfits; $fit++) {
  $ga_ave = $gamma[0][$fit];
  $ga_Sqd = $gamma[0][$fit]**2;
  $chi_ave = $chiSq[0][$fit];
  $N_ave = $Npts[0][$fit];
  for ($N = 1; $N < $meas; $N++) {
    $ga_ave += $gamma[$N][$fit];
    $ga_Sqd += $gamma[$N][$fit]**2;
    $chi_ave += $chiSq[$N][$fit];
    $N_ave += $Npts[$N][$fit];
  }
  $ga_ave /= $Nmeas;
  $ga_Sqd /= $Nmeas;
  $ga_err /= $Nmeas;
  $ga_err = sqrt(($Nmeas - 1) * ($ga_Sqd - $ga_ave**2));
  $chi_ave /= $Nmeas;
  $N_ave /= $Nmeas;

  $dat = ($la_start + $la_end) / 2;
  printf OUT "%.6g %.6g %.6g %d %.3g\n", $dat, $ga_ave, $ga_err, $N_ave, $chi_ave;
  $la_start += $shift;
  $la_end = $la_start + $range;
}

$runtime += time();
print OUT "# Runtime: $runtime seconds\n";

# Print list of files
print OUT "# Files analyzed:\n# ";
print OUT join("\n# ", @files), "\n";
close OUT;
exit(0);
# ------------------------------------------------------------------
