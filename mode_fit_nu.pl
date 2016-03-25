#!/usr/bin/perl -w
use strict;
use warnings;
# ------------------------------------------------------------------
# Run jackknifed linear fits of log(nu) vs. log(la_ave)
die "Usage: $0 <Nvec> <range> <shift> <tag> <files>\n"
  if (@ARGV < 6);         # Require multiple files

my $Nvec = shift;         # Number of eigenvalues to fit
my $range = int(shift);   # Number of modes per fit (integer)
my $shift = int(shift);   # Amount to shift fit range (integer)
my $name = shift;         # Tag base for output flies

my $tag = "${name}_D${range}_S${shift}_N$Nvec"; # Output tag
my $meas, my $N, my $i;                         # For loop indices
my $nu, my $dat, my $err, my $temp, my $junk;   # To hold split data

die "ERROR: NEED LARGER FIT RANGE THAN $range.\n"
  if ($range < 2);
die "ERROR: SHIFTING BY $shift WON'T GET YOU ANYWHERE.\n"
  if ($shift < 1);

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
my $minmax = 99;    # Minimum value of last eigenvalue
my @in = (), my @files = (), my @eig = (), my @all_eig = ();
while (@ARGV > 0) {
  my $file = shift;
  push @files, $file;
  open IN, "< $file" or die "Error opening $file ($!)\n";
  @in = <IN>;
  close IN;

  @eig = ();                # Reset array for each file
  for my $line (@in) {
    if ($line =~ /^EIGENVALUE /) {
      ($junk, $nu, $temp, $junk) = split /\s+/, $line;
      if ($nu >= $Nvec) {   # May cause problems with accumulator arrays
        next;               # And won't be used later
      }
      $dat = sqrt($temp);
      push @eig, $dat;
      $tot_la[$nu] += $dat;
      $tot_laSq[$nu] += $temp;

      # Find minimum of maximum lambda
      if ($nu == ($Nvec - 1) && $dat < $minmax) {
        $minmax = $dat;
      }
    }
  }
  push @all_eig, [@eig];    # All data from all files in 2d array
  $Nmeas++;
}

# See how many fit ranges we have to simplify future loops
my $Nfits = 0;
my $nu_start = 1;
for ($nu_start = 1;
     $nu_start + $range < $Nvec && $tot_la[$nu_start + $range - 1] / $Nmeas < $minmax;
     $nu_start += $shift) {
#  printf "Fit %d: [%d, %d]\n", $Nfits, $nu_start, $nu_start + $range - 1;
  $Nfits++;
}
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct jackknife samples through single elimination
my $lo = 0, my $hi = 0;                   # Fit range in log(nu)
my @fit, my $fit_pts, my $fitChi;         # More output from fitter
my @JKgamma = (), my @JKerr = (), my @JKchiSq = ();
my @gamma = (), my @all_err = (), my @chiSq = ();
for ($meas = 0; $meas < $Nmeas; $meas++) {
  for ($N = 0; $N < $Nvec; $N++) {
    $la_ave[$N] = $tot_la[$N] - $all_eig[$meas][$N];
    $la_err[$N] = $tot_laSq[$N] - $all_eig[$meas][$N]**2;
    $la_ave[$N] /= ($Nmeas - 1);
    $la_err[$N] /= ($Nmeas - 1);
    $la_err[$N] -= $la_ave[$N]**2;
    $la_err[$N] = sqrt($la_err[$N]);
#    printf "%.4g\n", $la_ave[$N]";
  }

  # Print file to be fit (error on log(la) is delta_lambda / lambda)
  open TOFIT, "> TOFIT_$tag" or die "Error opening TOFIT_$tag ($!)\n";
  for ($N = 0; $N < $Nvec; $N++) {
    printf TOFIT "%.6g %.6g %.6g\n",
    log($N + 1), log($la_ave[$N]), ($la_err[$N] / $la_ave[$N]);
  }
  close TOFIT;

  # Reset output arrays for each jackknife sample
  @JKgamma = ();
  @JKerr = ();

  # More robust loop over fit ranges in nu
  my $nu_start = 1;
  for ($i = 0; $i < $Nfits; $i++) {
    $lo = log($nu_start) - 0.00001;
    $hi = log($nu_start + $range - 1) + 0.00001;
#    printf "[%.4g, %.4g] --> [%.4g, %.4g]\n",
#           $nu_start, $nu_start + $range - 1, $lo, $hi;

    # Extract information from last three lines of fit results
    @fit = split /\n/, `fit2 $lo $hi < TOFIT_$tag | tail -n 3`;
    ($junk, $junk, $dat, $junk, $err) = split /\s+/, $fit[0];
#    ($junk, $junk, $junk, $junk, $junk, $fit_pts) = split /\s+/, $fit[1];
    ($junk, $fitChi) = split /\s+/, $fit[2];

    # CHECK that fit_pts equals $range -- comment out to increase performance
#    if ($fit_pts != $range) {
#      print "ERROR: fit_pts=$fit_pts but range=$range\n";
#      printf "[%d, %d] --> [%.4g, %.4g]\n",
#            $nu_start, $nu_start + $range - 1, $lo, $hi;
#      exit(0);
#    }

    push @JKgamma, $dat;
    push @JKerr, $err;
    push @JKchiSq, $fitChi;
    $nu_start += $shift;
  }

  # Store results from all jackknife samples in 2d arrays
  push @gamma, [@JKgamma];
  push @all_err, [@JKerr];
  push @chiSq, [@JKchiSq];
}
$junk = `rm -f TOFIT_$tag`;
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
open OUT, "> jack_lam_gamma_$tag" or die "Error opening jack_lam_gamma_$tag ($!)\n";
print OUT "# Jackknifed results from $Nmeas measurements";
print OUT " (listed at bottom)\n";
print OUT "# Delta_nu = $range\n";
print OUT "# lambda delta_lambda gamma delta_gamma chiSq/dof\n";

# Loop over fit ranges in nu in same way as above
my $ga_ave = 0, my $ga_Sqd = 0, my $ga_err = 0, my $chi_ave = 0;
$nu_start = 1;
for ($i = 0; $i < $Nfits; $i++) {
  # Reset averages using first jackknife sample
  $ga_ave = $gamma[0][$i];
  $ga_Sqd = $gamma[0][$i]**2;
  $ga_err = $all_err[0][$i];
  $chi_ave = $chiSq[0][$i];
  for ($N = 1; $N < $Nmeas; $N++) {
    $ga_ave += $gamma[$N][$i];
    $ga_Sqd += $gamma[$N][$i]**2;
    $ga_err += $all_err[$N][$i];
    $chi_ave += $chiSq[$N][$i];
  }
  $ga_ave /= $Nmeas;
  $ga_Sqd /= $Nmeas;
  $ga_err /= $Nmeas;
  $ga_err += sqrt($ga_Sqd - $ga_ave**2);    # Copied from bash script
  $chi_ave /= $Nmeas;

#  printf "[%.d, %.d]\n", $nu_start, $nu_start + $range - 1;
  $dat = ($la_ave[$nu_start + $range - 1] + $la_ave[$nu_start]) / 2;
  $err = ($la_ave[$nu_start + $range - 1] - $la_ave[$nu_start]) / 2;
  printf OUT "%.6g %.6g %.6g %.6g %.3g\n", $dat, $err, $ga_ave, $ga_err, $chi_ave;
  $nu_start += $shift;
}

$runtime += time();
print OUT "# Runtime: $runtime seconds\n";

# Print list of files
print OUT "# Files analyzed:\n# ";
print OUT join("\n# ", @files), "\n";
close OUT;
exit(0);
# ------------------------------------------------------------------
