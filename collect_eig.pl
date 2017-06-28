#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw[min];
# ------------------------------------------------------------------
# This script accumulates and sorts eigenvalues from the given files
# Since the files have DDdag eigenvalues, take square root
# Cut off by minimum of last calculated eigenvalue
# We ASSUME that all files print at least Nvec EIGENVALUE

die "Usage: $0 <L> <Nt> <Nvec> <tag> <files>\n"
  if (@ARGV < 6);             # Require multiple data files

my $L = shift;                # To normalize by the volume
my $Nt = shift;
my $Nvec = shift;             # Number of eigenvalues to fit
my $tag = shift;              # Tag for output flies

my $junk, my $N, my $dat;   # To hold split data
my @eig = ();                 # 1d array for all eigenvalues
my @max = ();                 # Array for $Nvec - 1 eigenvalues
my $iter = 0;                 # To count number of files
my $minmax = 99;              # Minimum value of last eigenvalue
FILE: while (@ARGV > 0) {
  my $file = shift;
  open IN, "< $file" or die "Error opening $file ($!)\n";
  my @in = <IN>;
  close IN;

  $iter++;              # Must have before hitting "next FILE" below
  for my $line (@in) {
    if ($line =~ /^EIGENVALUE /) {
      ($junk, $N, $dat, $junk)  = split /\s+/, $line;
      if ($N == $Nvec - 1) {
        push(@max, sqrt($dat));
      }
      elsif ($N >= $Nvec) {
        next FILE;      # Don't bother pushing data we won't use later
      }
      push(@eig, sqrt($dat));
    }
  }
}
$minmax = min(@max);
# Normalize by volume and number of eigenvalue measurements
my $norm = $iter * $L**3 * $Nt;
#print "$iter * $L^3 * $Nt = $norm...\n";

# Print out sorted eigenvalues and sum, up to min of max
my @sorted = sort { $a <=> $b } @eig;
my $i, my $sum = 0;
open OUT, "> eig.$tag" or die "Error opening eig.$tag ($!)\n";
for ($i = 0; $i <= $#sorted; $i++) {
  $sum += 2 / $norm;           # Factor of two for complex conjugate pairs
  if ($sorted[$i] > $minmax) {
    last;
  }
  # Including the $norm on every line is redundant,
  # but may make things easier to awk or plot
  printf OUT "%.6g %.8g %.4g\n", $sorted[$i], $sum, $norm;
}
close OUT;
exit(0);
# ------------------------------------------------------------------
