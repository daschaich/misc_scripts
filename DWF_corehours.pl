#!/usr/bin/perl -w
use strict;
use warnings;
# ------------------------------------------------------------------
# This script cycles through the given output files,
# printing out and summing up the number of core-hours they used.
#
# Last modified 7 June 2015 by David Schaich
# ------------------------------------------------------------------



# ------------------------------------------------------------------
die "Usage: $0 <files> \n"
  if (@ARGV < 1);

my $cores;
my $seconds;
my $TEMP;
my $junk;

my $ch = 0;
my $iter = 0;
while (@ARGV > 0) {
  my $file = shift;
  $TEMP = `grep cpus $file`;
  ($junk, $junk, $junk, $junk, $junk, $cores, $junk) = split /\s+/, $TEMP;
  $TEMP = `grep "Inversion time = " $file`;
  ($junk, $junk, $junk, $seconds, $junk) = split /\s+/, $TEMP;
#  $cores = `grep cpus $file | awk '{print $6}'`;
#  $seconds = `grep "Time = " $file | awk '{print $3}'`;
  $TEMP = $cores * $seconds / 3600;
  printf "$file: %.0f\n", $TEMP;
  $ch += $TEMP;
  $iter++;
#  print "$cores * $seconds / 3600 --> $ch\n"
}
printf "Total: %.0f (%.0f files)\n", $ch, $iter;

exit(0);
# ------------------------------------------------------------------
