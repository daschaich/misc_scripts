#!/usr/bin/perl -w
use strict;
use warnings;
# ------------------------------------------------------------------
# This script cycles through the given files,
# printing out and summing up the number of core-hours they used.

die "Usage: $0 <files> \n"
  if (@ARGV < 1);

my $cores;
my $seconds;
my $TEMP;
my $TEMP2;
my $junk;

my $ch = 0;
my $iter = 0;
while (@ARGV > 0) {
  my $file = shift;
  $TEMP = `grep Machine $file`;
  ($junk, $junk, $junk, $junk, $junk, $cores, $junk) = split /\s+/, $TEMP;
#  $TEMP = `egrep 'cpus|nodes' $file`;
#  ($junk, $junk, $junk, $junk, $junk, $cores, $junk) = split /\s+/, $TEMP;
  $TEMP = `grep "Time = " $file`;
  ($junk, $junk, $seconds, $junk) = split /\s+/, $TEMP;
#  $cores = `grep cpus $file | awk '{print $6}'`;
#  $seconds = `grep "Time = " $file | awk '{print $3}'`;
  $TEMP = $cores * $seconds;
  $TEMP2 = $cores * $seconds / 3600;
  printf "$file: %.2f / 3600 --> %.2f\n", $TEMP, $TEMP2;
  $ch += $TEMP2;
  $iter++;
#  print "$ch\n"
}
printf "Total: %.2f (%.0f files)\n", $ch, $iter;

exit(0);
# ------------------------------------------------------------------
