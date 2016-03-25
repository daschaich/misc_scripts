#!/usr/bin/perl -w
use strict;
use warnings;
# ------------------------------------------------------------------
# This script forces Wilson-flowed MCRG output files to print the
# correct flow time t corresponding to the MCRG blocking output
#
# Last modified 27 November 2012 by David Schaich
# ------------------------------------------------------------------



# ------------------------------------------------------------------
die "Usage: $0 <files>\n"
  if (@ARGV < 1);

my $t;
my $toCheck;
my $junk;

while (@ARGV > 0) {
  my $file = shift;

  open IN, "< $file" or die "Error opening $file ($!)\n";
  open OUT, "> temp" or die "Error opening temp ($!)\n";
  my @in = <IN>;
  close IN;

  for my $line (@in) {
    if ($line =~ /^WFLOW /) { # Increment t
      ($junk, $t, $junk, $junk, $junk, $junk, $junk, $junk) = split /\s+/, $line; # WFLOW  t  plaq  E  t^2*E  t*d(t^2*E)/dt  12t^2*(3-plaq)  top.charge
    }
    elsif ($line =~ /^LOOPS /) {
      ($junk, $toCheck, $junk, $junk, $junk, $junk, $junk) = split /\s+/, $line; # LOOPS  t  nloop  nrep  bl  alpha  value
      if ($toCheck != $t) {
        $line =~ s/LOOPS $toCheck/LOOPS $t/;
      }
    }
    elsif ($line =~ /^POLYA ORIG /) {
      ($junk, $junk, $toCheck, $junk, $junk, $junk, $junk) = split /\s+/, $line; # POLYA ORIG  t  plp.r  plp.i  xplp.r  xplp.i
      if ($toCheck != $t) {
        $line =~ s/POLYA ORIG $toCheck/POLYA ORIG $t/;
      }
    }
    elsif ($line =~ /^POLYA NHYP /) {
      ($junk, $junk, $toCheck, $junk, $junk, $junk, $junk, $junk, $junk) = split /\s+/, $line; # POLYA ORIG  t  bl  alpha  plp.r  plp.i  xplp.r  xplp.i
      if ($toCheck != $t) {
        $line =~ s/POLYA NHYP $toCheck/POLYA NHYP $t/;
      }
    }
    print OUT "$line";
  }
  close OUT;
  rename ("temp", $file);
}
exit(0);
# ------------------------------------------------------------------
