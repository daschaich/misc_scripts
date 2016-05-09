#!/usr/bin/perl -w
use strict;
use warnings;
# ------------------------------------------------------------------
# This scripts parses lists of MILC output files, to figure out how
# many ensembles there are, and which files below to which ensembles

die "Usage: $0 <dir> <input files (full path)>\n"
  if (@ARGV < 1);

my $dir = shift;
my $write_path = "/nfs/beowulf03/beowulf02/anna2/KSNHYP/FA-ks481216/Run$dir";

while (@ARGV > 0) {
  my $file = shift;

  # Extract information from file name
  my $file_number = $file;
  $file_number =~ s/^.*\.(\d+)$/$1/;

  my $tag = $file;
  $tag =~ s/^.*\/out_(.*)\.\d+$/$1/;

  # Construct list of distinct ensembles
  my @tags;
  if (!(grep {$_ eq $tag} @tags)) {
    push(@tags, $tag);
    open LIST, ">> $write_path/ensembles.txt" or die "Error opening $write_path/ensembles.txt ($!)\n";
    print LIST "$tag\n";
    close LIST;
  }

  # Construct lists of files in each ensemble
  open TEMP, ">> $write_path/$tag.txt" or die "Error opening $write_path/$tag.txt ($!)\n";
  print TEMP "$file_number\n";
  close TEMP;
} # Loop over files
exit(0);
# ------------------------------------------------------------------
