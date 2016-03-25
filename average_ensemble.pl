#!/usr/bin/perl -w
use POSIX;    # For floor
use strict;
use warnings;
# ------------------------------------------------------------------
# This script parses my dygraphs data files to construct averages
# and standard deviations given a thermalization cut.

# The data files must be ordered!
# The data files must be labelled "high", "low" or "mix"!

# TODO
# At least for eigenvalues, average into blocks of size $skip
# rather than simply selecting one measurement every $skip MDTU
# (Should also save resulting blocked data for jackknife analysis)

# !!! Eigenvalue analysis currently broken by name format change !!!
# ------------------------------------------------------------------



# ------------------------------------------------------------------
die "Usage: $0 <dir> <tag> <L> <Nt> <cut> <skip>\n"
  if (@ARGV != 6);

my $dir = shift;
my $tag = shift;
my $L = shift;
my $Nt = shift;
my $cut_in = shift;
my $skip = shift;
my $MDTU;
my $junk;
my $last_file;

print "$0 $dir $tag $L $Nt $cut_in $skip\n";

# Translate MDTU cut into corresponding output file
my $file_cut;
open KEY, "< $dir/data/key.$tag.csv" or die "Error opening $dir/data/key.$tag.csv ($!)\n";
my @key = <KEY>;
close KEY;
for my $line (@key) {
  ($MDTU, $file_cut)  = split /,/, $line;
  if ($MDTU eq "t") { next; }
  elsif ($MDTU <= $cut_in) { next; }
  else { last; }  # $file_cut now contains run corresponding to MDTU cut
}

# Figure out how many files there are in total
# Exit immediately if placeholder thermalization cut is after last MDTU
my $line = $key[$#key];
($MDTU, $last_file)  = split /,/, $line;
if ($MDTU < $cut_in) {
  exit(0);
}

# Figure out maximum number of RG blockings
if ($L < $Nt) {
  $junk = $L;
}
else {
  $junk = $Nt;
}

my $blMax = 0;
while ($junk % 2 == 0 && $junk > 2) {
  $junk /= 2;
  $blMax++;
}

# Now we are ready to go
# Will be overwritten -- need to reset for each file
my $cut = $cut_in;
my $next = $cut + $skip;

# Extract beta and mass from tag ~ volume_beta_ratio_mass
my $beta;
my $m;
if ($tag =~ /^.*_(\d+\.\d+)_.*\.\d+_\d+\.\d+$/) {
  $beta = $tag;
  $beta =~ s/^.*_(\d+\.\d+)_.*\.\d+_\d+\.\d+$/$1/;
  $m = $tag;
  $m =~ s/^.*_\d+\.\d+_.*\.\d+_(\d+\.\d+)$/$1/;
}
elsif ($tag =~ /^.*_(\d+\.\d+)_.*\.\d+_\d+$/) {
  $beta = $tag;
  $beta =~ s/^.*_(\d+\.\d+)_.*\.\d+_\d+$/$1/;
  $m = $tag;
  $m =~ s/^.*_\d+\.\d+_.*\.\d+_(\d+)$/$1/;
}
else {
  print STDERR "Problem: can't extract beta and mass from $tag\n";
  exit(1);
}

# Extract starting configuration from tag
my $init;
if ($tag =~ /high/) {
  $init = "hot";
}
elsif ($tag =~ /low/) {
  $init = "cold";
}
elsif ($tag =~ /mix/) {
  $init = "mix";
}
else {
  die;
}
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Plaquette -- two data (to be averaged) per line
my $count = 0;
my $ss;
my $st;
my $dat;
my $ave = 0;
my $err = 0;
open IN, "< $dir/data/plaq.$tag.csv" or die "Error opening $dir/data/plaq.$tag.csv ($!)\n";
my @in = <IN>;
close IN;
for my $line (@in) {
  ($MDTU, $ss, $st)  = split /,/, $line;
  if ($MDTU eq "MDTU") { next; }
  elsif ($MDTU < $cut) { next; }
  elsif ($cut <= $MDTU && $MDTU < $next) {
    $count++;
    $dat = ($ss + $st) / 2;
    $ave += $dat;
    $err += $dat**2;

    # Increment by $skip
    $cut = $next;
    $next = $cut + $skip;
    next;   # Make sure next elsif isn't caught, too
  }
  # This next case should happen at most once per ensemble,
  # if the initial data are missing (as for blocked Polyakov loops)
  elsif ($MDTU >= $next) {
    while ($MDTU >= $next) {
      $cut = $next;
      $next = $cut + $skip;
    }
    $count++;
    $dat = ($ss + $st) / 2;
    $ave += $dat;
    $err += $dat**2;

    # Increment by $skip
    $cut = $next;
    $next = $cut + $skip;
  }
}

# Average and print plaquette if there is data
# Store even a single measurement
if ($count >= 1) {
  open PLAQ, ">> $dir/results/plaq_m$m.$init" or die "Error opening $dir/results/plaq_m$m.$init ($!)\n";

  if ($count > 1) {
    $ave /= $count;
    $err /= $count;
    $err -= $ave**2;
    $err = sqrt($err / ($count - 1));    # Standard deviation
    printf PLAQ "%.4g %.4g %.4g %d\n", $beta, $ave, $err, $count;
  }
  elsif ($count == 1) {
    printf PLAQ "%.4g %.4g 0 1\n", $beta, $ave;
  }
  close PLAQ;
}
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Chiral condensate, Polyakov loop and spatial Wilson line
# Only one datum per line
my @files = ("pbp", "poly_r", "poly_mod", "xpoly_r", "xpoly_mod");
foreach my $file (@files) {
  $count = 0;
  $cut = $cut_in;
  $next = $cut + $skip;
  $ave = 0;
  $err = 0;
  open IN, "< $dir/data/$file.$tag.csv" or die "Error opening $dir/data/$file.$tag.csv ($!)\n";
  @in = <IN>;
  close IN;

  for my $line (@in) {
    ($MDTU, $dat)  = split /,/, $line;
    if ($MDTU eq "MDTU") { next; }
    elsif ($MDTU < $cut) { next; }
    elsif ($cut <= $MDTU && $MDTU < $next) {
      $count++;
      $ave += $dat;
      $err += $dat**2;

      # Increment by $skip
      $cut = $next;
      $next = $cut + $skip;
      next;   # Make sure next elsif isn't caught, too
    }
    elsif ($MDTU >= $next) {
      while ($MDTU >= $next) {
        $cut = $next;
        $next = $cut + $skip;
      }
      $count++;
      $ave += $dat;
      $err += $dat**2;

      # Increment by $skip
      $cut = $next;
      $next = $cut + $skip;
    }
  }

  # Average and print (only if there is data!)
  # Store even a single measurement
  if ($count >= 1) {
    open OUT, ">> $dir/results/${file}_m$m.$init" or die "Error opening $dir/results/${file}_m$m.$init ($!)\n";

    if ($count > 1) {
      $ave /= $count;
      $err /= $count;
      $err -= $ave**2;
      $err = sqrt($err / ($count - 1));    # Standard deviation
      printf OUT "%.4g %.4g %.4g %d\n", $beta, $ave, $err, $count;
    }
    elsif ($count == 1) {
      printf OUT "%.4g %.4g 0 1\n", $beta, $ave;
    }
    close OUT;
  }
}
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Blocked data files -- data per line depend on number of blockings
@files = ("plaqB", "poly_rB", "poly_modB", "xpoly_rB", "xpoly_modB");
my $i;
foreach my $file (@files) {
  $count = 0;
  $cut = $cut_in;
  $next = $cut + $skip;
  my @dat;
  my @ave;
  my @err;
  open IN, "< $dir/data/$file.$tag.csv" or die "Error opening $dir/data/$file.$tag.csv ($!)\n";
  @in = <IN>;
  close IN;

  # Initialize arrays
  for ($i = 0; $i <= $blMax; $i++) {
    $dat[$i] = 0;
    $ave[$i] = 0;
    $err[$i] = 0;
  }

  for my $line (@in) {
    ($MDTU, @dat)  = split /,/, $line;
    if ($MDTU eq "MDTU") { next; }
    elsif ($MDTU < $cut) { next; }
    elsif ($cut <= $MDTU && $MDTU < $next) {
      $count++;
      for ($i = 0; $i <= $blMax; $i++) {
        $ave[$i] += $dat[$i];
        $err[$i] += $dat[$i]**2;
      }

      # Increment by $skip
      $cut = $next;
      $next = $cut + $skip;
      next;   # Make sure next elsif isn't caught, too
    }
    elsif ($MDTU >= $next) {
      while ($MDTU >= $next) {
        $cut = $next;
        $next = $cut + $skip;
      }
      $count++;
      for ($i = 0; $i <= $blMax; $i++) {
        $ave[$i] += $dat[$i];
        $err[$i] += $dat[$i]**2;
      }

      # Increment by $skip
      $cut = $next;
      $next = $cut + $skip;
    }
  }

  # Average and print (only if there is data!)
  # Store even a single measurement
  if ($count >= 1) {
    open OUT, ">> $dir/results/${file}_m$m.$init" or die "Error opening $dir/results/${file}_m$m.$init ($!)\n";
    printf OUT "%.4g", $beta;
    if ($count > 1) {
      for ($i = 0; $i <= $blMax; $i++) {
        $ave[$i] /= $count;
        $err[$i] /= $count;
        $err[$i] -= $ave[$i]**2;
        $err[$i] = sqrt($err[$i] / ($count - 1));    # Standard deviation
        printf OUT " %.4g %.4g", $ave[$i], $err[$i];
      }
    }
    elsif ($count == 1) {
      for ($i = 0; $i <= $blMax; $i++) {
        printf OUT " %.4g 0", $ave[$i];
      }
    }
    printf OUT " %d\n", $count;
    close OUT;
  }
} # Done cycling over blocked data files
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Order parameters -- only grab euclidean norm on each line
# Normalize link diff by (half)-volume
@files = ("plaq_diff", "link_diff");
foreach my $file (@files) {
  $count = 0;
  $cut = $cut_in;
  $next = $cut + $skip;
  $ave = 0;
  $err = 0;
  open IN, "< $dir/data/$file.$tag.csv" or die "Error opening $dir/data/$file.$tag.csv ($!)\n";
  @in = <IN>;
  close IN;

  for my $line (@in) {
    ($MDTU, $junk, $junk, $junk, $junk, $dat)  = split /,/, $line;
    if ($MDTU eq "MDTU") { next; }
    elsif ($MDTU < $cut) { next; }
    elsif ($cut <= $MDTU && $MDTU < $next) {
      $count++;
      if ($file eq "link_diff") {
        $dat = $dat * 2 / ($L**3 * $Nt);    # Normalization missing from link diff
      }
      $ave += $dat;
      $err += $dat**2;

      # Increment by $skip
      $cut = $next;
      $next = $cut + $skip;
      next;   # Make sure next elsif isn't caught, too
    }
    elsif ($MDTU >= $next) {
      while ($MDTU >= $next) {
        $cut = $next;
        $next = $cut + $skip;
      }
      $count++;
      if ($file eq "link_diff") {
        $dat = $dat * 2 / ($L**3 * $Nt);    # Normalization missing from link diff
      }
      $ave += $dat;
      $err += $dat**2;

      # Increment by $skip
      $cut = $next;
      $next = $cut + $skip;
    }
  }

  # Average and print (only if there is data!)
  # Store even a single measurement
  if ($count >= 1) {
    open OUT, ">> $dir/results/${file}_m$m.$init" or die "Error opening $dir/results/${file}_m$m.$init ($!)\n";

    if ($count > 1) {
      $ave /= $count;
      $err /= $count;
      $err -= $ave**2;
      $err = sqrt($err / ($count - 1));    # Standard deviation
      printf OUT "%.4g %.4g %.4g %d\n", $beta, $ave, $err, $count;
    }
    elsif ($count == 1) {
      printf OUT "%.4g %.4g 0 1\n", $beta, $ave;
    }
    close OUT;
  }
} # Done cycling over order parameters
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Eigenvalues -- one datum per line
# !!! Currently broken by name format change !!!
print "Skipping eigenvalues..."
exit(0);
# Here don't average, but accumulate, sort and sum (density and mode number)
# Trickiness: use output files instead of time series data
# Now grab eigenvalues from all files following cut
# Cut off by minimum of last calculated eigenvalue
my $num;
my $max = -1; # Assume to be the same for all files
my $check;
my $iter = 0;
my $temp;
my @eig;
my @max_eig;
for ($i = $file_cut; $i <= $last_file; $i++) {
  $check = open IN, "< $dir/Out/eig_100_$tag.$i";
  if (!$check) {
    next;
  }
  $iter++;
  @in = <IN>;
  close IN;

  for my $line (@in) {
    if ($line =~ /^EIGENVALUE /) {
      ($junk, $num, $temp)  = split /\s+/, $line;
      push(@eig, sqrt($temp));
      if ($num == $max) {
        push(@max_eig, sqrt($temp));
      }
      elsif ($num > $max) {
        $max = $num;
        $max_eig[0] = sqrt($temp);
      }
    }
  }
}
# Normalize by volume and number of eigenvalue measurements
my $norm = $iter * $L**3 * $Nt;

# Print out sorted eigenvalues and sum, up to min of max
my @sorted = sort { $a <=> $b } @eig;
my @max_sort = sort { $a <=> $b } @max_eig;
my $sum = 0;
open OUT, ">> $dir/results/eig.$tag" or die "Error opening $dir/results/eig.$tag ($!)\n";
for ($i = 0; $i <= $#sorted; $i++) {
#  printf OUT "SORTED: %.4g, MAX %.4g\n", $sorted[$i], $max_sort[$i];
  $sum += 2 / $norm;           # Factor of two for complex conjugate pairs
  if ($sorted[$i] > $max_sort[0]) {
    last;
  }
  printf OUT "%.6g %.8g %.4g\n", $sorted[$i], $sum, $norm;   # Redundant, but easy
}
close OUT;
exit(0);
# ------------------------------------------------------------------
