#!/usr/bin/env perl
# EPN, Mon May 18 11:34:28 2015
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

# The definition of $usage explains the script and usage:
my $usage = "\ndnaorg_compare_genomes.pl <file with list of directories> <output directory>\n";
$usage .= "\n"; 
$usage .= " This script compares genomes from the same species based.\n";
$usage .= " mostly on files containing parsed information from their\n";
$usage .= " 'feature tables' which must already exist. Those files are\n";
$usage .= " created by running 'dnaorg_fetch_dna_wrapper.pl -ftable' and\n";
$usage .= " subsequently 'parse-ftable.pl'.\n";
$usage .= "\n";
$usage .= " BASIC OPTIONS:\n";
$usage .= "  -f      : force; if dir <outdir> exists, overwrite it\n";
#$usage .= "  -v      : be verbose; output commands to stdout as they're run\n";
#$usage .= "\n";

my ($seconds, $microseconds) = gettimeofday();
my $start_secs = ($seconds + ($microseconds / 1000000.));
my $executable = $0;
my $be_verbose = 1;

my $out_dir        = undef; 
my $dirlist_file   = undef; 

# initialize variables that can be changed with cmdline options
my $do_force        = 0;     # set to '1' with -f, overwrite output dir if it exists

&GetOptions( "f"        => \$do_force);

if(scalar(@ARGV) != 2) { die $usage; }
($dirlist_file, $out_dir) = (@ARGV);

# store options used, so we can output them 
my $opts_used_short = "";
my $opts_used_long  = "";
if($do_force) { 
  $opts_used_short .= "-f ";
  $opts_used_long  .= "# option:  forcing overwrite of $out_dir dir/file [-f]\n"; 
}

# check for incompatible option combinations:
# NONE YET

###############
# Preliminaries
###############
# check if the $dirlist_file exists
if(! -s $dirlist_file) { die "ERROR $dirlist_file does not exist"; }

# check if our output dir $out_dir exists
if($out_dir !~ m/\/$/) { $out_dir .= "/"; } # add '/'
if(-d $out_dir) { 
  if($do_force) { RunCommand("rm -rf $out_dir", $be_verbose, undef); }
  else          { die "ERROR directory named $out_dir already exists. Remove it, or use -f to overwrite it."; }
}

my $tbl_type = "gene";

########################
# read the dirlist file
########################
open(IN, $dirlist_file) || die "ERROR unable to open $dirlist_file";
my @dir_A = (); # array of the directories
my @len_A = (); # array of lengths
my @tbl_A = (); # array of table files
my $ctr = 0;
while(my $dir = <IN>) { 
  chomp $dir;
  if(! -d $dir) { die "ERROR no directory named $dir exists"; }
  my $dir_tail = $dir;
  $dir_tail =~ s/^.+\///; # remove all but last dir
  my $len = $dir . "/" . $dir_tail . ".length";
  my $tbl = $dir . "/" . $dir_tail . "." . $tbl_type . ".tbl";
  if(! -s $len)                { die "ERROR no length file in ($len does not exist or is empty)"; }
  if($ctr == 0 && (! -s $tbl)) { die "ERROR no table file in first directory ($tbl does not exist or is empty)"; }
  # if $ctr > 0, $tbl may not exist, that's okay
  push(@dir_A, $dir);
  push(@len_A, $len);
  push(@tbl_A, $tbl);
  $ctr++;
}
close(IN);

##################################
# parse the table and length files
##################################
my @values_AHA = ();  # array of hashes of arrays, 
                      # 1D: per-genome listed in $dirlist_file
                      # 2D: key: column name in tbl file
                      # 3D: per-row values for each column
my @totlen_A = (); # lengths of each sequence (in order or @dir_A)

for(my $d = 0; $d < scalar(@dir_A); $d++) { 
  my $len = $len_A[$d];
  parseLength($len, \@totlen_A);

  my $tbl = $tbl_A[$d];
  %{$values_AHA[$d]} = ();
  parseTable($tbl, \%{$values_AHA[$d]});
}

########
# output 
########
# now create the output
# the verbose output
for(my $d = 0; $d < scalar(@dir_A); $d++) { 
  my $dir = $dir_A[$d];
  if(! exists $values_AHA[$d]{"accession"}) { die("ERROR didn't read accession information for for genome number %d\n", $d+1); }
  my $acc = $values_AHA[$d]{"accession"}[0];
  my ($ngenes, $npos, $nneg, $nunc, $nbth, $strand_str) = getGeneStats(\@values_AHA, $d);
  my @genelen_A = ();
  my ($tot_len) = getLengthStats(\@values_AHA, $d, \@genelen_A);
  printf("%-30s  %5d  %5d  %5d  %5d  %5d  %s  %d  ", $acc, $ngenes, $npos, $nneg, $nbth, $nunc, $strand_str, $totlen_A[$d]);
  for(my $i = 0; $i < scalar(@genelen_A); $i++) { 
    printf("  %5d", $genelen_A[$i]);
  }
  printf("\n");
}

printf("\n\n");
# the concise output
my ($ngenes0, $npos0, $nneg0, $nunc0, $nbth0, $strand_str0) = getGeneStats(\@values_AHA, 0);
for(my $d = 0; $d < scalar(@dir_A); $d++) { 
  my $dir = $dir_A[$d];
  my ($ngenes, $npos, $nneg, $nunc, $nbth, $strand_str) = getGeneStats(\@values_AHA, $d);
  my $acc = $values_AHA[$d]{"accession"}[0];
  printf("%-30s  %s%s%s%s%s%s\n", $acc,
         $ngenes     == $ngenes0 ? "*" : "!",
         $npos       == $npos0 ? "*" : "!",
         $nneg       == $nneg0 ? "*" : "!",
         $nbth       == $nbth0 ? "*" : "!",
         $nunc       == $nunc0 ? "*" : "!",
         $strand_str eq $strand_str0 ? "*" : "!");
}


#############
# SUBROUTINES
#############
# Subroutine: RunCommand()
# Args:       $cmd:            command to run, with a "system" command;
#             $be_verbose:     '1' to output command to stdout before we run it, '0' not to
#
# Returns:    amount of time the command took, in seconds
# Dies:       if $cmd fails

sub RunCommand {
  my $sub_name = "RunCommand()";
  my $nargs_exp = 2;

  my ($cmd, $be_verbose) = @_;

  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));
  system($cmd);
  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));

  if($? != 0) { die "ERROR command failed:\n$cmd\n"; }

  return ($stop_time - $start_time);
}

# Subroutine: parseLength()
# Synopsis:   Parses a length file and stores the lengths read
#             into $len_AR.
# Args:       $lenfile: full path to a length file
#             $len_AR:  ref array of lengths
#
# Returns:    void; fills @{$len_AR}
#
# Dies:       if problem parsing $lenfile

sub parseLength {
  my $sub_name = "parseLength()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($lenfile, $len_AR) = @_;

#HM448898.1	2751

  open(LEN, $lenfile) || die "ERROR unable to open $lenfile for reading";

  while(my $line = <LEN>) { 
    chomp $line;
    my ($accn, $length) = split(/\s+/, $line);
    if($length !~ m/^\d+$/) { die "ERROR couldn't parse length file line: $line\n"; } 
    push(@{$len_AR}, $length);
  }
  close(LEN);

  return;
}

# Subroutine: parseTable()
# Synopsis:   Parses a table file and stores the relevant info in it 
#             into $values_HAR.
# Args:       $tblfile:     full path to a table file
#             $values_HAR:  ref to hash of arrays
#
# Returns:    void; fills @{$values_HAR}
#
# Dies:       if problem parsing $tblfile

sub parseTable {
  my $sub_name = "parseTable()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($tblfile, $values_HAR) = @_;

##full-accession	accession	coords	strand	min-coord	gene
#gb|HM448898.1|	HM448898.1	129..476	+	129	AV2

  open(TBL, $tblfile) || die "ERROR unable to open $tblfile for reading";

  # get column header line:
  my $line_ctr = 0;
  my @colnames_A = ();
  my $line = <TBL>;
  my $ncols = undef;
  $line_ctr++;
  if(! defined $line) { die "ERROR did not read any lines from file $tblfile"; }
  chomp $line;
  if($line =~ s/^\#//) { 
    @colnames_A = split(/\s+/, $line);
    $ncols = scalar(@colnames_A);
  }
  else { 
    die "ERROR first line of $tblfile did not start with \"#\"";
  }

  # read remaining lines
  while($line = <TBL>) { 
    chomp $line;
    $line_ctr++;
    if($line =~ m/^\#/) { die "ERROR, line $line_ctr of $tblfile begins with \"#\""; }
    my @el_A = split(/\s+/, $line);
    if(scalar(@el_A) != $ncols) { 
      die "ERROR, read wrong number of columns in line $line_ctr of file $tblfile";
    }
    my $prv_min_coord = 0;
    for(my $i = 0; $i < $ncols; $i++) { 
      my $colname = $colnames_A[$i];
      my $value   = $el_A[$i];
      if($colname eq "min-coord") { 
        if($value < $prv_min_coord) { 
          die "ERROR, minimum coordinates out of order at line $line_ctr and previous line of file $tblfile"; 
        }
        $prv_min_coord = $value; 
        # printf("prv_min_coord: $prv_min_coord\n");
      }
      if(! exists $values_HAR->{$colname}) { 
        @{$values_HAR->{$colname}} = ();
      }
      push(@{$values_HAR->{$colname}}, $el_A[$i]);
    }
  }
  close(TBL);
  return;
}

# Subroutine: getGeneStats()
# Synopsis:   Retreive gene stats for a genome.
# Args:       $values_AHAR:  ref to array of hash of arrays
#             $idx:          first dim index in $values_AHAR to print
#
# Returns:    6 values:
#             $ngenes:     number of genes
#             $npos:       number of genes with all segments on positive strand
#             $nneg:       number of genes with all segmenst on negative strand
#             $nunc:       number of genes with all segments on unknown strand 
#             $nbth:       number of genes with that don't fit above 3 categories
#             $strand_str: strand string, summarizing strand of all genes, in order
#
sub getGeneStats {
  my $sub_name = "getGeneStats()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($values_AHAR, $idx) = @_;

  my $acc;      # accession of this genome
  my $ngenes;   # number of genes in this genome
  my $npos = 0; # number of genes on positive strand 
  my $nneg = 0; # number of genes on negative strand 
  my $nbth = 0; # number of genes with >= 1 segment on both strands (usually 0)
  my $nunc = 0; # number of genes with >= 1 segments that are uncertain (usually 0)
  my $strand_str = "";

  if(! exists $values_AHAR->[$idx]{"strand"})    { die("ERROR didn't read strand information for for genome number %d\n", $idx+1); }

  $acc    = $values_AHAR->[$idx]{"accession"}[0];
  $ngenes = scalar(@{$values_AHAR->[$idx]{"accession"}});
  if ($ngenes > 0) { 
    for(my $i = 0; $i < $ngenes; $i++) { 
      if($values_AHAR->[$idx]{"accession"}[0] ne $acc) { die("ERROR accession mismatch for genome number %d\n", $idx+1); }
      if   ($values_AHAR->[$idx]{"strand"}[$i] eq "+") { $npos++; }
      elsif($values_AHAR->[$idx]{"strand"}[$i] eq "-") { $nneg++; }
      elsif($values_AHAR->[$idx]{"strand"}[$i] eq "!") { $nbth++; }
      elsif($values_AHAR->[$idx]{"strand"}[$i] eq "?") { $nunc++; }
      else { die("ERROR unable to parse strand for gene %d for genome number %s\n", $i+1, $idx+1); }
      $strand_str .= $values_AHAR->[$idx]{"strand"}[$i];
    }
  }

  return ($ngenes, $npos, $nneg, $nunc, $nbth, $strand_str);
}


# Subroutine: getLengthStats()
# Synopsis:   Retreive length stats for a genome, including 
#             the length of all annotated genes.
# Args:       $values_AHAR:  ref to array of hash of arrays
#             $idx:          first dim index in $values_AHAR to print
#             $len_AR:       ref to array of lengths to add to in this subroutine
#
# Returns:    $tot_len:    total length of the genome
#             And @{$len_AR} is filled.
#
sub getLengthStats { 
  my $sub_name = "getLengthStats()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($values_AHAR, $idx, $len_AR) = @_;

  if(! exists $values_AHAR->[$idx]{"strand"})    { die("ERROR didn't read strand information for for genome number %d\n", $idx+1); }

  my $acc    = $values_AHAR->[$idx]{"accession"}[0];
  my $ngenes = scalar(@{$values_AHAR->[$idx]{"accession"}});
  my $tot_len = 0;

  if ($ngenes > 0) { 
    for(my $i = 0; $i < $ngenes; $i++) { 
      if($values_AHAR->[$idx]{"accession"}[0] ne $acc) { die("ERROR accession mismatch for genome number %d\n", $idx+1); }

      my $length = 0;
      if(exists $values_AHAR->[$idx]{"coords"}[$i]) { 
        $length = lengthFromCoords($values_AHAR->[$idx]{"coords"}[$i]);
      }
      push(@{$len_AR}, $length);
    }
  }

  return ($tot_len);
}


# Subroutine: lengthFromCoords()
# Synopsis:   Determine the length of a region give its coords in NCBI format.
#
# Args:       $coords:  the coords string
#
# Returns:    length in nucleotides implied by $coords  
#
sub lengthFromCoords { 
  my $sub_name = "lengthFromCoords()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($coords) = @_;

  my $orig_coords = $coords;
  # Examples:
  # complement(2173412..2176090)
  # complement(join(226623..226774, 226854..229725))

  # remove 'complement('  ')'
  $coords =~ s/^complement\(//;
  $coords =~ s/\)$//;

  # remove 'join('  ')'
  $coords =~ s/^join\(//;
  $coords =~ s/\)$//;

  my @el_A = split(/\s*\,\s*/, $coords);

  my $length = 0;
  foreach my $el (@el_A) { 
    # rare case: remove 'complement(' ')' that still exists:
    $el =~ s/^complement\(//;
    $el =~ s/\)$//;
    $el =~ s/\<//; # remove '<'
    $el =~ s/\>//; # remove '>'
    if($el =~ m/^(\d+)\.\.(\d+)$/) { 
      my ($start, $stop) = ($1, $2);
      $length += abs($start - $stop) + 1;
    }
    else { 
      die "ERROR unable to parse $orig_coords in $sub_name"; 
    }
  }

  # printf("in lengthFromCoords(): orig_coords: $orig_coords returning length: $length\n");
  return $length;
}


