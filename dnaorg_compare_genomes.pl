#!/usr/bin/env perl
# EPN, Mon May 18 11:34:28 2015
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

# The definition of $usage explains the script and usage:
my $usage = "\ndnaorg_compare_genomes.pl <file with list of directories> <output directory>:\n";
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
my @dir_A = ();
my @tbl_A = ();
my $ctr = 0;
while(my $dir = <IN>) { 
  chomp $dir;
  if(! -d $dir) { die "ERROR no directory named $dir exists"; }
  my $tbl = $dir . "/" . $dir . "." . $tbl_type . ".tbl";
  if($ctr == 0 && (! -s $tbl)) { die "ERROR no table file in first directorty ($tbl does not exist or is empty)"; }
  # if $ctr > 0, $tbl may not exist, that's okay
  push(@dir_A, $dir);
  push(@tbl_A, $tbl);
  $ctr++;
}
close(IN);

my @values_AHA = ();  # array of hashes of arrays, 
                      # 1D: per-genome listed in $dirlist_file
                      # 2D: key: column name in tbl file
                      # 3D: per-row values for each column

# parse the tables
my $i = 0;
foreach my $tbl (@tbl_A) {  
  %{$values_AHA[$i]} = ();
  parseTable($tbl, \%{$values_AHA[$i]});
  $i++;
}

# now create the output
for(my $d = 0; $d < scalar(@dir_A); $d++) { 
  my $dir = $dir_A[$d];
  if(! exists $values_AHA[$d]{"accession"}) { 
    die "ERROR no 'accession field for $d\n";
  }
  printf("$dir %d\n", scalar(@{$values_AHA[$d]{"accession"}}));
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

#
# Subroutine: parseTable()
# Synopsis:   Parses a table file and stores the relevant info in it 
#             into $values_HAR.
# Args:       $tblfile:      full path to a table file
#             $values_HAR:  ref to array of hash of arrays
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
        printf("prv_min_coord: $prv_min_coord\n");
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

