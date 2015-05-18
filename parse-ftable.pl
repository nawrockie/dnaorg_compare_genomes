#!/usr/bin/env perl
# EPN, Mon May 11 10:32:44 2015
#
use strict;
use warnings;
use Getopt::Long;

## Example 'ftable':, my comments added after the second # in each line
##
#
##################################################
#>Feature ref|NC_001359.1|                          # accession line, denotes a new DB record
#230	985	gene                                # coords_feature line, gives coordinates and feature name
#			gene	AR1                 # quals line, gives qualifier and (usually) value
#			locus_tag	PhyvvsAgp1  # ditto
#			db_xref	GeneID:988143       # ditto
#230	985	CDS                                 # coords_feature line, gives coordinates and feature name
#1      102                                         # coords_only line, gives only coordinates, feature was on last seen coords_feature line
#			product	hypothetical protein 
#			prot_desc	AR1
#			transl_table	1
#			protein_id	ref|NP_040320.1|
#			db_xref	GOA:Q06912
#			db_xref	InterPro:IPR000263
#			db_xref	InterPro:IPR000650
#			db_xref	UniProtKB/Swiss-Prot:Q06912
#
#>Feature ref|NC_001265.2|
#70	2679	gene
#			locus_tag	CarMVgp1
#			db_xref	GeneID:1724753
##################################################
# For each unique 'feature' name ('gene' and 'CDS' in the example above), this 
# script will create a table that contains all the information for that feature
# in all accessions in the ftable file.
# 
# For example:
#
# --------------------------------------------------------------------
# > perl parse-ftable.pl example.ftable example
# Output 1 data lines to file example.CDS.tbl.
# Output 2 data lines to file example.gene.tbl.
# --------------------------------------------------------------------
#
# > cat example.CDS.tbl
##full-accession	accession	coords	product	prot_desc	transl_table	protein_id	db_xref
#ref|NC_001359.1|	NC_001359.1	join(230..985,1..102)	hypothetical protein 	AR1	1	ref|NP_040320.1|	GOA:Q06912;;InterPro:IPR000263;;InterPro:IPR000650;;UniProtKB/Swiss-Prot:Q06912
# 
# > cat example.gene.tbl 
##full-accession	accession	coords	gene	locus_tag	db_xref
# ref|NC_001359.1|	NC_001359.1	230..985	AR1	PhyvvsAgp1	GeneID:988143
# ref|NC_001265.2|	NC_001265.2	70..2679	-	CarMVgp1	GeneID:1724753
#
###################################################
# The definition of $usage explains the script and usage:
my $usage = "\nparse-ftable.pl:\n";
$usage .= "\n"; 
$usage .= " Create separate tables for each feature listed in feature table\n";
$usage .= " output from dnaorg-fetch-dna-wrapper.pl\n";
$usage .= "\n";
$usage .= "\nusage: parse-ftable.pl <ftable file> <root for naming output files>\n";
$usage .= "\n";
$usage .= " OPTIONS:\n";
$usage .= "  -f <s>: only output table for feature <s>\n";
$usage .= "  -d <s>: create output files in directory <s>, not the cwd\n";
$usage .= "\n";

# initialize variables that can be changed with cmdline options
my $only_feature     = undef; # set to a value if -f used
my $cap_only_feature = undef; # capitalized version of $only_feature
my $do_only_feature  = 0;     # set to '1' below if -f used
my $out_dir          = undef; # set to a value if -d used
&GetOptions( "f=s" => \$only_feature, 
             "d=s" => \$out_dir);

if(scalar(@ARGV) != 2) { die $usage; }
my ($ftable_in, $out_root) = (@ARGV);

open(IN, $ftable_in) || die "ERROR unable to open $ftable_in for reading."; 

if(defined $only_feature) { 
  $cap_only_feature = $only_feature;
  $cap_only_feature =~ tr/a-z/A-Z/;
  $do_only_feature = 1; 
}

# define other variables we'll use
my %feature_H       = ();      # keys: name of features; we'll output a separate table for each
my $faccn           = undef;   # full accession, e.g. ref|NC_000883.2|
my $fac             = undef;   # full accession and coords, accession and each segment coords separated by $fac_sep
my $faccn2          = undef;   # another full accession
my $accn            = undef;   # full accession with extraneous characters removed: e.g. NC_000883.2
my $start           = undef;   # start position
my $stop            = undef;   # stop position
my $qname           = undef;   # a qualifier name,  e.g. 'organism'
my $qval            = undef;   # a qualifier value, e.g. 'Paramecia bursaria Chlorella virus 1'
my $feature         = undef;   # a feature name, e.g. "CDS", or "gene"
my $cap_feature     = undef;   # a feature name capitalized, necessary so -f option can be case insensitive
my $coords          = undef;   # coordinates
my $qnqv_sep        = "!!";    # string that separates 'qualifier name' and and 'qualifier_values' 
my $fac_sep         = "::";    # strings separating full accessions, start and stop coords
my $qval_sep        = ";;";    # strings separating qualifier values  
my $line_ctr        = 0;       # count of number of lines read in ftable
my $accn_ctr        = 0;       # count of number of accessions read

# more complicated data structures:
my %quals_HHA       = ();      # values that will go into a table
                               # key 1: feature name, e.g. "CDS"
                               # key 2: 'fac' string: <full_accession><fac_seq><coordinates>
                               # value: array of 'qnqv' strings: <qualifier_name><qnqv_sep><qualifier_value>
my @faccn_A         = ();      # array of all full accessions read
my %fac_HHA         = ();      # used to easily determine list of 2D keys ('fac's) in quals_HHA 
                               # for a given feature and faccn.
                               # key 1: feature name: e.g. "CDS"
                               # key 2: 'faccn', full accession
                               # value: array of 'fac' strings: <full_accession><fac_seq><coordinates>
my %faccn2accn_H    = ();      # key: full accession, value: short accession, used for convenience
my %column_HA       = ();      # key: feature, e.g. "CDS"
                               # value: array of qualifiers that exist for that feature in at least 1 faccn
my %column_HH       = ();      # existence 2D hash, used to aide construction of %column_HA
                               # key 1: feature, e.g. "CDS"
                               # key 2: qualifier
                               # value: '1' (always)

# variables used to make sure input feature table is in 
# expected format w.r.t to order of line types
my $prv_was_accn           = 0; # set to '1' if previous line was an accession line
my $prv_was_coords_feature = 0; # set to '1' if previous line was a coordinates line with a feature name
my $prv_was_coords_only    = 0; # set to '1' if previous line was a coordinates line without a feature name
my $prv_was_quals          = 0; # set to '1' if previous line was a qualifier_name qualifier value line

#########
# INPUT # 
#########
while(my $line = <IN>) { 
  $line_ctr++;
  chomp $line;
  if($line =~ m/\w/) { 
    # parse each of the 4 line types differently
    # -------------------------------------------------------
    if($line =~ /Feature\s+(\S+)$/) { 
      # ACCESSION LINE
      # example:
      #>Feature ref|NC_001359.1|    
      $faccn = $1;
      $accn_ctr++;
      # does line order make sense?
      if($accn_ctr == 1   ||          # first accession of the file
         $prv_was_quals ||          # previous line was a quals line (common)
         $prv_was_coords_feature || # previous line was a coords feature line (rare)
         $prv_was_coords_only    || # previous line was a coords line without a feature (rare)
         $prv_was_accn) {           # previous line was an accession line (rare)
        # line order makes sense, keep going...

        # determine short version of the accession, e.g. NC_001359 in above example
        if($faccn =~ /[^\|]*\|([^\|]*)\|/) { 
          $accn = $1;
          $faccn2accn_H{$faccn} = $accn;
          push(@faccn_A, $faccn);
        }
        else { 
          die "ERROR unable to parse Feature line $line"; 
        }
        $feature     = undef; 
        $cap_feature = undef;
        $fac         = undef;

        # update '$prv_*' values that we use to make sure line order makes sense
        $prv_was_accn           = 1;
        $prv_was_coords_feature = 0;
        $prv_was_coords_only    = 0;
        $prv_was_quals          = 0;
        #printf("set prv_was_accn\n");
      }
      else { # line order is unexpected
        die "ERROR unexpected line order (accession line) at line $line_ctr: $line\n";
      }
    }
    # -------------------------------------------------------
    elsif($line =~ /^(\<?\d+\^?)\s+(\>?\d+)\s+(\S+)$/) { 
      # COORDINATES LINE WITH A FEATURE NAME
      # example:
      #230	985	gene
      ($start, $stop, $feature) = ($1, $2, $3);

      # does line order make sense?
      if($prv_was_accn  ||           # previous line was accession line (common)
         $prv_was_quals ||           # previous line was quals line (common)
         $prv_was_coords_feature ||  # previous line was coords line with a feature (rare)
         $prv_was_coords_only)    {  # previous line was coords line without a feature (rarer)
        # line order makes sense, keep going...
        $cap_feature = $feature;
        $cap_feature =~ tr/a-z/A-Z/;
        $fac = $faccn . $fac_sep . $start . $fac_sep . $stop;

        # update '$prv_*' values that we use to make sure line order makes sense
        $prv_was_accn           = 0;
        $prv_was_coords_feature = 1;
        $prv_was_coords_only    = 0;
        $prv_was_quals          = 0;
        #printf("set prv_was_coords_feature\n");
      }
      else { # line order is unexpected
        die "ERROR unexpected line order (coords_feature) at line $line_ctr: $line\n";
      }
    }
    # -------------------------------------------------------
    elsif($line =~ /^(\<?\d+)\s+(\>?\d+)$/) { 
      # COORDINATES LINE WITHOUT A FEATURE NAME
      # example:
      #1	54
      ($start, $stop) = ($1, $2);

      # does line order make sense?
      if($prv_was_coords_feature || # previous line was a coords line with a feature (common)
         $prv_was_coords_only) {    # previous line was a coords line without a feature (common)
        # line order makes sense, keep going...

        $fac .= $fac_sep . $start . $fac_sep . $stop;
        
        # update '$prv_*' values that we use to make sure line order makes sense
        $prv_was_accn           = 0;
        $prv_was_coords_feature = 0;
        $prv_was_coords_only    = 1;
        $prv_was_quals          = 0;
        #printf("set prv_was_coords_only\n");
      }
      else { # line order is unexpected
        die "ERROR unexpected line order (coords_only line) at line $line_ctr: $line\n";
      }
    }
    # -------------------------------------------------------
    elsif(($line =~ /^\s+\S+\s+.+$/) || 
          ($line =~ /^\s+\S+$/)) { 
      # QUALIFIER LINE
      # examples:
      #			gene	AR1
      #			locus_tag	PhyvvsAgp1

      # before parsing it, do two sanity checks
      if(! defined $fac)     { die "ERROR coordinates undefined at a qualifier line"; }
      if(! defined $feature) { die "ERROR didn't read feature line before line: $line\n"; }
      # and determine if we even care about this feature
      if((! $do_only_feature) || ($cap_feature eq $cap_only_feature)) { 
        # does line order make sense?
        if($prv_was_coords_feature || 
           $prv_was_coords_only    ||
           $prv_was_quals) { 
          # line order makes sense, keep going...


          if(! $prv_was_quals) { 
            # first quals line for this feature
            # at this point, we know that we have the full coordinates for the feature
            # so initialize the information
            if(! exists $fac_HHA{$feature}) {
              %{$fac_HHA{$feature}} = (); 
            }
            if(! exists $fac_HHA{$feature}{$faccn}) { 
              @{$fac_HHA{$feature}{$faccn}} = ();
            }
            push(@{$fac_HHA{$feature}{$faccn}}, $fac);
            # printf("feature: $feature\n");
            if(! exists $feature_H{$feature}) { 
              if((! $do_only_feature) || 
                 ($cap_feature eq $cap_only_feature)) { 
                $feature_H{$feature} = 1;
              }
            }
            if(! exists $quals_HHA{$feature}) { 
              %{$quals_HHA{$feature}} = ();
              @{$quals_HHA{$feature}{$fac}} = ();
            }
          } # end of 'if(! $prv_was_quals)'

          # now parse the line;
          # examples:
          #			gene	AR1
          #			locus_tag	PhyvvsAgp1
          if($line =~ /^\s+(\S+)\s+(.+)$/) { 
            ($qname, $qval) = ($1, $2);
          }
          elsif($line =~ /^\s+(\S+)$/) { 
            $qname = $1;
            $qval = "<no_value>";
          }
          else { 
            die "ERROR didn't parse quals line on second pass: $line\n"; 
          }
          if($qname =~ m/\Q$qnqv_sep/)   { die "ERROR qualifier_name $qname has the string $qnqv_sep in it"; }
          if($qval  =~ m/\Q$qnqv_sep/)   { die "ERROR qualifier_value $qval has the string $qnqv_sep in it"; }
          my $qnqv = $qname . $qnqv_sep . $qval; # this is how we store the qualifier name and value, as a concatenated string in values_HHA
          push(@{$quals_HHA{$feature}{$fac}}, $qnqv);

          # and update the column data structures which just keep info on names and order of columns
          if(! exists $column_HH{$feature}) { 
            %{$column_HH{$feature}} = ();
            @{$column_HA{$feature}} = (); 
          }
          if(! exists $column_HH{$feature}{$qname}) { 
            push(@{$column_HA{$feature}}, $qname);
            $column_HH{$feature}{$qname} = 1;
          }
        }
        else { # unexpected line order
          die "ERROR unexpected line order (quals line) at line $line_ctr: $line\n";
        }          
      } # end of 'if((! $do_only_feature) || ($cap_feature ne $cap_only_feature))'
      # update '$prv_*' values that we use to make sure line order makes sense
      $prv_was_accn           = 0;
      $prv_was_coords_feature = 0;
      $prv_was_coords_only    = 0;
      $prv_was_quals          = 1;
      #printf("set prv_was_quals\n");
    }
    # -------------------------------------------------------
    else { 
      die "ERROR unable to parse line $line_ctr: $line\n"; 
    }
    # -------------------------------------------------------
  }
}
close(IN);

##########
# OUTPUT # 
##########
my $column;
my $strand; 
my $sort_coord;
# one table per feature
foreach $feature (sort keys %feature_H) { 
  my $cur_line_ctr = 0;
  # name the file
  my $tblout;
  if(defined $out_dir) { 
    $tblout = $out_dir . "/" . $out_root . "." . $feature . ".tbl";
  }
  else { 
    $tblout = $out_root . "." . $feature . ".tbl";
  }

  open(OUT, ">" . $tblout) || die "ERROR unable to open $tblout for writing"; 

  # print header line with names of column
  print OUT "#full-accession\taccession\tcoords\tstrand\tfirst-start";
  foreach $column (@{$column_HA{$feature}}) { 
    print OUT "\t$column"; 
  }
  print OUT "\n";

  # go through all full accessions in order they were read from feature table
  foreach $faccn (@faccn_A) { 
    if(exists $fac_HHA{$feature}{$faccn}) { # if this accession has >= 1 qualifiers for this feature
      foreach $fac (@{$fac_HHA{$feature}{$faccn}}) { # foreach 'fac', accession + set of coords
        my @output_A = ();
        ($faccn2, $coords, $sort_coord, $strand) = breakdownFac($fac, $fac_sep);
        if($faccn ne $faccn2) { die "ERROR inconsistent fac value: $faccn ne $faccn2"; }

        if(exists $quals_HHA{$feature}{$fac}) { # if there's any qualifiers for this fac
          # printf("values_HHA feature: $feature fac: $fac exists!\n"); 
          print OUT $faccn. "\t" . $faccn2accn_H{$faccn} . "\t" . $coords . "\t" . $strand . "\t" . $sort_coord;

          # for all columns in the table
          foreach $column (@{$column_HA{$feature}}) {
            ### printf("\n\n");
            my $column_str = ""; 

            # for all qualifier names and values 
            foreach my $qnqv (@{$quals_HHA{$feature}{$fac}}) { 
              ($qname, $qval) = split($qnqv_sep, $qnqv);
              ### printf("faccn: $faccn qnqv: $qnqv split into $qname $qval\n");

              # if this qname matches this column, then it's the appropriate value to output here
              if($qname eq $column) { 
                if($column_str eq "") { # first value in this cell
                  $column_str = $qval;  
                }
                else { 
                  if($qval =~ m/\Q$qval_sep/) { die "ERROR qualifier_name $qval has the string $qval_sep in it"; }
                  $column_str .= $qval_sep . $qval; # not first value, concatenate onto previous values
                }
              }
              ### printf("\tcolumn_str: $column_str\n");
            }
            # if there's no value for this qualifier, put '-'
            if($column_str eq "") { $column_str = "-"; }
            print OUT "\t$column_str";
          }
          print OUT "\n";
          $cur_line_ctr++;
        }
      }
    }
  }
  close(OUT);
  printf("Output %6d data lines to file $tblout.\n", $cur_line_ctr);
}
  

#############
# SUBROUTINES
#############
#
# Subroutine: breakdownFac()
# Args:       $fac:            full accession and coordinates, concatenated together
#             $fac_sep:        character in between each token in $fac
#
# Returns:    Three values:
#             $faccn:       full accession
#             $ncbi_coords: NCBI format for coordinates
#             $sort_coord:  start position of first segment, possibly useful for ordering multiple features
#             $strand:      '+' if all segments are on fwd strand
#                           '-' if all segments are on neg strand
#                           '?' if all segments are 1 nucleotide (strand is consequently uncertain)
#                           '!' if >= 1 segment on two or more of following: fwd strand, rev strand, uncertain
#
# Dies:       if there's not an odd number of tokens

sub breakdownFac {
  my $sub_name = "breakdownFac()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($fac, $fac_sep) = @_;

  my @elA = split($fac_sep, $fac);
  my $nel = scalar(@elA);
  if($nel % 2 == 0) { die "ERROR unexpectedly there's an even number of tokens in breakdownFac(): fac: $fac\n"; }

  my $faccn = $elA[0];
  my ($istart, $istop) = (1, 2);
  my $is_rev = 0; # true if current segment is on reverse strand
  my $nfwd   = 0; # number of segments that are on forward strand
  my $nrev   = 0; # number of segments that are on reverse strand
  my $nunc   = 0; # number of segments for which strand is uncertain (start == stop)
  my $some_revcomp  = undef;
  my @ncbi_coords_A = (); # [0..$nsegments-1] strings to concatenate to make ncbi coords
  my @is_rev_A      = (); # [0..$nsegments-1] '1' if current segment is on reverse strand
  my $nsegments     = 0;
  my $have_carrot   = 0; # set to '1' if $start has a '^' at the end, e.g. 804^ 805
  my $sort_coord = undef; # start position of first segment
  while($istop <= $nel) { 
    my ($orig_start, $orig_stop) = ($elA[$istart], $elA[$istop]);
    $start = $orig_start;
    $stop  = $orig_stop;
    if($start =~ s/\^$//) { 
      $have_carrot = 1;
    }
    $start =~ s/^<//;
    $start =~ s/^>//;
    $stop  =~ s/^<//;
    $stop  =~ s/^>//;
    if($istart == 1) { $sort_coord = $start; }

    if($have_carrot) { 
      if(abs($start-$stop) > 1) { die "ERROR found carrot but range is not exactly 1 nt: fac: $fac\n"; }
    }

    if($start == $stop) { # special case, we can't tell if we're in reverse complement or not
      $nunc++;
      push(@is_rev_A, 0);
    }
    elsif($start > $stop) {
      $nrev++;
      $is_rev = 1;
      push(@is_rev_A, 1);
    }
    else { # $start < $stop, not reverse complement
      $nfwd++;
      push(@is_rev_A, 0);
    }
    if($have_carrot) { 
      push(@ncbi_coords_A, $orig_start . $orig_stop);
    }
    else { 
      if($is_rev) { 
        push(@ncbi_coords_A, $orig_stop . ".." . $orig_start);
      }
      else { 
        push(@ncbi_coords_A, $orig_start . ".." . $orig_stop);
      }      
    }
    $nsegments++;
    $istart += 2; 
    $istop  += 2;
  }
  if($have_carrot) { 
    if($nsegments > 1) { die "ERROR found carrot but more than one segment: fac: $fac\n"; }
  }

  my $ncbi_coords = "";
  if($nfwd > 0 && $nrev > 0) { 
    # special case, we need to put 'complement(' ')' around each rev segment separately
    # e.g: join(161990..162784,complement(88222..88806),complement(86666..87448))
    for(my $i = 0; $i < $nsegments; $i++) { 
      if($i > 0) { $ncbi_coords .= ","; }
      if($is_rev_A[$i]) { $ncbi_coords .= "complement(" . $ncbi_coords_A[$i] . ")" }
      else              { $ncbi_coords .= $ncbi_coords_A[$i]; }
    }
  }
  else { 
    for(my $i = 0; $i < $nsegments; $i++) { 
      if($i > 1) { $ncbi_coords .= ","; }
      $ncbi_coords .= $ncbi_coords_A[$i];
    }
    if($nrev > 0) { # if we get here all segments are on reverse strand
      $ncbi_coords = "complement(" . $ncbi_coords . ")";
    }
  }
  if($nsegments > 1) { 
    $ncbi_coords = "join(" . $ncbi_coords . ")";
  }
  # printf("in breakdownFac() input: $fac, returning $faccn $coords ncbi_coords:$ncbi_coords\n");

  # determine strand
  my $ret_strand = undef;
  if   ($nfwd == $nsegments) { $ret_strand = "+"; }
  elsif($nrev == $nsegments) { $ret_strand = "-"; }
  elsif($nunc == $nsegments) { $ret_strand = "?"; }
  else                       { $ret_strand = "!"; }

  return($faccn, $ncbi_coords, $sort_coord, $ret_strand);
}
