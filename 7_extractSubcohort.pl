#!/usr/bin/perl

# 18/08/2019
# NTM

# Takes 1 argument: $subcohortFile, a text file with a list of grexomes
# of interest, one per line;
# read on stdin a cohort TSV as produced by extractCohorts.pl, possibly
# filtered, possibly gone through addPatientIDs (but NOT gzipped, just 
# "gunzip -c |" if needed);
# print to stdout a sub-cohort file, identical to the infile except it
# only has the lines where one of the samples of the sub-cohort is HET or HV.
#
# NOTE: the fact that HV et al are the last 6 columns (allowing for 
# one more column with filter logs) is hard-coded.
# If this isn't true the script dies.

use strict;
use warnings;

(@ARGV == 1) || die "needs 1 arg: a text file with one grexomeID per line\n";
my ($subcohortFile) = @ARGV;

#########################################################
# parse subcohort file

# key==grexomeID, value is 1 if grexome is of interest
my %subcohort = ();

(-f $subcohortFile) ||
    die "E: the supplied subcohortFile $subcohortFile doesn't exist or isn't a file\n";
open(SUBC, "$subcohortFile") ||
    die "E: cannot open subcohortFile $subcohortFile for reading\n";
while (my $line = <SUBC>) {
    chomp($line);
    # skip blank lines
    ($line =~ /^\s*$/) && next;
    ($line =~ /^\s*(grexome\d\d\d\d)\s*$/) ||
	die "E: cannot find grexomeID in line $line from subcohort file\n";
    $subcohort{$1} = 1;
}
 

#########################################################
# read cohort file on stdin

# header
my $header = <STDIN>;
chomp($header);
($header =~ /\tHV\tNEGCTRL_HV\tHET\tNEGCTRL_HET\tOTHER\tNEGCTRL_OTHER$/) ||
    ($header =~ /\tHV\tNEGCTRL_HV\tHET\tNEGCTRL_HET\tOTHER\tNEGCTRL_OTHER\tmax_ctrl_hv=[^\t]+$/) ||
    die "E: HV et al are not the last 6 columns of infile header:\n$header\n";
print "$header\n";

# data
 LINE: while (my $line = <STDIN>) {
     chomp($line);
     my @fields = split(/\t/, $line, -1) ;
     foreach my $i ($#fields-5,$#fields-3) {
	 if ($fields[$i]) {
	     ($fields[$i] =~ /^[^~]+~([^~\|]+)$/) || 
		 die "cannot parse HV/HET data $fields[$i]\n";
	     my $samples = $1;
	     foreach my $sample (split(/,/,$samples)) {
		 # no ending $ so we're fine with [DP;AF] and/or if infile went through addPatientIDs.pl
		 ($sample =~ /^(grexome\d\d\d\d)/) ||
		     die "E: inFile has a genotype call for a sample I can't recognize: $sample\n";
		 if ($subcohort{$1}) {
		     print "$line\n";
		     next LINE;
		 }
	     }
	 }
     }
}
