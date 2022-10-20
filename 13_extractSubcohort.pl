#!/usr/bin/perl

# 18/08/2019
# NTM

# Takes 1 argument: $subcohortFile, a text file with a list of samples
# of interest, one per line;
# read on stdin one of:
# - a cohort TSV as produced by extractCohorts.pl, probably filtered
#   and reordered, probably gone through addPatientIDs and requireUndiagnosed
#   (but NOT gzipped, just "gunzip -c |" if needed);
# - a transcript TSV as produced by extractTranscripts.pl, probably gone
#   through addPatientIDs.
# "Probably" here means that this script should work without, but I
# currently test and use it that way.
# We automatically detect whether stdin is a cohort or a transcript file.
# Print to stdout a sub-cohort file, identical to the infile except it
# only has the lines where one of the samples of the sub-cohort is:
# - Cohorts -> HV or HET (OTHERCAUSE is not enough)
# - Transcripts -> HV or BIALLELIC (including OTHERCAUSE)


use strict;
use warnings;
use File::Basename qw(basename);
use POSIX qw(strftime);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


(@ARGV == 1) || die "E $0, needs 1 arg: a text file with one sampleID per line\n";
my ($subcohortFile) = @ARGV;

my $now = strftime("%F %T", localtime);
my $subCpretty = basename($subcohortFile);
warn "I $now: $0 - starting to run with $subCpretty, don't yet know if STDIN is cohort or transcripts file\n";


#########################################################
# parse subcohort file

# key==sampleID, value is 1 if sample is of interest
my %subcohort = ();

(-f $subcohortFile) ||
    die "E $0: the supplied subcohortFile $subcohortFile doesn't exist or isn't a file\n";
open(SUBC, "$subcohortFile") ||
    die "E $0: cannot open subcohortFile $subcohortFile for reading\n";
while (my $line = <SUBC>) {
    chomp($line);
    # skip blank lines
    ($line =~ /^\s*$/) && next;
    # remove leading or trailing blank chars and make sure we have a reasonable ID (no blanks,
    # no ( or [)
    ($line =~ /^\s*(\S+)\s*$/) ||
	die "E $0: cannot find reasonable sampleID in line $line from subcohort file\n";
    ($1 =~ /[(\[]/) && 
	die "E $0: sampleID $1 from subcohort file contains ( or [, illegal\n";
    $subcohort{$1} = 1;
}


#########################################################
# read cohort file on stdin and detect type

# header
my $header = <STDIN>;
chomp($header);

# detect type:
# $type == 1 if stdin is a cohort file, 2 if it's a transcript file
my $type = 0;
if ($header =~ /\sCOUNTSAMPLES_/) {
    # transcripts
    $type = 2;
}
if ($header =~ /\sCOUNT_/) {
    ($type) && die "E $0: stdin is detected both as a transcripts and a cohorts file...\n";
    $type = 1;
}
($type) || die "E $0: stdin doesn't seem to be either a cohorts or a transcripts file\n";


my @header = split(/\t/,$header);
# indexes of columns we want to look at
my @genoCols = ();
foreach my $i (0..$#header) {
    if ($type == 1) {
	# start by ignoring COUNT*, NEGCTRL, COMPAT, and OTHERCAUSE
	if (($header[$i] =~ /^COUNT_/) || ($header[$i] =~ /^NEGCTRL_/) ||
	    ($header[$i] =~ /^COMPAT_/) || ($header[$i] =~ /_OTHERCAUSE_/)) {
	    next;
	}
	# remaining _HV or _HET column should be good
	elsif (($header[$i] =~ /_HV$/) || ($header[$i] =~ /_HET$/)) {
	    push(@genoCols, $i);
	}
    }
    elsif ($type==2) {
	# transcripts
	if (($header[$i] =~ /^HV_/) || ($header[$i] =~ /^BIALLELIC_/) ||
	    ($header[$i] =~ /^OTHERCAUSE_/)) {
	    push(@genoCols, $i);
	}
    }
}
# we should have found exactly 2 columns (HV, HET) for cohorts
# or 7 for transcripts
(($type==1) && (@genoCols == 2)) ||
    (($type==2) && (@genoCols == 7)) ||
    die "E $0: type==$type but can't find the right number of geno headers in:\n$header\n";

print "$header\n";

# data
 LINE: while (my $line = <STDIN>) {
     chomp($line);
     my @fields = split(/\t/, $line, -1) ;
     foreach my $i (@genoCols) {
	 if ($fields[$i]) {
	     # clean up: remove genotype eg '1/1~' in cohortfiles
	     if ($type==1) {
		 ($fields[$i] =~ s/^[^~]+~([^~\|]+)$/$1/) || 
		     die "E $0: cohortfile but cannot parse genoData $fields[$i]\n";
	     }
	     foreach my $sample (split(/,/,$fields[$i])) {
		 # valid sampleIDs cannot have '(' or '[' ,  so with the below regexp we're fine 
		 # with [DP:AF] / [BF:RR] and/or if infile went through addPatientIDs.pl
		 ($sample =~ /^([^(\[]+)/) ||
		     die "E $0: inFile has a genoData for a sample I can't recognize: $sample\n";
		 if ($subcohort{$1}) {
		     print "$line\n";
		     next LINE;
		 }
	     }
	 }
     }
}

$now = strftime("%F %T", localtime);
my $typeString = "Cohort";
($type == 2) && ($typeString = "Transcripts");
warn "I $now: $0 - ALL DONE, completed successfully with $subCpretty on a $typeString file!\n";
