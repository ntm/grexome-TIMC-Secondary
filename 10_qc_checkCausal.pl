#!/usr/bin/perl

# 25/06/2021
# NTM

# QC script: look for severe biallelic variants affecting  the canonical
# transcript of each "causal gene" from the samples metadata file.


use strict;
use warnings;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use Getopt::Long;
use POSIX qw(strftime);

use lib "$RealBin";
use grexome_metaParse qw(parseSamples);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## options / params from the command-line

# samples XLSX file, no default
my $samplesFile = "";

# dir containing samples CSV files produced by grexome-TIMC-Secondary,
# default to current dir
my $inDir = ".";

# help: if true just print $USAGE and exit
my $help = '';


my $USAGE = "\nFor each sample that has a causal gene in the metadata XLSX, examine the SAMPLES 
results CSV files in indir and report if and how the causal gene's canonical transcript is hit.
Print our findings to stdout.
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--samplesFile : samples metadata xlsx file, with path
--indir [$inDir] : dir containing samples CSV files produced by grexome-TIMC-Secondary
--help : print this USAGE";

GetOptions ("samplesFile=s" => \$samplesFile,
	    "indir=s" => \$inDir,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($samplesFile) || die "E $0: you must provide a samplesFile file\n$USAGE\n";
(-f $samplesFile) || die "E $0: the supplied samplesFile file doesn't exist\n";

(-d $inDir) ||
    die "E $0: inDir $inDir doesn't exist or isn't a directory\n$USAGE\n";
opendir(INDIR, $inDir) ||
    die "E $0: cannot opendir inDir $inDir\n";

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


#########################################################
# parse useful info from samples metadata file:
# hashref, key==sample id, value is the $cohort this sample belongs to (ie pathologyID)
my $sample2cohortR;
# hashref, key==sample id, value == HGNC gene name of the known causal gene
my $sample2causalR;

{
    my @parsed = &parseSamples($samplesFile);
    $sample2cohortR = $parsed[0];
    $sample2causalR = $parsed[3];
}

#########################################################
# examine each sample CSV file

# total number of samples
my $nbSamples = scalar(keys(%$sample2cohortR));
# number of samples with a causal gene in metadata
my $nbCausal = scalar(keys(%$sample2causalR));
# number of samples whose "causal gene" is hit by biallelic HIGH / MODHIGH
my ($nbCausalH,$nbCausalMH) = (0,0);
# arrays of strings to print (grep commands for investigating) for samples with:
# only LOW or MODERATE biallelic variants / NO biallelic variants / NO variants at all
my @noMH = ();
my @noBA = ();
my @noVar = ();

foreach my $inFile (sort(readdir(INDIR))) {
    ($inFile =~ /^\./) && next;
    my $sample;
    foreach my $s (keys(%$sample2cohortR)) {
	if ($inFile =~ /^$sample2cohortR->{$s}\.$s\./) {
	    $sample = $s;
	    last;
	}
    }
    ($sample) || 
	((print "W: inFile $inFile doesn't seem to correspond to any sample, skipping it\n") && next);

    # empty $sample2cohortR as we go, for sanity testing and speed
    delete($sample2cohortR->{$sample});

    # if no known causal gene for sample $sample, skip $inFile
    (defined($sample2causalR->{$sample})) || next;

    my $causal = $sample2causalR->{$sample};

    # find lines affecting a CANONICAL transcript of gene $causal in $inFile
    open(INFILE, "$inDir/$inFile") ||
	die "E $0: cannot open inFile $inDir/$inFile\n";

    # indexes of columns of interest
    my ($symbolCol, $biallelCol, $canonCol) = (-1,-1,-1);
    # header: grab column indexes of interest
    my $header = <INFILE>;
    chomp($header);
    my @header = split(/\t/,$header);
    foreach my $i (0..$#header) {
	($header[$i] eq 'SYMBOL') && ($symbolCol = $i);
	($header[$i] eq 'BIALLELIC') && ($biallelCol = $i);
	($header[$i] eq 'CANONICAL') && ($canonCol = $i);
	($biallelCol > -1) && ($symbolCol > -1) && ($canonCol > -1) && last;
    }
    (($biallelCol > -1) && ($symbolCol > -1) && ($canonCol > -1)) || 
	die "E $0: cannot find required column headers in $inFile:\n$header\n";

    # find lines affecting a CANONICAL transcript of $causal in $inFile,
    # if no lines $foundCausal stays false
    my $foundCausal = 0;
    while (my $line = <INFILE>) {
	chomp($line);
	my @line = split(/\t/,$line);
	($line[$symbolCol] eq "' $causal") || next;
	($line[$canonCol] eq 'YES') || next;

	if ($line[$biallelCol] eq 'HIGH') {
	    $nbCausalH++;
	}
	elsif ($line[$biallelCol] eq 'MODHIGH') {
	    $nbCausalMH++;
	}
	elsif ($line[$biallelCol] eq 'NO') {
	    push(@noBA, "grep -P ' $causal\\t' $inDir/$inFile");
	}
	else {
	    push(@noMH, "grep -P ' $causal\\t' $inDir/$inFile");
	}
	$foundCausal = 1;
	last;
    }
    close(INFILE);
    ($foundCausal) ||
	push(@noVar, "grep -P ' $causal\\t' $inDir/$inFile");
}

closedir(INDIR);

# sanity: every sample from metadata should have been seen
(keys(%$sample2cohortR)) &&
    (print "W: some samples from metadata have no samples.csv files! ".join(" ",sort(keys(%$sample2cohortR)))."\n\n");

print "Total number of samples in metadata: $nbSamples\n";
print "Samples with a causal gene in metadata: $nbCausal\n\n";
print "Examining the samples CSV files in $inDir (canonical transcripts only), we found:\n";
print "samples whose causal gene is BIALLELIC HIGH: $nbCausalH\n";
print "samples whose causal gene is BIALLELIC MODHIGH: $nbCausalMH\n";
print "samples whose causal gene is BIALLELIC MODERATE or LOW: ".scalar(@noMH)."\n";
print "samples whose causal gene is only hit by MONOALLELIC: ".scalar(@noBA)."\n";
print "samples whose causal gene has NO VARIANT: ".scalar(@noVar)."\n";
print"\n";

if (@noMH) {
    print "To investigate the samples whose causal gene is BIALLELIC MODERATE or LOW:\n";
    print join("\n", @noMH);
    print "\n\n";
}

if (@noBA) {
    print "To investigate samples whose causal gene is only hit by MONOALLELIC:\n";
    print join("\n", @noBA);
    print "\n\n";
}

if (@noVar) {
    print "To confirm samples whose causal gene has NO VARIANT:\n";
    print join("\n", @noVar);
    print "\n\n";
}

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";
