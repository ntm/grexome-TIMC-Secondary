#!/usr/bin/perl


# 26/03/2018
# NTM

# Takes as args: $inDir $outDir [$pick]
# - $inDir must contain cohort or sample TSVs (possibly gzipped) as 
#   produced by extractCohorts.pl or extractSamples.pl;
# - $outDir doesn't exist, it will be created and filled with one TSV
#   per infile (never gzipped), adding .filtered or .filtered.pick to 
#   the filename;
# - $pick is optional, it's value must be 0 or 1, default is 1, if true
#   we apply --pick as a filter.
# Every infile is filtered and columns are reordered.

use strict;
use warnings;
use POSIX qw(strftime);
use FindBin qw($RealBin);
use Parallel::ForkManager;


#############################################
## hard-coded stuff that shouldn't change much

# number of parallel jobs to run
my $numJobs = 8;

# names of the scripts that this wrapper actually calls
my ($filterBin,$reorderBin) = ("$RealBin/7_filterVariants.pl","$RealBin/7_reorderColumns.pl");


#############################################
## options / params from the command-line

(@ARGV == 2) || (@ARGV == 3) ||
    die "E $0: needs 2 or 3 args: an inDir, a non-existant outDir, and optionally PICK (value 0 or 1, default 1)\n";
my ($inDir, $outDir, $pick) = (@ARGV,1);

($pick == 0) || ($pick == 1) ||
    die "E $0: last arg, if present, must be 0 or 1\n";

(-f "$filterBin") || 
    die "E $0: filterBin $filterBin not found\n";

(-f "$reorderBin") || 
    die "E $0: reorderBin $reorderBin not found\n";

(-d $inDir) ||
    die "E $0: inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "E $0: cannot opendir inDir $inDir\n";

(-e $outDir) && 
    die "E $0: found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) ||
    die "E $0: cannot mkdir outDir $outDir\n";


#############################################


my $now = strftime("%F %T", localtime);
warn "I $0: $now - starting to run: ".join(" ", $0, @ARGV)."\n";

my $pm = new Parallel::ForkManager($numJobs);

while (my $inFile = readdir(INDIR)) {
    ($inFile =~ /^\./) && next;
    $pm->start && next;
    my ($fileStart,$gz);
    if ($inFile =~ (/^(.+)\.csv$/)) {
	$fileStart = $1;
    }
    elsif ($inFile =~ (/^(.+)\.csv\.gz$/)) {
	$fileStart = $1;
	$gz = 1;
    }
    else {
	warn "W $0: cannot parse filename of inFile $inFile, skipping it\n";
	$pm->finish;
    }

    my $outFile = $fileStart ;
    ($outFile =~ /\.filtered/) || ($outFile .= ".filtered");
    ($pick) && ($outFile .= ".pick");
    $outFile .= ".csv";
    
    my $com = "perl $filterBin --max_ctrl_hv 3 --max_ctrl_het 10 --no_mod";
    # using defaults for AFs 
    # $com .= " --max_af_gnomad 0.01 --max_af_1kg 0.03 --max_af_esp 0.05"
    ($pick) && ($com .= " --pick");
    if ($gz) {
	$com = "gunzip -c $inDir/$inFile | $com ";
    }
    else {
	$com = "cat $inDir/$inFile | $com ";
    }
    
    $com .= " | perl $reorderBin ";
    $com .= " > $outDir/$outFile";
    my $now = strftime("%F %T", localtime);
    warn "I $0: $now - starting $com\n";
    system($com);
    $now = strftime("%F %T", localtime);
    warn "I $0: $now - Finished $com\n";
    $pm->finish;
}
closedir(INDIR);

$pm->wait_all_children;

$now = strftime("%F %T", localtime);
warn "I: $now - DONE running: ".join(" ", $0, @ARGV)."\n";
