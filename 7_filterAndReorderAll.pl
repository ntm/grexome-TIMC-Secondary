#!/usr/bin/perl


# 26/03/2018
# NTM

# Takes as args: $inDir $outDir $secAnalPath [$pick]
# - $inDir must contain cohort or sample TSVs (possibly gzipped) as 
#   produced by extractCohorts.pl or extractSamples.pl;
# - $outDir doesn't exist, it will be created and filled with one TSV
#   per infile (never gzipped), adding .filtered or .filtered.pick to 
#   the filename;
# - $secAnalPath is the path to the secondary analysis scripts
# - $pick is optional, it's value must be 0 or 1, default is 1, if true
#   we apply --pick as a filter.
# Every infile is filtered and columns are reordered.

use strict;
use warnings;
use POSIX qw(strftime);
use Parallel::ForkManager;


#############################################
## hard-coded stuff that shouldn't change much

# number of parallel jobs to run
my $numJobs = 8;

# names of the scripts that this wrapper actually calls
my ($filterBin,$reorderBin) = ("7_filterVariants.pl","7_reorderColumns.pl");


#############################################
## options / params from the command-line

(@ARGV == 3) || (@ARGV == 4) ||
    die "E: needs 3 or 4 args: an inDir, a non-existant outDir, the path to the secondary analysis scripts, ".
    "and optionally PICK (value 0 or 1, default 1)\n";
my ($inDir, $outDir,$binDir,$pick) = (@ARGV,1);

($pick == 0) || ($pick == 1) ||
    die "E: last arg, if present, must be 0 or 1\n";

(-d $binDir) ||
    die "E: binDir $binDir doesn't exist or isn't a directory\n";
(-f "$binDir/$filterBin") || 
    die "E: binDir $binDir doesn't contain filterBin $filterBin\n";

(-f "$binDir/$reorderBin") || 
    die "E: binDir $binDir doesn't contain reorderBin $reorderBin\n";

(-d $inDir) ||
    die "E: inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "E: cannot opendir inDir $inDir\n";

(-e $outDir) && 
    die "E: found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) ||
    die "cannot mkdir outDir $outDir\n";


#############################################


my $now = strftime("%F %T", localtime);
warn "I: $now - starting to run: ".join(" ", $0, @ARGV)."\n";

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
	warn "W: cannot parse filename of inFile $inFile, skipping it\n";
	$pm->finish;
    }

    my $outFile = $fileStart ;
    ($outFile =~ /\.filtered/) || ($outFile .= ".filtered");
    ($pick) && ($outFile .= ".pick");
    $outFile .= ".csv";
    
    my $com = "perl $binDir/$filterBin --max_ctrl_hv 3 --max_ctrl_het 10 --no_mod";
    # using defaults for AFs 
    # $com .= " --max_af_gnomad 0.01 --max_af_1kg 0.03 --max_af_esp 0.05"
    ($pick) && ($com .= " --pick");
    if ($gz) {
	$com = "gunzip -c $inDir/$inFile | $com ";
    }
    else {
	$com = "cat $inDir/$inFile | $com ";
    }
    
    $com .= " | perl $binDir/$reorderBin ";
    $com .= " > $outDir/$outFile";
    my $now = strftime("%F %T", localtime);
    warn "I: $now - starting $com\n";
    system($com);
    $now = strftime("%F %T", localtime);
    warn "I: $now - Finished $com\n";
    $pm->finish;
}
closedir(INDIR);

$pm->wait_all_children;

$now = strftime("%F %T", localtime);
warn "I: $now - DONE running: ".join(" ", $0, @ARGV)."\n";
