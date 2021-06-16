#!/usr/bin/perl


# 26/03/2018
# NTM

# Every infile from $inDir is filtered and columns are reordered,
# results go into $outDir (which must not pre-exist)
# - $inDir must contain cohort or sample TSVs (possibly gzipped) as 
#   produced by extractCohorts.pl or extractSamples.pl;
# - $outDir doesn't exist, it will be created and filled with one TSV
#   per infile (never gzipped), adding .filtered or .filtered.pick to 
#   the filename.

use strict;
use warnings;
use File::Basename qw(basename);
use POSIX qw(strftime);
use FindBin qw($RealBin);
use Parallel::ForkManager;
use Getopt::Long;

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## hard-coded stuff that shouldn't change much

# names of the scripts that this wrapper actually calls
my ($filterBin,$reorderBin) = ("$RealBin/7_filterVariants.pl","$RealBin/7_reorderColumns.pl");


# check hard-coded stuff
(-f "$filterBin") || 
    die "E $0: filterBin $filterBin not found\n";
(-f "$reorderBin") || 
    die "E $0: reorderBin $reorderBin not found\n";


#############################################
## options / params from the command-line

my ($inDir, $outDir);

# number of parallel jobs to run
my $jobs = 8;

# if true, restrict to PICKed transcripts
my $pick = '';

# min number of HR calls at a position to consider it,
# default 100 works well for us with a heterogeneous cohort of 500
my $min_hr = 100;

# other arguments to $filterBin could be added here if we want to
# occasionally use non-defaults:
## $max_ctrl_hv,$max_ctrl_het,$min_cohort_hv,
## $no_mod, $no_low
## $max_af_gnomad,$max_af_1kg,$max_af_esp

GetOptions ("indir=s" => \$inDir,
	    "outdir=s" => \$outDir,
	    "jobs=i" => \$jobs,
	    "pick" => \$pick,
	    "min_hr=i" => \$min_hr)
    or die("E $0: Error in command line arguments\n");

(-d $inDir) ||
    die "E $0: inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "E $0: cannot opendir inDir $inDir\n";

(-e $outDir) && 
    die "E $0: found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) ||
    die "E $0: cannot mkdir outDir $outDir\n";

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


#############################################

my $pm = new Parallel::ForkManager($jobs);

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

    # using some hard-coded filter params here, add them as options if
    # they need to be customized
    my $com = "perl $filterBin --max_ctrl_hv 3 --max_ctrl_het 10 --no_mod";
    # using defaults for AFs 
    # $com .= " --max_af_gnomad 0.01 --max_af_1kg 0.03 --max_af_esp 0.05"
    $com .= " --min_hr $min_hr";
    ($pick) && ($com .= " --pick");

    if ($gz) {
	$com = "gunzip -c $inDir/$inFile | $com ";
    }
    else {
	$com = "cat $inDir/$inFile | $com ";
    }
    
    $com .= " | perl $reorderBin ";
    $com .= " > $outDir/$outFile";
    # my $now = strftime("%F %T", localtime);
    # warn "I $now: $0 - starting $com\n";
    system($com);
    # $now = strftime("%F %T", localtime);
    # warn "I $now: $0 - Finished $com\n";
    $pm->finish;
}
closedir(INDIR);

$pm->wait_all_children;

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";
