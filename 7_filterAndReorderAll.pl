#!/usr/bin/perl


# 26/03/2018
# NTM

# Every infile from $inDir is filtered and columns are optionally reordered,
# results go into $outDir (which must not pre-exist)
# - $inDir must contain cohort or sample TSVs (possibly gzipped) as 
#   produced by extractCohorts.pl or extractSamples.pl;
# - $outDir doesn't exist, it will be created and filled with one TSV
#   per infile (never gzipped), adding .canon to the filename if called
#   with --canonical.

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

# $reorder: if enabled we reorder columns with $reorderBin
my $reorder = '';

# for $reorderBin: comma-separated list of favorite tissues (with default), 
# should be the same as what was passed to 5_addGTEX.pl
my $favoriteTissues = "testis,ovary";


# arguments for filtering: no default values, all filters disabled by default
my $max_ctrl_hv; # COUNT_NEGCTRL_HV <= $x
my $max_ctrl_het; # COUNT_NEGCTRL_HET <= $x
my $min_cohort_hv; # COUNT_$cohort_HV >= $x
my $min_hr; # COUNT_HR >= $x
my $no_mod = ''; # if enabled, filter out MODIFIER impacts
my $no_low = ''; # if enabled, filter out LOW impacts
my $canon = ''; # if enabled, only keep lines with CANONICAL==YES
my $max_af_gnomad; # gnomADe_AF <= $x AND gnomADg_AF <= $x
my $max_af_1kg; # AF <= $x, this is 1KG phase 3

GetOptions ("indir=s" => \$inDir,
	    "outdir=s" => \$outDir,
	    "jobs=i" => \$jobs,
	    "reorder" => \$reorder,
	    "favoriteTissues=s" => \$favoriteTissues,
	    "max_ctrl_hv=i" => \$max_ctrl_hv,
	    "max_ctrl_het=i" => \$max_ctrl_het,
	    "min_cohort_hv=i" => \$min_cohort_hv,
	    "min_hr=i" => \$min_hr,
	    "no_mod" => \$no_mod,
	    "no_low" => \$no_low,
	    "canonical" => \$canon,
	    "max_af_gnomad=f" => \$max_af_gnomad,
	    "max_af_1kg=f" => \$max_af_1kg)
    or die("E $0: Error in command line arguments\n");

($inDir) ||
    die "E $0: you must provide an indir containing files to filter and reorder\n";
(-d $inDir) ||
    die "E $0: inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "E $0: cannot opendir inDir $inDir\n";

($outDir) ||
    die "E $0: you must provide a non-existing outdir\n";
(-e $outDir) && 
    die "E $0: found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) ||
    die "E $0: cannot mkdir outDir $outDir\n";

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


#############################################

my $pm = new Parallel::ForkManager($jobs);

# $childFailed will become non-zero if at least one child died
my $childFailed = 0;
# Set up a callback so the parent knows if a child dies
$pm->run_on_finish( sub {
    my ($pid, $exit_code) = @_;
    ($exit_code) && ($childFailed=1);
		    });

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
    ($canon) && ($outFile .= ".canon");
    $outFile .= ".csv";

    # filtering options
    my $com = "perl $filterBin";
    ($max_ctrl_hv) && ($com .= " --max_ctrl_hv=$max_ctrl_hv");
    ($max_ctrl_het) && ($com .= " --max_ctrl_het=$max_ctrl_het");
    ($min_cohort_hv) && ($com .= " --min_cohort_hv=$min_cohort_hv");
    ($min_hr) && ($com .= " --min_hr=$min_hr");
    ($no_mod) && ($com .= " --no_mod");
    ($no_low) && ($com .= " --no_low");
    ($canon) && ($com .= " --canonical");
    ($max_af_gnomad) && ($com .= " --max_af_gnomad=$max_af_gnomad");
    ($max_af_1kg) && ($com .= " --max_af_1kg=$max_af_1kg");
    
    if ($gz) {
	$com = "gunzip -c $inDir/$inFile | $com ";
    }
    else {
	$com = "cat $inDir/$inFile | $com ";
    }

    # reorder columns if requested
    ($reorder) && ($com .= " | perl $reorderBin $favoriteTissues ");

    $com .= " > $outDir/$outFile";
    # my $now = strftime("%F %T", localtime);
    # warn "I $now: $0 - starting $com\n";
    system($com) && die "E $0: filter and/or reorder failed for $inFile\n";
    # $now = strftime("%F %T", localtime);
    # warn "I $now: $0 - Finished $com\n";
    $pm->finish;
}
closedir(INDIR);

$pm->wait_all_children;

$now = strftime("%F %T", localtime);
if ($childFailed) {
    die "E $now: $0 FAILED\n";
}
else {
    warn "I $now: $0 - ALL DONE, completed successfully!\n";
}
