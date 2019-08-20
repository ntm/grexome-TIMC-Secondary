#!/usr/bin/perl


# 26/03/2018
# NTM

# Takes 3 arguments: $pick $inDir $outDir
# - $pick is 0 or 1, if true we apply --pick as a filter;
# - $inDir must contain cohort or sample TSVs (possibly gzipped) as 
#   produced by extractCohorts.pl or extractSamples.pl;
# - $outDir doesn't exist, it will be created and filled with one TSV
#   per infile (never gzipped), adding .filtered or .filtered.pick to 
#   the filename.
#

use strict;
use warnings;
use POSIX qw(strftime);
use Parallel::ForkManager;

# number of parallel jobs to run
my $numJobs = 8;

# full path to finalFilter.pl, unfortunately hard-coded for now
# but with several possibbilities
my $filterBin;
my @possibleFilterBins = ("/home/nthierry/PierreRay/Grexome/SecondaryAnalyses/8_finalFilters.pl",
			  "/home/nthierry/VariantCalling/GrexomeFauve/SecondaryAnalyses/8_finalFilters.pl");
foreach my $f (@possibleFilterBins) {
    (-f $f) && ($filterBin = "perl $f") && last;
}
($filterBin) || 
    die "Sorry, can't find 8_finalFilters.pl, update \@possibleFilterBins\n";

(@ARGV == 3) || die "needs 3 args: PICK (boolean), an inDir and a non-existant outDir\n";
my ($pick, $inDir, $outDir) = @ARGV;
($pick == 0) || ($pick == 1) || die "first arg must be 0 or 1\n";
(-d $inDir) ||
    die "inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "cannot opendir inDir $inDir\n";
(-e $outDir) && 
    die "found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "cannot mkdir outDir $outDir\n";


# create fork manager
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
    
    my $com = "$filterBin --max_ctrl_hv 3 --max_ctrl_het 10 --no_mod";
    # using defaults for AFs 
    # $com .= " --max_af_gnomad 0.01 --max_af_1kg 0.03 --max_af_esp 0.05"
    ($pick) && ($com .= " --pick");
    if ($gz) {
	$com = "gunzip -c $inDir/$inFile | $com ";
    }
    else {
	$com = "cat $inDir/$inFile | $com ";
    }
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
