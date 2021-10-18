#!/usr/bin/perl

# 25/06/2021
# NTM

# QC script: examine variant calls on the X and Y chromosomes in the
# samples CSV files.
# Taking into account the Sex column from the samples metadata file, 
# print to stdout some summary statistics of the aberrant calls (any calls on
# the Y in women, HET calls on the X or Y in men), and identify putative
# annotation errors (outliers).


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


my $USAGE = "\nGrab the sex of each patient in the metadata XLSX, and examine the variants 
called on the sex chromosomes in each SAMPLES results CSV files in indir.
Print summary statistics to stdout in CSV format, including info on putative errors (outliers).
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
# hashref, key==sample id, value is the sex
my $sample2sexR;
{
    my @parsed = &parseSamples($samplesFile);
    (@parsed == 5) ||
	die "E $0: need to know the sex of each patient but there's no Sex column in $samplesFile";
    $sample2cohortR = $parsed[0];
    $sample2sexR = $parsed[4];
}

#########################################################
# examine each sample CSV file

# Save relevant metadata and summary stats in @results, one hashref per sample storing:
# patho sample sex nbHetX nbHetY nbHomoY
my @results;

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

    my $patho = $sample2cohortR->{$sample};
    # empty $sample2cohortR as we go, for sanity testing and speed
    delete($sample2cohortR->{$sample});

    my $sex = $sample2sexR->{$sample};

    # number of UNIQUE POSITIONS with HET calls on the X / HET or HV calls on the Y chromosome
    my ($nbHetX,$nbHetY,$nbHomoY) = (0,0,0);
    
    open(INFILE, "$inDir/$inFile") ||
	die "E $0: cannot open inFile $inDir/$inFile\n";

    # indexes of columns of interest
    my ($posCol, $genoCol) = (-1,-1);
    # header: grab column indexes of interest
    my $header = <INFILE>;
    chomp($header);
    my @header = split(/\t/,$header);
    foreach my $i (0..$#header) {
	($header[$i] eq 'POSITION') && ($posCol = $i);
	($header[$i] eq 'GENOTYPE') && ($genoCol = $i);
	($posCol > -1) && ($genoCol > -1) && last;
    }
    (($posCol > -1) && ($genoCol > -1)) || 
	die "E $0: cannot find required column headers in $inFile:\n$header\n";

    # for speed we assume below that POSITION is the first column, check it
    ($posCol == 0) ||
	die "E $0: speed optimization expects POSITION to be the first column but it isn't in $inFile, fix the inFile or the code\n";
    # only count unique positions -> need to remember position of previous line
    my $prevPos = "";
    while (my $line = <INFILE>) {
	# we checked above that POSITION is the first column, immediately skip lines
	# for non-sex chromosomes
	($line =~ /^chr([XY]):/) || next;
	my $chr = $1;

	chomp($line);
	my @line = split(/\t/,$line);
	($line[$posCol] eq $prevPos) && next;
	$prevPos = $line[$posCol];

	if ($chr eq 'X') {
	    if ($line[$genoCol] eq 'HET') {
		$nbHetX++;
	    }
	    # else HOMO on X, skip ie NOOP
	}
	elsif ($chr eq 'Y') {
	    if ($line[$genoCol] eq 'HET') {
		$nbHetY++;
	    }
	    elsif ($line[$genoCol] eq 'HV') {
		$nbHomoY++;
	    }
	    else {
		die "E $0: expecting genotype to be HET or HV, found $line[$genoCol] in $inFile\n";
	    }
	}
	else {
	    die "E $0: expecting chromosome to be X or Y, found $chr in $inFile, impossible!\n";
	}
    }
    
    close(INFILE);
    # patho sample sex nbHetX nbHetY nbHomoY
    push(@results, [$patho,$sample,$sex,$nbHetX,$nbHetY,$nbHomoY]);
}

closedir(INDIR);

# sanity: every sample from metadata should have been seen, assuming
# we didn't analyze a subcohort/incomplete indir
(keys(%$sample2cohortR)) &&
    (print "W $0: ".scalar(keys(%$sample2cohortR))." samples from metadata have no samples.csv files\n");

# calculate mean and stddev of relevant stats, per sex
my ($hetXMeanM, $hetXSdM) = (0,0);
my ($hetXMeanF, $hetXSdF) = (0,0);
my ($hetYMeanM, $hetYSdM) = (0,0);
my ($hetYMeanF, $hetYSdF) = (0,0);
my ($homYMeanM, $homYSdM) = (0,0);
my ($homYMeanF, $homYSdF) = (0,0);

# numbers of individuals of each sex
my ($nbM,$nbF) = (0,0);

# means
foreach my $res (@results) {
    if ($res->[2] eq "M") {
	$nbM++;
	$hetXMeanM += $res->[3];
	$hetYMeanM += $res->[4];
	$homYMeanM += $res->[5];
    }
    elsif ($res->[2] eq "F") {
	$nbF++;
	$hetXMeanF += $res->[3];
	$hetYMeanF += $res->[4];
	$homYMeanF += $res->[5];
    }
    else {
	die "E $0: sex of ".$res->[1]." is neither M or F, it's ".$res->[2]."\n";
    }
}
if ($nbM > 0) {
    $hetXMeanM /= $nbM;
    $hetYMeanM /= $nbM;
    $homYMeanM /= $nbM;
}
if ($nbF > 0) {
    $hetXMeanF /= $nbF;
    $hetYMeanF /= $nbF;
    $homYMeanF /= $nbF;
}

# std devs
foreach my $res (@results) {
    if ($res->[2] eq "M") {
	$hetXSdM += ($res->[3] - $hetXMeanM)**2;
	$hetYSdM += ($res->[4] - $hetYMeanM)**2;
	$homYSdM += ($res->[5] - $homYMeanM)**2;
    }
    else {
	$hetXSdF += ($res->[3] - $hetXMeanF)**2;
	$hetYSdF += ($res->[4] - $hetYMeanF)**2;
	$homYSdF += ($res->[5] - $homYMeanF)**2;
    }
}
if ($nbM > 0) {
    $hetXSdM = sqrt($hetXSdM / $nbM);
    $hetYSdM = sqrt($hetYSdM / $nbM);
    $homYSdM = sqrt($homYSdM / $nbM);
}
if ($nbF > 0) {
    $hetXSdF = sqrt($hetXSdF / $nbF);
    $hetYSdF = sqrt($hetYSdF / $nbF);
    $homYSdF = sqrt($homYSdF / $nbF);
}


# print results as TSV, identifying outliers on-the-fly
# outliers here are defined by abs(X-MEAN) > 3*SD
print "pathologyID\tsampleID\tSex\tnbHetX\tnbHetY\tnbHomoY\tOutlier\n";
foreach my $res (@results) {
    my $outlier = "";
    if ($res->[2] eq "M") {
	if (abs($res->[3] - $hetXMeanM) > 3 * $hetXSdM) {
	    $outlier .= "HetXM:";
	}
	if (abs($res->[4] - $hetYMeanM) > 3 * $hetYSdM) {
	    $outlier .= "HetYM:";
	}
	if (abs($res->[5] - $homYMeanM) > 3 * $homYSdM) {
	    $outlier .= "HomYM:";
	}
    }
    else {
	if (abs($res->[3] - $hetXMeanF) > 3 * $hetXSdF) {
	    $outlier .= "HetXF:";
	}
	if (abs($res->[4] - $hetYMeanF) > 3 * $hetYSdF) {
	    $outlier .= "HetYF:";
	}
	if (abs($res->[5] - $homYMeanF) > 3 * $homYSdF) {
	    $outlier .= "HomYF:";
	}
    }
    # remove trailing ':'
    ($outlier) && (chop($outlier));

    print join("\t", @$res)."\t$outlier\n";
}

# print summary stats at the end
print "\n\n";
print "SUMMARY STATS\n";
print "TYPE MEAN SD\n";
printf("HetXM %.2f %.2f\n", $hetXMeanM, $hetXSdM);
printf("HetYM %.2f %.2f\n", $hetYMeanM, $hetYSdM);
printf("HomYM %.2f %.2f\n", $homYMeanM, $homYSdM);
printf("HetXF %.2f %.2f\n", $hetXMeanF, $hetXSdF);
printf("HetYF %.2f %.2f\n", $hetYMeanF, $hetYSdF);
printf("HomYF %.2f %.2f\n", $homYMeanF, $homYSdF);



$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";
