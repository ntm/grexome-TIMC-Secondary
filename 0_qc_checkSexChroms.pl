#!/usr/bin/perl

# 26/10/2021
# NTM

# QC script: examine variant calls on the X and Y chromosomes in
# single-sample FILTERED GVCFs; taking into account the sex of
# each sample, print to stdout some summary statistics of the
# calls and identify putative annotation errors, ie outliers in
# the numbers of aberrant calls (any calls on the Y in women,
# HET calls on the X or Y in men).
# You can provide the output of this script from a previous run
# with --prevQC, to avoid parsing every GVCF whenever you get
# a few new samples (n+1 scenario).

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

# samples XLSX file
my $samplesFile = "";

# dir containing the single-sample filtered GVCFs to analyze
my $inDir = "";

# previous QC file (if any), produced by this script for some
# of the same samples. If provided, the counts of HET and HV calls
# on X and Y are taken from $prevQC for samples already in that file.
my $prevQC = "";

# tabix binary, with path if needed
my $tabix = "tabix";

# force: if true any pre-existing lines from $prevQC that seem obsolete
# compared to $samplesFile (nonexistant sampleID, or pathologyID or sex
# mismatch) are discarded and rebuilt from scratch;
# without --force we die if discrepencies are detected
my $force = '';

# help: if true just print $USAGE and exit
my $help = '';


my $USAGE = "\nGrab the sex of each patient in the metadata XLSX, and count the variants 
called on the sex chromosomes in each patient's FILTERED GVCF in indir.
Print summary statistics to stdout in TSV format, including info on putative errors (outliers).
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--samplesFile : samples metadata xlsx file, with path
--indir : dir containing single-sample filtered GVCF files
--prevQC : previous QC file (if any) produced by this script for some of the same samples,
	   existing counts will be re-used instead of parsing every GVCF again
--tabix [$tabix] : name of tabix binary, with path if it's not in PATH
--force : discard lines from prevQC if they disagree with samplesFile (eg pathologyID mismatch)
--help : print this USAGE";

GetOptions ("samplesFile=s" => \$samplesFile,
	    "indir=s" => \$inDir,
	    "prevQC=s" => \$prevQC,
	    "tabix=s" => \$tabix,
	    "force" => \$force,
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

#$prevQC is optional but if provided the file must exist
(! $prevQC) || (-f $prevQC) ||
    die "E $0: --prevQC specified but $prevQC isn't a file / doesn't exist\n";

(`which $tabix` =~ /$tabix/) ||
    die "E $0: the tabix executable $tabix can't be found\n";


my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";

#########################################################
# We will save relevant metadata and summary stats in %results, 
# key is $sample, value is an arrayref storing: 
# [patho sample sex nbHetX nbHomoX nbHetY nbHomoY]
my %results;

# same structure, populated from prevQC
my %resultsPrev;

#########################################################

# parse useful info from samples metadata file:
# hashref, key==sampleID, value==pathologyID
my $sample2cohortR;
# hashref, key==sampleID, value is the sex
my $sample2sexR;
{
    my @parsed = &parseSamples($samplesFile);
    (@parsed == 5) ||
	die "E $0: need to know the sex of each patient but there's no Sex column in $samplesFile";
    $sample2cohortR = $parsed[0];
    $sample2sexR = $parsed[4];
}

# parse counts from $prevQC if provided, save in %resultsPrev
if ($prevQC) {
    open(PREV, $prevQC) ||
	die "E $0: prevQC specified but cannot open file $prevQC: $!\n";

    # $ok: if false after parsing prevQC, there was at least one
    # metadata mismatch between prevQC and samplesFile
    my $ok=1;

    <PREV>; # skip header

    while (my $line = <PREV>) {
	chomp($line);
	# file has per-sample counts, then a blank line, then global summary stats;
	# stop at the blank line:
	($line) || last;

	# -1 to grab "outlier" column even if it's empty
	my @fields = split(/\t/, $line, -1);
	(@fields == 8) ||
	    die "E $0: line from prevQC doesn't have correct number of fields: $line\n";
	# last field "Outlier" is always rebuilt, discard previous value
	pop(@fields);

	my $sample = $fields[1];
	# check that metadata matches $samplesFile
	my $thisOK = 1;
	if ($fields[0] ne $sample2cohortR->{$sample}) {
	    warn "W $0: pathologyID for $sample in prevQC differs from samplesFile $samplesFile\n";
	    $thisOK = 0;
	}
	if ($fields[2] ne $sample2sexR->{$sample}) {
	    warn "W $0: sex for $sample in prevQC differs from samplesFile $samplesFile\n";
	    $thisOK = 0;
	}

	if ($thisOK) {
	    ($resultsPrev{$sample}) &&
		die "E $0: found 2 lines in prevQC $prevQC with the same sampleID $sample\n";
	    $resultsPrev{$sample} = \@fields;
	}
	else {
	    $ok = 0;
	}
    }

    close(PREV);
    
    if (!$ok) {
	if ($force) {
	    warn "W $0: metadata mismatch between prevQC and samplesFile, --force mode => stale prevQC lines will be discarded\n";
	}
	else {
	    die "E $0:  metadata mismatch between prevQC and samplesFile, dying. Use --force to discard stale prevQC lines\n";
	}
    }
}

#########################################################
# examine each sample VCF file

foreach my $inFile (sort(readdir(INDIR))) {
    # silently skip files that don't look like bgzipped (g)vcf
    ($inFile =~ /vcf\.gz$/) || next;
    # this script requires tabix-indexed VCFs
    (-e "$inDir/$inFile.tbi") ||
	((warn "W $0: can't find tabix index $inDir/$inFile.tbi , skipping (G)VCF $inFile\n") && next);

    my ($sample,$patho,$sex);

    # parse sampleID from VCF header
    my $header = `$tabix -H $inDir/$inFile | tail -n 1`;
    chomp($header);
    ($header =~ /^#CHROM\t/) ||
	die "E $0: last header line of gvcf $inDir/$inFile doesn't start with CHROM?\n$header\n";
    my @headers = split(/\t/,$header);
    (@headers == 10) ||
	die "E $0: gvcf $inDir/$inFile is not single-sample according to header:\n$header\n";
    $sample = $headers[9];

    if ($sample2cohortR->{$sample}) {
	$patho = $sample2cohortR->{$sample};
	$sex = $sample2sexR->{$sample};
	# empty $sample2cohortR as we go, for sanity testing
	delete($sample2cohortR->{$sample});
    }
    else {
	warn "W $0: inFile $inFile has data for sample $sample which isn't in $samplesFile, skipping it\n";
	next;
    }

    if ($resultsPrev{$sample}) {
	# we already have counts from prevQC
	$results{$sample} = $resultsPrev{$sample};
	delete($resultsPrev{$sample});
	next;
    }
    else {
	# otherwise parse $inFile to count HET/HV calls on X/Y
	$now = strftime("%F %T", localtime);
	warn "I $now: $0 - starting to parse $inDir/$inFile\n";
	my @res = ($patho,$sample,$sex, &countCalls("$inDir/$inFile",$tabix));
	$results{$sample} = \@res;
    }
}

# sanity
if (my @prevMissing = sort(keys %resultsPrev)) {
    if ($force) {
	warn "W $0: prevQC has data for samples that don't have files in inDir,".
	    " --force mode => these prevQC lines will be discarded. Samples:\n".join(" ",@prevMissing)."\n";
    }
    else {
	die "E $0: prevQC has data for samples that don't have files in inDir,".
	    " dying. Use --force to discard stale prevQC lines. Samples:\n".join(" ",@prevMissing)."\n";
    }
}
# sanity: every sample from metadata should have been seen, assuming
# we didn't analyze a subcohort/incomplete indir
if (my @metaMissing = sort(keys %$sample2cohortR)) {
    warn "W $0: some samples from samplesFile have no GVCF in $inDir, ignoring them:\n".
	join(" ",@metaMissing)."\n";
}
closedir(INDIR);


# calculate mean and stddev of relevant stats, per sex
my ($hetXMeanM, $hetXSdM) = (0,0);
my ($hetXMeanF, $hetXSdF) = (0,0);
my ($homXMeanM, $homXSdM) = (0,0);
my ($homXMeanF, $homXSdF) = (0,0);
my ($hetYMeanM, $hetYSdM) = (0,0);
my ($hetYMeanF, $hetYSdF) = (0,0);
my ($homYMeanM, $homYSdM) = (0,0);
my ($homYMeanF, $homYSdF) = (0,0);

# numbers of individuals of each sex
my ($nbM,$nbF) = (0,0);

# means
foreach my $res (values(%results)) {
    if ($res->[2] eq "M") {
	$nbM++;
	$hetXMeanM += $res->[3];
	$homXMeanM += $res->[4];
	$hetYMeanM += $res->[5];
	$homYMeanM += $res->[6];
    }
    elsif ($res->[2] eq "F") {
	$nbF++;
	$hetXMeanF += $res->[3];
	$homXMeanF += $res->[4];
	$hetYMeanF += $res->[5];
	$homYMeanF += $res->[6];
    }
    else {
	die "E $0: sex of ".$res->[1]." is neither M or F, it's ".$res->[2]."\n";
    }
}
if ($nbM > 0) {
    $hetXMeanM /= $nbM;
    $homXMeanM /= $nbM;
    $hetYMeanM /= $nbM;
    $homYMeanM /= $nbM;
}
if ($nbF > 0) {
    $hetXMeanF /= $nbF;
    $homXMeanF /= $nbF;
    $hetYMeanF /= $nbF;
    $homYMeanF /= $nbF;
}

# std devs
foreach my $res (values(%results)) {
    if ($res->[2] eq "M") {
	$hetXSdM += ($res->[3] - $hetXMeanM)**2;
	$homXSdM += ($res->[4] - $homXMeanM)**2;
	$hetYSdM += ($res->[5] - $hetYMeanM)**2;
	$homYSdM += ($res->[6] - $homYMeanM)**2;
    }
    else {
	$hetXSdF += ($res->[3] - $hetXMeanF)**2;
	$homXSdF += ($res->[4] - $homXMeanF)**2;
	$hetYSdF += ($res->[5] - $hetYMeanF)**2;
	$homYSdF += ($res->[6] - $homYMeanF)**2;
    }
}
if ($nbM > 0) {
    $hetXSdM = sqrt($hetXSdM / $nbM);
    $homXSdM = sqrt($homXSdM / $nbM);
    $hetYSdM = sqrt($hetYSdM / $nbM);
    $homYSdM = sqrt($homYSdM / $nbM);
}
if ($nbF > 0) {
    $hetXSdF = sqrt($hetXSdF / $nbF);
    $homXSdF = sqrt($homXSdF / $nbF);
    $hetYSdF = sqrt($hetYSdF / $nbF);
    $homYSdF = sqrt($homYSdF / $nbF);
}


# print results as TSV, identifying outliers on-the-fly
# outliers here are defined by abs(X-MEAN) > 3*SD
print "pathologyID\tsampleID\tSex\tnbHetX\tnbHomoX\tnbHetY\tnbHomoY\tOutlier\n";
foreach my $res (values(%results)) {
    my $outlier = "";
    if ($res->[2] eq "M") {
	if (abs($res->[3] - $hetXMeanM) > 3 * $hetXSdM) {
	    $outlier .= "HetXM:";
	}
	if (abs($res->[4] - $homXMeanM) > 3 * $homXSdM) {
	    $outlier .= "HomXM:";
	}
	if (abs($res->[5] - $hetYMeanM) > 3 * $hetYSdM) {
	    $outlier .= "HetYM:";
	}
	if (abs($res->[6] - $homYMeanM) > 3 * $homYSdM) {
	    $outlier .= "HomYM:";
	}
    }
    else {
	if (abs($res->[3] - $hetXMeanF) > 3 * $hetXSdF) {
	    $outlier .= "HetXF:";
	}
	if (abs($res->[4] - $homXMeanF) > 3 * $homXSdF) {
	    $outlier .= "HomXF:";
	}
	if (abs($res->[5] - $hetYMeanF) > 3 * $hetYSdF) {
	    $outlier .= "HetYF:";
	}
	if (abs($res->[6] - $homYMeanF) > 3 * $homYSdF) {
	    $outlier .= "HomYF:";
	}
    }
    # remove trailing ':'
    ($outlier) && (chop($outlier));

    print join("\t", @$res)."\t$outlier\n";
}

# print summary stats at the end after blank line
print "\n\n";
print "SUMMARY STATS\n";
print "TYPE MEAN SD\n";
printf("HetXM %.2f %.2f\n", $hetXMeanM, $hetXSdM);
printf("HomXM %.2f %.2f\n", $homXMeanM, $homXSdM);
printf("HetYM %.2f %.2f\n", $hetYMeanM, $hetYSdM);
printf("HomYM %.2f %.2f\n", $homYMeanM, $homYSdM);
printf("HetXF %.2f %.2f\n", $hetXMeanF, $hetXSdF);
printf("HomXF %.2f %.2f\n", $homXMeanF, $homXSdF);
printf("HetYF %.2f %.2f\n", $hetYMeanF, $hetYSdF);
printf("HomYF %.2f %.2f\n", $homYMeanF, $homYSdF);



$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";



#########################################################
# subs

# Process one single-sample filtered (by filterBadCalls.pl), bgzipped and
# tabix-indexed (G)VCF, containing at least all non-HR calls on chrX and chrY
# (any HR calls and non X|Y chroms are skipped).
# Return the number of HET or HV variant calls on the X and Y chromosomes,
# as a list: ($nbHetX,$nbHomoX,$nbHetY,$nbHomoY)
# Die on errors.
# Pre-conditions: $tabix exists, and $gvcf is single-sample and tabix-indexed
sub countCalls {
    (@_ == 2) || die "E $0: countCalls needs 2 args\n";
    my ($gvcf,$tabix) = @_;

    (-e $gvcf) || die "E $0: countCalls called with GVCF that doesn't exist: $gvcf\n";
    (-e "$gvcf.tbi") || die "E $0: countCalls called with GVCF lacking a tbi: $gvcf.tbi\n";
    (`which $tabix` =~ /$tabix/) ||
	die "E $0: countCalls called with bad tabix binary: $tabix\n";

    # number of HET or HV calls on the X or Y chromosomes
    my ($nbHetX,$nbHomoX,$nbHetY,$nbHomoY) = (0,0,0,0);

    open(my $GVCF, "$tabix $gvcf chrX chrY |") ||
	die "E $0: countCalls cannot tabix-open GVCF $gvcf\n";

    while (my $line = <$GVCF>) {
	chomp($line);
	# grab chrom
	($line =~ /^chr([XY])\t/) ||
	    die "E $0: countCalls cannot grab chrom from $gvcf in line:\n$line\n";
	my $chr = $1;
	# grab geno: single-sample => last column
	($line =~ /\t(\d+)\/(\d+):[^\t]+$/) ||
	    die "E $0: cannot grab genotype in:\n$line\n";
	my ($gt1,$gt2) = ($1,$2);
	if (($gt1 == $gt2) && ($gt1 != 0)) {
	    # HV, non-HR
	    if ($chr eq "X") { $nbHomoX++; }
	    elsif ($chr eq "Y") { $nbHomoY++; }
	}
	elsif ($gt1 != $gt2) {
	    # HET
	    if ($chr eq "X") { $nbHetX++; }
	    elsif ($chr eq "Y") { $nbHetY++; }
	}
	# else: HR => NOOP
    }

    close($GVCF);
    return($nbHetX,$nbHomoX,$nbHetY,$nbHomoY);
}

