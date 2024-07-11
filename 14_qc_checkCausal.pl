#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2024
#
# This file is part of grexome-TIMC-Secondary, written by Nicolas Thierry-Mieg
# (CNRS, France) Nicolas.Thierry-Mieg@univ-grenoble-alpes.fr
#
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.
############################################################################################


# 25/06/2021
# NTM

# QC script: for diagnosed patients (ie with a "causal gene" in the samples metadata file),
# look for severe biallelic variants affecting the causal gene's canonical transcript;
# for undiagnosed patients, look for severe biallelic variants affecting a known candidate gene.


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

# samples metadata file
my $samplesFile = "";

# dir containing samples CSV files produced by grexome-TIMC-Secondary,
# default to current dir
my $inDir = ".";

# help: if true just print $USAGE and exit
my $help = '';


my $USAGE = "
Examine the SAMPLES results CSV files in indir, and:
- for each sample that has a causal gene in samplesFile, report if and how the causal gene's
canonical transcript is hit;
- for each sample that DOESN'T have a causal gene in samplesFile, report if the sample has any
severe biallelic variant(s) affecting the canonical transcript of a known candidate gene
(for the sample's pathology, or for another pathology (HV only)).
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

# arrays of strings to print for samples whose "causal gene" is hit by
# biallelic HIGH, biallelic MODHIGH, biallelic MODERATE/LOW, monoallelic,
# or no variant at all
my @causalHigh = ();
my @causalModHigh = ();
my @causalLow = ();
my @mono = ();
my @noVar = ();

# arrays of strings to print for samples without a "causal gene" but where
# the canonical transcript of a candidateGene is hit by HV-HIGH , HV-MODHIGH
# or at least 2 HET-MODHIGH+ variants
my @candidateHVHigh = ();
my @candidateHVModHigh = ();
my @candidateHets= ();

# same as @candidate* above but for candidateGenes of other pathologies, and
# only for HV variants
my @candOtherPathosHVHigh = ();
my @candOtherPathosHVModHigh = ();


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

    open(INFILE, "$inDir/$inFile") ||
	die "E $0: cannot open inFile $inDir/$inFile\n";
    
    # if $sample has a causal gene: see if and how it is hit
    if (defined($sample2causalR->{$sample})) {
        my $causal = $sample2causalR->{$sample};

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
		push(@causalHigh, "$sample\t$patho\t$causal\tHIGH");
	    }
	    elsif ($line[$biallelCol] eq 'MODHIGH') {
		push(@causalModHigh, "$sample\t$patho\t$causal\tMODHIGH");
	    }
	    elsif ($line[$biallelCol] eq 'NO') {
		push(@mono, "$sample\t$patho\t$causal\tMONO-ALLELIC");
	    }
	    else {
		push(@causalLow, "$sample\t$patho\t$causal\tLOW-MODER");
	    }
	    $foundCausal = 1;
	    last;
        }
        ($foundCausal) ||
	    push(@noVar, "$sample\t$patho\t$causal\tNOVARIANT");
    }

    # else $sample doesn't have a causal gene: see if any known candidate gene is severely hit
    else {
        # indexes of columns of interest
        my ($symbolCol,$candCol,$genoCol,$impactCol,$biallelCol,$canonCol) = (-1,-1,-1,-1,-1,-1);
        # header: grab column indexes of interest
        my $header = <INFILE>;
        chomp($header);
        my @header = split(/\t/,$header);
        foreach my $i (0..$#header) {
	    ($header[$i] eq 'SYMBOL') && ($symbolCol = $i);
	    ($header[$i] eq 'KNOWN_CANDIDATE_GENE') && ($candCol = $i);
	    ($header[$i] eq 'GENOTYPE') && ($genoCol = $i);
	    ($header[$i] eq 'IMPACT') && ($impactCol = $i);
	    ($header[$i] eq 'BIALLELIC') && ($biallelCol = $i);
	    ($header[$i] eq 'CANONICAL') && ($canonCol = $i);
	    ($impactCol > -1) && ($candCol > -1) && ($biallelCol > -1) && ($genoCol > -1) &&
		($symbolCol > -1) && ($canonCol > -1) && last;
        }
        (($impactCol > -1) && ($candCol > -1) && ($biallelCol > -1) && ($genoCol > -1) &&
	 ($symbolCol > -1) && ($canonCol > -1)) || 
	    die "E $0: cannot find required column headers for undiagnosed sample in $inFile:\n$header\n";

	# candidate genes (for $sample's patho) that are hit in $sample:
	# key == "$sample\t$patho\t$gene\t$pathosLevels", where $patho is $sample's pathology, and
	# $pathosLevels is a string of comma-separated "$pat:$level" pairs, $gene being a known
	# candidate gene of confidence level $level for pathology $pat,
	# value == ref to an array of 3 ints == numbers of HV_HIGH, HV_MODHIGH
	# and HET_MODHIGH+ variants affecting $gene's canonical transcript
	my %hitCands = ();

	# candidate genes for other pathologies that are hit in $sample:
	# same key as as %hitCands,
	# value is a ref to an array of 2 ints == numbers of HV_HIGH, HV_MODHIGH variants
	# affecting $gene's canonical transcript
	my %hitCandsOtherPathos = ();
        
        while (my $line = <INFILE>) {
	    chomp($line);
	    my @line = split(/\t/,$line);

	    ($line[$canonCol] eq 'YES') || next;
	    ($line[$biallelCol] eq 'HIGH') || ($line[$biallelCol] eq 'MODHIGH') || next;

	    if ($line[$candCol] ne "") {
		# OK we have a severely hit candidate gene, build the key for %hitCands / %hitCandsOtherPathos
		my $key = "$sample\t$patho\t";
		($line[$symbolCol] =~ /^' (\S+)$/) ||
		    die "E $0: cannot extract gene name in:\n$line\n";
		$key .= "$1\t";
		$key .= $line[$candCol];
		# $samePatho: true iff $gene is a candidate for $patho
		my $samePatho = 0;
		if ($line[$candCol] =~ /($patho:[^,]+)/) {
		    # candidate is for this sample's patho
		    $samePatho = 1;
		}

		if ($samePatho) {
		    (defined $hitCands{$key}) || ($hitCands{$key} = [0,0,0]);
		    if (($line[$impactCol] eq 'HIGH') && ($line[$genoCol] eq 'HV')) {
			$hitCands{$key}->[0]++;
		    }
		    elsif (($line[$impactCol] eq 'MODHIGH') && ($line[$genoCol] eq 'HV')) {
			$hitCands{$key}->[1]++;
		    }
		    elsif ((($line[$impactCol] eq 'HIGH') || ($line[$impactCol] eq 'MODHIGH')) &&
			   ($line[$genoCol] eq 'HET')) {
			$hitCands{$key}->[2]++;
		    }
		}
		elsif (($line[$genoCol] eq 'HV') && (($line[$impactCol] eq 'HIGH') ||
						     ($line[$impactCol] eq 'MODHIGH'))) {
		    # only count HVs if gene is a candidate for another patho
		    (defined $hitCandsOtherPathos{$key}) || ($hitCandsOtherPathos{$key} = [0,0]);
		    if ($line[$impactCol] eq 'HIGH') {
			$hitCandsOtherPathos{$key}->[0]++;
		    }
		    elsif ($line[$impactCol] eq 'MODHIGH') {
			$hitCandsOtherPathos{$key}->[1]++;
		    }
		}
	    }
	}

	# done parsing INFILE, now prepare output: for each gene, only print info
	# for the most severe category of variants
	foreach my $key (sort keys(%hitCands)) {
	    if ($hitCands{$key}->[0] > 0) {
		push(@candidateHVHigh, "$key\tHV-HIGH\t".$hitCands{$key}->[0]);
	    }
	    elsif ($hitCands{$key}->[1] > 0) {
		push(@candidateHVModHigh, "$key\tHV-MODHIGH\t".$hitCands{$key}->[1]);
	    }
	    elsif ($hitCands{$key}->[2] > 1) {
		push(@candidateHets, "$key\tHET-MODHIGH+\t".$hitCands{$key}->[2]);
	    }
	    else {
		die "E $0: impossible count vector for key $key: ".@{$hitCands{$key}};
	    }
	}
	foreach my $key (sort keys(%hitCandsOtherPathos)) {
	    if ($hitCandsOtherPathos{$key}->[0] > 0) {
		push(@candOtherPathosHVHigh, "$key\tHV-HIGH\t".$hitCandsOtherPathos{$key}->[0]);
	    }
	    elsif ($hitCandsOtherPathos{$key}->[1] > 0) {
		push(@candOtherPathosHVModHigh, "$key\tHV-MODHIGH\t".$hitCandsOtherPathos{$key}->[1]);
	    }
	    else {
		die "E $0: impossible count vector (otherPathos) for key $key: ".@{$hitCandsOtherPathos{$key}};
	    }
	}
    }
    close(INFILE);
}
closedir(INDIR);


# sanity: every sample from metadata should have been seen, assuming
# we didn't analyze a subcohort/incomplete infile
(keys(%$sample2cohortR)) &&
    (print "W: ".scalar(keys(%$sample2cohortR))." samples from metadata have no samples.csv files\n");

print "Total number of samples in metadata: $nbSamples\n";
print "Examining the samples CSV files in $inDir (canonical transcripts only), we find:\n\n";
print "Samples WITH a causal gene in metadata ($nbCausal), where this causal gene is hit by:\n";
print "a BIALLELIC HIGH variant: ".scalar(@causalHigh)."\n";
print "a BIALLELIC MODHIGH variant: ".scalar(@causalModHigh)."\n";
print "a BIALLELIC MODERATE or LOW variant: ".scalar(@causalLow)."\n";
print "a MONOALLELIC variant: ".scalar(@mono)."\n";
print "NO VARIANT: ".scalar(@noVar)."\n";
print "\n";
print "Samples WITHOUT a causal gene in metadata but where a candidate gene is hit by:\n";
print "a HOMOZYGOUS HIGH variant: ".scalar(@candidateHVHigh)."\n";
print "a HOMOZYGOUS MODHIGH variant: ".scalar(@candidateHVModHigh)."\n";
print "two or more HETEROZYGOUS MODHIGH+ variants: ".scalar(@candidateHets)."\n";
print "a HOMOZYGOUS HIGH variant hitting a candidate gene for another pathology: ".scalar(@candOtherPathosHVHigh)."\n";
print "a HOMOZYGOUS MODHIGH variant hitting a candidate gene for another pathology: ".scalar(@candOtherPathosHVModHigh)."\n";
print"\n";

print "\n#############################################################\n";
if (@causalHigh) {
    print "Samples whose causal gene is BIALLELIC HIGH:\n";
    print join("\n", @causalHigh);
    print "\n\n";
}
if (@causalModHigh) {
    print "Samples whose causal gene is BIALLELIC MODHIGH:\n";
    print join("\n", @causalModHigh);
    print "\n\n";
}
if (@causalLow) {
    print "Samples whose causal gene is BIALLELIC LOW or MODERATE:\n";
    print join("\n", @causalLow);
    print "\n\n";
}
if (@mono) {
    print "Samples whose causal gene is only hit by ONE MONOALLELIC variant:\n";
    print join("\n", @mono);
    print "\n\n";
}
if (@noVar) {
    print "Samples whose causal gene has NO VARIANT:\n";
    print join("\n", @noVar);
    print "\n\n";
}

print "\n#############################################################\n";
if (@candidateHVHigh){
    print "Undiagnosed samples where a candidate gene is hit by N HV-HIGH variant(s):\n";
    print join("\n", @candidateHVHigh);
    print "\n\n";
}

if (@candidateHVModHigh){
    print "Undiagnosed samples where a candidate gene is hit by N HV-MODHIGH variant(s):\n";
    print join("\n", @candidateHVModHigh);
    print "\n\n";
}

if (@candidateHets){
    print "Undiagnosed samples where a candidate gene is hit by N>1 HET MODHIGH+ variants:\n";
    print join("\n", @candidateHets);
    print "\n\n";
}

if (@candOtherPathosHVHigh){
    print "Undiagnosed samples where a candidate gene for another pathology is hit by N HV-HIGH variant(s):\n";
    print join("\n", @candOtherPathosHVHigh);
    print "\n\n";
}

if (@candOtherPathosHVModHigh){
    print "Undiagnosed samples where a candidate gene for another pathology is hit by N HV-MODHIGH variant(s):\n";
    print join("\n", @candOtherPathosHVModHigh);
    print "\n\n";
}

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";
