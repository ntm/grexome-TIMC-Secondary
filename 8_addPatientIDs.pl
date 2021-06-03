#!/usr/bin/perl

# 12/08/2019
# NTM

# Take 3 arguments: $samplesFile $inDir $outDir
# $samplesFile is the samples metadata xlsx file;
# $inDir must contain cohort TSVs as produced by extractCohorts.pl,
# possibly filtered and reordered with 7_filterAndReorderAll.pl,
# and possibly gzipped;
# $outDir doesn't exist, it will be created and filled with 
# similar TSVs (gzipped if infiles were gzipped), but where every
# $sampleID identifier in the genoData columns becomes "$sampleID($patientID)",
# with $patientID taken from patientID column if it's not empty, 
# specimenID otherwise.
# Filenames get ".patientIDs" added before .csv.

use strict;
use warnings;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use POSIX qw(strftime);

use lib "$RealBin";
use grexome_metaParse qw(parseSamples);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


(@ARGV == 3) || die "E $0: needs 3 args: a samples XLSX, an inDir and a non-existant outDir\n";
my ($samplesFile, $inDir, $outDir) = @ARGV;
(-d $inDir) ||
    die "E $0: inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "E $0: cannot opendir inDir $inDir\n";
(-e $outDir) && 
    die "E $0: found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "E $0: cannot mkdir outDir $outDir\n";

my $now = strftime("%F %T", localtime);
warn "I $0: $now - starting to run\n";


#########################################################
# parse samples file to populate $sample2patientR:
# hashref, key==sampleID, value is patientID || specimenID
my $sample2patientR;
{
    # we want the third hashref returned by parseSamples
    my @parsed = &parseSamples($samplesFile);
    $sample2patientR = $parsed[2];
}

#########################################################
# read infiles

while (my $inFile = readdir(INDIR)) {
    ($inFile =~ /^\./) && next;
    my ($fileStart,$gz);
    if ($inFile =~ /^(.+)\.csv$/) {
	$fileStart = $1;
    }
    elsif ($inFile =~ /^(.+)\.csv\.gz$/) {
	$fileStart = $1;
	$gz = 1;
    }
    else {
	warn "W $0: cannot parse filename of inFile $inDir/$inFile, skipping it\n";
    }

    my $inFull = "$inDir/$inFile";
    ($gz) && ($inFull = "gunzip -c $inFull | ");
    open(IN, $inFull) ||
	die "E $0: cannot (gunzip-?)open cohort datafile $inDir/$inFile (as $inFull)\n";

    my $outFile = "$outDir/$fileStart.patientIDs.csv";
    my $outFull = " > $outFile";
    ($gz) && ($outFull = " | gzip -c $outFull.gz");
    open (OUT, $outFull) || 
	die "E $0: cannot (gzip-?)open $outFile for writing (as $outFull)\n";

    while (my $line = <IN>) {
	chomp($line);
        # add trailing ',' so we know sampleID is always followed by some char
        $line .= ',';
	# chuck norris style: brutal but it works...
	foreach my $sample (keys %$sample2patientR) {
	    $line =~ s/$sample([\[,\s|])/$sample($sample2patientR->{$sample})$1/g ;
	}
	# remove the trailing ,
	($line =~ s/,$//) || die "E $0: cannot remove trailing , in:\n$line\n";
	print OUT "$line\n";
    }
    close(IN);
    close(OUT);
}
closedir(INDIR);

$now = strftime("%F %T", localtime);
warn "I $0: $now - ALL DONE, completed successfully!\n";

