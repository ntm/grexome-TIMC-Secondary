#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2025
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


# 12/08/2019
# NTM

# Append ($patientID) to each $sampleID (present in --SOIs if provided),
# for each file in --inDir.

use strict;
use warnings;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use Getopt::Long;
use POSIX qw(strftime);
use Parallel::ForkManager;

use lib "$RealBin";
use grexome_metaParse qw(parseSamples);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## options / params from the command-line

# samples metadata xlsx file:
my $samplesFile;

# $inDir can contain Cohorts/Transcripts/Samples TSVs, possibly filtered and
# reordered with 10_filterAndReorderAll.pl, and possibly gzipped;
my $inDir = "";

# $outDir must not pre-exist; it will be created and filled with one file for
# each inFile, similar name with ".patientIDs" prepended before .csv, gzipped
# if infiles were gzipped
my $outDir = "";

# optional txt file listing sampleIDs of interest, one per line. If provided,
# only these sampleIDs will be enriched with their patientIDs.
my $SOIsFile;
    
# number of files to process in parallel
my $jobs = 8;

GetOptions ("samples=s" => \$samplesFile,
            "inDir=s" => \$inDir,
            "outDir=s" => \$outDir,
            "SOIs=s" => \$SOIsFile,
            "jobs=i" => \$jobs)
    or die("E: $0 - Error in command line arguments");

# make sure required options were provided and sanity check them
(($samplesFile) && (-f $samplesFile)) ||
    die "E $0: need a samples file with --samples\n";

(($inDir) && (-d $inDir)) ||
    die "E $0: --inDir required and $inDir must be a directory\n";
opendir(INDIR, $inDir) || die "E $0: cannot opendir inDir $inDir\n";

(($outDir) && (!-e $outDir)) ||
    die "E $0: --outDir must be provided and $outDir must not pre-exist\n";
mkdir($outDir) || die "E $0: cannot mkdir outDir $outDir\n";

(!$SOIsFile) || (-f $SOIsFile) ||
    die "E $0: --SOIs is optional but if provided it must exist\n";

($jobs >= 1) || die "E $0: optional --jobs must be an integer >= 1\n";

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


#########################################################
# parse samples file to populate $sample2patientR:
# hashref, key==sampleID, value is patientID || specimenID as returned
# by parseSamples()
my $sample2patientR;
{
    # we want the third hashref returned by parseSamples
    my @parsed = &parseSamples($samplesFile);
    $sample2patientR = $parsed[2];
}

# samples of interest, if provided
my %SOIs;
if ($SOIsFile) {
    open(my $SOIsFH, $SOIsFile) || die "E $0: cannot open --SOIs file $SOIsFile: $!\n";
    while (my $soi = <$SOIsFH>) {
        chomp($soi);
        $SOIs{$soi} = 1;
    }
    close($SOIsFH);
}

# precompile the regexps we will use to search for each sample,
# in this way RAM usage drops from +50G to 21M !!
# key==sampleID, value is a precompiled regexp 
my %sample2re;
foreach my $sample (keys %$sample2patientR) {
    ($SOIsFile) && (! $SOIs{$sample}) && next;
    $sample2re{$sample} = qr/$sample([\[,\s|])/;
}


#########################################################
# read infiles

my $pm = new Parallel::ForkManager($jobs);

# $childFailed will become non-zero if at least one child died
my $childFailed = 0;
# Set up a callback so the parent knows if a child dies
$pm->run_on_finish( sub { ($_[1]) && ($childFailed=1) });

while (my $inFile = readdir(INDIR)) {
    ($inFile =~ /^\./) && next;
    if ($childFailed) {
        $now = strftime("%F %T", localtime);
        die "E $now: $0 FAILED - some child died, no point going on\n";
    }
    $pm->start && next;
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
        $pm->finish;
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
            $line =~ s/$sample2re{$sample}/$sample($sample2patientR->{$sample})$1/g ;
        }
        # remove the trailing ,
        ($line =~ s/,$//) || die "E $0: cannot remove trailing , in:\n$line\n";
        print OUT "$line\n";
    }
    close(IN);
    close(OUT);
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
