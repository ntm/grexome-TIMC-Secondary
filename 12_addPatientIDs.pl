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

# Take 3 or 4 arguments: $samplesFile $inDir $outDir [$jobs]
# $samplesFile is the samples metadata xlsx file;
# $inDir must contain cohort TSVs as produced by extractCohorts.pl,
# possibly filtered and reordered with 10_filterAndReorderAll.pl,
# and possibly gzipped;
# $outDir doesn't exist, it will be created and filled with 
# similar TSVs (gzipped if infiles were gzipped), but where every
# $sampleID identifier in the genoData columns becomes "$sampleID($patientID)",
# with $patientID taken from patientID column if it's not empty, 
# specimenID otherwise.
# Filenames get ".patientIDs" added before .csv.
# $jobs is the number of cohort infiles to process in parallel.

use strict;
use warnings;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use POSIX qw(strftime);
use Parallel::ForkManager;

use lib "$RealBin";
use grexome_metaParse qw(parseSamples);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

# number of cohorts to process in parallel
my $jobs = 8;

(@ARGV == 3) || (@ARGV == 4) ||
    die "E $0: needs 3 or 4 args: a samples XLSX, an inDir, a non-existant outDir, and optionally the number of jobs (default=$jobs)\n";
(@ARGV == 4) && ($jobs = $ARGV[3]);
my ($samplesFile, $inDir, $outDir) = @ARGV;
(-d $inDir) ||
    die "E $0: inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "E $0: cannot opendir inDir $inDir\n";
(-e $outDir) && 
    die "E $0: found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "E $0: cannot mkdir outDir $outDir\n";

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


#########################################################
# parse samples file to populate $sample2patientR:
# hashref, key==sampleID, value is patientID || specimenID
my $sample2patientR;
{
    # we want the third hashref returned by parseSamples
    my @parsed = &parseSamples($samplesFile);
    $sample2patientR = $parsed[2];
}

# precompile the regexps we will use to search for each sample,
# in this way RAM usage drops from +50G to 21M !!
# key==sampleID, value is a precompiled regexp 
my %sample2re;
foreach my $sample (keys %$sample2patientR) {
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
