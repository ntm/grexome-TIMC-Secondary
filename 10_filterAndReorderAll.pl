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


# 26/03/2018
# NTM

# Every infile from $inDir is filtered and columns are optionally reordered,
# results go into $outDir (which must not pre-exist)
# - $inDir must contain cohort or sample TSVs (possibly gzipped) as 
#   produced by extractCohorts.pl or extractSamples.pl;
# - $outDir doesn't exist, it will be created and filled with one TSV
#   per infile (never gzipped), adding .canon to the filename if called
#   with --canonical.
# This script can only apply filters on COUNTs and --canonical, other filters
# from $filterBin should have been applied earlier (when we still had a single 
# file for all cohorts).

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
my ($filterBin,$reorderBin) = ("$RealBin/07_filterVariants.pl","$RealBin/10_reorderColumns.pl");


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
# should be the same as what was passed to 08_addGTEX.pl
my $favoriteTissues = "testis,ovary";


# arguments for filtering: no default values, all filters disabled by default
my $max_ctrl_hv; # COUNT_NEGCTRL_HV <= $x
my $max_ctrl_het; # COUNT_NEGCTRL_HET <= $x
my $min_cohort_hv; # COUNT_$cohort_HV >= $x
my $min_hr; # COUNT_HR >= $x
my $canon = ''; # if enabled, only keep lines with CANONICAL==YES

GetOptions ("indir=s" => \$inDir,
            "outdir=s" => \$outDir,
            "jobs=i" => \$jobs,
            "reorder" => \$reorder,
            "favoriteTissues=s" => \$favoriteTissues,
            "max_ctrl_hv=i" => \$max_ctrl_hv,
            "max_ctrl_het=i" => \$max_ctrl_het,
            "min_cohort_hv=i" => \$min_cohort_hv,
            "min_hr=i" => \$min_hr,
            "canonical" => \$canon)
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
$pm->run_on_finish( sub { ($_[1]) && ($childFailed=1) });

# verbose only for the first inFile
my $verbose = 1;

while (my $inFile = readdir(INDIR)) {
    ($inFile =~ /^\./) && next;
    if ($childFailed) {
        $now = strftime("%F %T", localtime);
        die "E $now: $0 FAILED - some child died, no point going on\n";
    }
    if ($pm->start) {
        $verbose = 0;
        next;
    }
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
    ($canon) && ($com .= " --canonical");
    
    if ($gz) {
        $com = "gunzip -c $inDir/$inFile | $com ";
    }
    else {
        $com = "cat $inDir/$inFile | $com ";
    }

    # reorder columns if requested
    ($reorder) && ($com .= " | perl $reorderBin $verbose $favoriteTissues ");

    $com .= " > $outDir/$outFile";

    # fail if any component of the pipe fails
    $com =~ s/"/\\"/g;
    $com = "bash -o pipefail -c \" $com \"";
    system($com) && die "E $0: filter and/or reorder failed for $inFile\n";
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
