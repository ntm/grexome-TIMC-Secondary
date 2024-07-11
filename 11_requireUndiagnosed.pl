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


# 25/03/2020
# NTM

# Takes 2 arguments: $inDir $outDir
# - $inDir must contain cohort TSVs as produced by extractCohorts.pl,
#   possibly filtered and reordered by filterVariants.pl and reorderColumns.pl,
#   possibly with patientIDs;
# - $outDir doesn't exist, it will be created and filled with one TSV
#   per infile named $cohort.final.csv
#
# The outfiles are identical to the infiles except we remove lines where
# COUNT_$cohort_HV == 0 AND COUNT_$cohort_HET == 0 : previous steps
# kept such lines if they had COUNT_*_OTHERCAUSE_* >=1, this was
# necessary for extractSamples (for samples with a casual gene).

use strict;
use warnings;
use File::Basename qw(basename);
use POSIX qw(strftime);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


(@ARGV == 2) || die "E: $0 - needs 2 args: an inDir and a non-existant outDir\n";
my ($inDir, $outDir) = @ARGV;
(-d $inDir) ||
    die "E: $0 - inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "E: $0 - cannot opendir inDir $inDir\n";
(-e $outDir) && 
    die "E: $0 - found argument outDir $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "E: $0 - cannot mkdir outDir $outDir\n";


#########################################################

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";

while (my $inFile = readdir(INDIR)) {
    ($inFile =~ /^\./) && next;
    my $cohort;
    if ($inFile =~ /^(\w+)\..*csv$/) {
	$cohort = $1;
    }
    else {
	warn "W: $0 - cannot parse filename of inFile $inFile, skipping it\n";
	next;
    }

    open(INFILE, "$inDir/$inFile") ||
	die "E: $0 - cannot open infile $inDir/$inFile\n";

    my $outFile = "$cohort.final.csv" ;
    open(OUTFILE, "> $outDir/$outFile") ||
	die "E: $0 - cannot open outfile $outDir/$outFile: $!\n";


    ###################################
    # header line
    my $header = <INFILE>;
    chomp($header);
    my @headers = split(/\t/,$header);

    # columns of COUNT_$cohort_HV and HET
    my ($hvCol,$hetCol) = (-1,-1);

    foreach my $hi (0..$#headers) {
	if ($headers[$hi] eq "COUNT_$cohort"."_HV") {
	    $hvCol = $hi;
	}
	elsif ($headers[$hi] eq "COUNT_$cohort"."_HET") {
	    $hetCol = $hi;
	}
    }
    # sanity check
    (($hvCol >= 0) && ($hetCol >= 0)) || 
	die "E: $0 couldn't find one of HV/HET for $cohort\n";

    print OUTFILE "$header\n";

    # data lines
    while (my $line = <INFILE>) {
	chomp($line);
	my @fields = split(/\t/, $line, -1) ;
	if (($fields[$hvCol] == 0) && ($fields[$hetCol] == 0)) {
	    next;
	}
	else {
	    print OUTFILE "$line\n";
	}
    }
    close(INFILE);
    close(OUTFILE);
}

closedir(INDIR);

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";
