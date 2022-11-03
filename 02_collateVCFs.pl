#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2022
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


# 07/04/2022
# NTM

# Parse 2 VCFs sorted by CHR and POS, print to STDOUT a new VCF with the
# samples from the STDIN VCF, and comprising all variants from both files
# in sorted CHR-POS order.
# See $USAGE.

use strict;
use warnings;
use File::Basename qw(basename);
use Getopt::Long;
use POSIX qw(strftime);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## options / params from the command-line

# full path to the VCF with new variants, possibly (b)gzipped
my $vcf;

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = '
Parse on STDIN a multi-sample VCF, and take as arg a secondary (possibly gzipped) 
multi-sample VCF. Both input VCFs must be sorted by CHR and POS (not checked).
Print to STDOUT a multi-sample VCF comprising all variants from both files, 
sorted by CHR and POS, but respecting the STDIN sample list and order.
Samples absent from STDIN are ignored, and samples missing in the secondary VCF
get ./. as genotype for variants from that file.
This script doesn\'t do anything fancy like normalizing or merging overlapping variants:
the goal is to collate VCFs with different variants (eg SNVs with CNVs) but shared samples.
Arguments (all can be abbreviated to shortest unambiguous prefix):
--vcf : secondary VCF file with new variants
--help : print this USAGE
';

GetOptions ("vcf=s" => \$vcf,
	    "help" => \$help)
    or die("E: $0 - Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n";

($vcf) || die "E $0: you must provide a secondary VCF file. Try $0 --help\n";
(-f $vcf) || die "E $0: the supplied secondary VCF file doesn't exist\n";

($vcf =~ /\.gz$/) && ($vcf = "gunzip -c $vcf | ");
open(VCF, $vcf) ||
    die "E $0: cannot (gunzip-?)open secondary VCF file (as \"$vcf\")\n";

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";

#############################################
# process headers:
# copy all headers from STDIN and ignore those from $vcf,
# but construct a mapping between sampleIDs:
# main2sec[$i] == $j means that data column $i in STDOUT must be
# copied from secondary VCF column $j (for variants from secondary VCF),
# main2sec[$i] == -1 means that sample $i from STDIN is missing from
# the secondary VCF (we will use './.' as GT)
my @main2sec;

{
    # samples in STDIN, indexed by column index
    my @samples;
    # samples in secondary VCF: key==sampleID, value==column index in $vcf
    my %secSamples;
    while(my $line = <STDIN>) {
	if ($line =~ /^##/) {
	    print($line);
	}
	elsif ($line =~ /^#CHROM/) {
	    print($line);
	    chomp($line);
	    @samples = split(/\t/, $line);
	    last;
	}
	else {
	    die "E $0: parsing headers from STDIN, #CHROM line not seen yet but found a non-header line:\n$line";
	}
    }
    while(my $line = <VCF>) {
	if ($line =~ /^##/) {
	    # NOOP
	}
	elsif  ($line =~ /^#CHROM/) {
	    chomp($line);
	    my @secSamps = split(/\t/, $line);
	    foreach my $i (0..$#secSamps) {
		$secSamples{$secSamps[$i]} = $i;
	    }
	    last;
	}
	else {
	    die "E $0: parsing headers from $vcf, #CHROM line not seen yet but found a non-header line:\n$line";
	}
    }
    foreach my $i (0..$#samples) {
	if (defined $secSamples{$samples[$i]}) {
	    $main2sec[$i] = $secSamples{$samples[$i]};
	}
	else {
	    $main2sec[$i] = -1;
	}
    }
}

#############################################
# collate data

# next data lines from each VCF
my ($nextMain, $nextSec);
# CHROM from each next line, removing leading 'chr' if present and replacing
# X Y M/MT by 23-25 for easy sorting (undef if no more data)
my ($nextMainChr, $nextSecChr);
# POS from each next line (undef if no more data)
my ($nextMainPos, $nextSecPos);

# prime with first lines
$nextMain = <STDIN>;
$nextSec = <VCF>;
if ($nextMain && ($nextMain =~ /^([^\t]+)\t(\d+)\t/)) {
    ($nextMainChr,$nextMainPos)=($1,$2);
    $nextMainChr = &chr2num($nextMainChr);
}
elsif ($nextMain) {
    die"E $0: STDIN has a first data line but I can't parse it:\n$nextMain\n";
}
if ($nextSec && ($nextSec =~ /^([^\t]+)\t(\d+)\t/)) {
    ($nextSecChr,$nextSecPos)=($1,$2);
    $nextSecChr = &chr2num($nextSecChr);
}
elsif ($nextSec) {
    die"E $0: secondary VCF has a first data line but I can't parse it:\n$nextSec\n";
}

#######################
# as long as there is data...
while ($nextMain || $nextSec) {
    if (($nextMain) && ((! $nextSec) ||
			($nextMainChr < $nextSecChr) ||
			(($nextMainChr == $nextSecChr) && ($nextMainPos <= $nextSecPos)))) {
	# still have main data and it must be printed now
	print $nextMain;
	# read next line fom main==STDIN
	$nextMain = <STDIN>;
	if ($nextMain && ($nextMain =~ /^([^\t]+)\t(\d+)\t/)) {
	    ($nextMainChr,$nextMainPos)=($1,$2);
	    $nextMainChr = &chr2num($nextMainChr);
	}
	elsif ($nextMain) {
	    die"E $0: STDIN has a data line but I can't parse it:\n$nextMain\n";
	}
	next;
    }
    else {
	# $nextSec line must be processed and printed
	chomp($nextSec);
	my @nextData = split(/\t/, $nextSec);
	my $toPrint = $nextData[0];
	foreach my $i (1..$#main2sec) {
	    if ($main2sec[$i] == -1) {
		$toPrint .= "\t./.";
	    }
	    else {
		$toPrint .= "\t".$nextData[$main2sec[$i]];
	    }
	}
	print $toPrint."\n";
	# read next line fom Sec
	$nextSec = <VCF>;
	if ($nextSec && ($nextSec =~ /^([^\t]+)\t(\d+)\t/)) {
	    ($nextSecChr,$nextSecPos)=($1,$2);
	    $nextSecChr = &chr2num($nextSecChr);
	}
	elsif ($nextSec) {
	    die"E $0: Secondary VCF has a data line but I can't parse it:\n$nextSec\n";
	}
	next;
    }
}


close(VCF);

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";


##########################################################################
### SUBS
##########################################################################

#######################
# take as arg a CHROM string, possibly starting with 'chr' and non-numerical;
# return an int representation of the chrom:
# remove leading 'chr' if present and replace X Y M/MT by 90-92
# [TODO fix this so it's not broken on organisms with other naming conventions]
sub chr2num {
    (@_ == 1) || die "E $0: chr2num needs one arg";
    my ($chr) = @_;
    $chr =~ s/^chr//;
    if ($chr eq 'X') {
	$chr = 90;
    }
    elsif ($chr eq 'Y') {
	$chr = 91;
    }
    elsif (($chr eq 'M') || ($chr eq 'MT')) {
	$chr = 92;
    }
    elsif ($chr !~ /^\d+$/) {
	die "E $0 in chr2num: CHROM is not numeric or X/Y/M/MT, need to fix the code to deal with this";
    }
    return($chr);
}

