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


# 19/10/2022
# NTM

# Read on stdin a TSV file as produced by vcf2tsv.pl, print the exact same file to stdout
# but along the way check whether each candidateGene or CausalGene is seen, reporting
# missing genes on stderr.
# See $USAGE.

use strict;
use warnings;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use Getopt::Long;
use POSIX qw(strftime);

use lib "$RealBin";
use grexome_metaParse qw(parseCandidateGenes);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## options / params from the command-line

# samples and candidateGene XSLX files
my ($samplesFile, $candidatesFiles) = ("","");

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "\nParse on STDIN a TSV file as produced by vcf2tsv.pl, and take as arg samples and candidateGenes xlsx files;
print to STDOUT exactly the same TSV, but also check whether each candidate or causal gene is seen in the TSV's SYMBOL column, 
reporting any missing genes (likely typoes) to stderr.\n
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--samples string : samples metadata xlsx file, with path
--candidateGenes string [optional, if not provided we still search for causal genes from samples file] : comma-separated list of candidateGene xlsx files, with paths
--help : print this USAGE";

GetOptions ("samples=s" => \$samplesFile,
	    "candidateGenes=s" => \$candidatesFiles,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($samplesFile) || die "E $0: you must provide a samples file\n";
(-f $samplesFile) || die "E $0: the supplied samples file doesn't exist\n";


my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


#########################################################
# parse the provided metadata files

# %knownCandidateGenes: key==$gene, value is a comma-separated list
# of "$patho:$score" pairs.
# This hash is initially filled from the XLSX files, and will then be emptied as
# we parse stdin. At the end any remaining genes in the hash were never seen.
my %knownCandidateGenes;

{
    # $parseCandidateGenes actually returns a hashref, key==$cohort, value is
    # a hashref whose keys are gene names and values are the "Confidence scores"
    my $knownCandsR = &parseCandidateGenes($candidatesFiles, $samplesFile);

    foreach my $c (keys(%$knownCandsR)) {
	foreach my $gene (keys(%{$knownCandsR->{$c}})) {
	    my $score = $knownCandsR->{$c}->{$gene};
	    if ($knownCandidateGenes{$gene}) {
		$knownCandidateGenes{$gene} .= ",$c:$score";
	    }
	    else {
		$knownCandidateGenes{$gene} = "$c:$score";
	    }
	}
    }
}

#############################################
# parse infile

# index of SYMBOL column
my $symbolCol = -1;

# header
my $header = <STDIN>;
# immediately print back to stdout
print $header;
chomp($header);
my @headers = split(/\t/,$header);
foreach my $i (0..$#headers) {
    ($headers[$i] eq "SYMBOL") || next;
    $symbolCol = $i;
    last;
}
($symbolCol != -1) || die "E $0: could not find SYMBOL in headers\n";

# max number of fields we want when splitting each data line
my $maxSplit = $symbolCol + 2;

# data lines
while (my $line = <STDIN>) {
    # immediately print
    print $line;
    chomp($line);
    my @fields = split(/\t/, $line, $maxSplit) ;
    ($#fields >= $symbolCol) ||
	die "E $0: data line doesn't have enough fields:\n$line\n";
    my $gene = $fields[$symbolCol];
    (defined $knownCandidateGenes{$gene}) && (delete $knownCandidateGenes{$gene});
}

#############################################
# any remaining key of %knownCandidateGenes was never seen
foreach my $gene (sort(keys(%knownCandidateGenes))) {
    warn "W $0: \"known candidate/causal gene\" $gene for ".$knownCandidateGenes{$gene}.
	" was never seen! typo in samples or candidateGenes xlsx files?\n";
}

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE\n";
