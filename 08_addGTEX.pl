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


# 25/03/2018
# NTM

# Take as args:
# - a GTEX file with full path, eg E-MTAB-5214-query-results.tpms.tsv
#   from the GTEX_Data/ subdir;
# - a list of "favorite tissues", at least one, these should be tissue names
#   as they appear in the header of $gtexFile but with spaces replaced with underscores.
# Read on stdin a TSV file as produced by vcf2tsv.pl, 
# print to stdout a TSV file with added columns holding
# GTEX TPM values taken from $gtexFile.
# The tissue names found in $gtexFile are changed: s/ /_/g , and "GTEX_"
# is prepended.
# We use the ENSG identifiers (they are non-redundant in $gtexFile).
# 06/04/2018: We now expect (and check) that the Gene column in inFile
# contains a single ENSG.
#
# New columns GTEX_$favoriteTissue_RATIO are added before other GTEX
# columns, holding GTEX_$favTiss / (sum of all gtex expression values)
# for each $favTiss in @ARGV.
# The GTEX columns are placed before the $insertBefore column (must exist), with 
# GTEX_$favoriteTissue(s) placed first.

use strict;
use warnings;
use File::Basename qw(basename);
use Getopt::Long;
use POSIX qw(strftime);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## hard-coded stuff that shouldn't change much

# the GTEX columns will be placed just before the $insertBefore column,
# which must exist (this is checked)
my $insertBefore = "HV";


#############################################
## options / params from the command-line

# full path to the GTEX datafile, eg GTEX_Data/E-MTAB-5214-query-results.tpms.tsv
my $gtexFile;

# comma-separated list of favorite tissues with a default (working on infertility here)
my $favoriteTissues = "testis,ovary";

# help: if true just print $USAGE and exit
my $help = '';


my $USAGE = "\nParse on STDIN a TSV file as produced by steps 1-7 of this secondaryAnalysis pipeline; print to STDOUT a similar file with additional columns holding the GTEX expression data.\n
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--gtex string [no default] : TSV file holding the GTEX TPM values, with path
--favoriteTissues string [default=$favoriteTissues] : comma-separated list of tissues of interest
--help : print this USAGE";

GetOptions ("gtex=s" => \$gtexFile,
            "favoriteTissues=s" => \$favoriteTissues,
            "help" => \$help)
    or die("E: $0 - Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($gtexFile) || die "E $0: you must provide a GTEX file. Try $0 --help\n";
(-f $gtexFile) || die "E $0: the supplied GTEX file doesn't exist\n";

my @favoriteTissues = split(/,/,$favoriteTissues);
(@favoriteTissues) || 
    die "E $0: we expect at least one favoriteTissue, just use defaults and ignore them if you don't care\n";

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";

#############################################
## parse GTEX file

# GTEX data will be stored in hash:
# key is ENSG, value is a ref to an array of strings (possibly "")
# containing the new GTEX RATIOs, then the GTEX TPMs with
# @favoriteTissues first.
my %gtex;


open(GTEX, $gtexFile) || die "E $0: cannot open gtex file $gtexFile for reading\n";
# skip header
foreach my $i (1..4) {
    my $line = <GTEX>;
}

# next line should be column names, check and grab tissue names
# and index of $favoriteTissue in @tissues
my $line = <GTEX>;
chomp $line;
($line =~ s/^Gene ID\tGene Name\t//) || 
    die "E $0: line should be GTEX header but can't parse it:\n$line\n" ;
my @tissues = split(/\t/, $line);
my @favTissIndex = (-1) x @favoriteTissues;
# improve tissue name strings, for our header line
foreach my $i (0..$#tissues) {
    $tissues[$i] =~ s/ /_/g ;
    foreach my $fti (0..$#favoriteTissues) {
        if ($tissues[$i] eq $favoriteTissues[$fti]) {
            ($favTissIndex[$fti] != -1) && die "E $0: found favorite tissue $tissues[$i] twice, WTF!\n";
            $favTissIndex[$fti] = $i;
            last;
        }
    }
    $tissues[$i] = "GTEX_".$tissues[$i];
}

foreach my $fti (0..$#favoriteTissues) {
    ($favTissIndex[$fti] == -1) && die "E $0: could not find favorite tissue $favoriteTissues[$fti] column in header.\n";
}

# now parse GTEX datalines
while ($line=<GTEX>) {
    chomp $line;
    # set LIMIT=-1 so we still produce (empty) array elements if we have trailing empty fields
    my @data = split(/\t/, $line, -1);
    (@data == @tissues+2) || die "E $0: wrong number of fields in:\n$line\n";
    # grab ENSG, ignore gene name
    my $ensg = shift(@data);
    shift(@data);
    ($gtex{$ensg}) &&
        die "E $0: ENSG $ensg present twice in GTEX file\n";
    
    # @thisGtex: array of strings holding expression values, one per tissue
    my @thisGtex =("") x @tissues ;
    # calculated gtex ratio for each favorite tissue
    my @favTissRatios = ("") x @favoriteTissues;

    # for calculating GTEX_*_RATIO:
    # sum of all gtex values
    my $sumOfGtex = 0;
    foreach my $i (0..$#data) {
        ($data[$i]) || next;
        $thisGtex[$i] = $data[$i];
        $sumOfGtex += $data[$i] ;
    }

    # favExp / averageExp == favExp / (sumExp / nbTissues) == favExp * nbTissues / sumExp
    # so make sure we can divide by $sumOfGtex
    ($sumOfGtex) || die "E $0: Sum of GTEX values is zero for gene $ensg, impossible?\n$line\n";
    foreach my $ii (0..$#favTissIndex) {
        ($thisGtex[$favTissIndex[$ii]]) && 
            ($favTissRatios[$ii] = $thisGtex[$favTissIndex[$ii]] * @tissues / $sumOfGtex) ;
    }
    
    # OK build array of strings with GTEX_RATIOs first, then favorites, then others
    my @toPrint = ();
    foreach my $favTissRatio (@favTissRatios) {
        # print max 2 digits after decimal
        my $favTR_toPrint = $favTissRatio;
        ($favTissRatio) && ($favTR_toPrint = sprintf("%.2f",$favTissRatio));
        push(@toPrint, $favTR_toPrint);
    }
    push(@toPrint, @thisGtex[@favTissIndex]);
    foreach my $i (0..$#tissues) {
        (grep(/^$i$/, @favTissIndex) == 0) && push(@toPrint, $thisGtex[$i]);
    }
    $gtex{$ensg} = \@toPrint;
}

close(GTEX);


#############################################
# parse infile

# header: add tissue names just before $insertBefore
my $header = <STDIN>;
chomp($header);
my @headers = split(/\t/,$header);
# find the index in @headers of $insertBefore and of "Gene" columns
my $insertBeforeIndex;
my $geneIndex;
foreach my $i (0..$#headers) {
    ($headers[$i] eq "Gene") && ($geneIndex = $i);
    ($headers[$i] eq $insertBefore) && ($insertBeforeIndex = $i);
}
# make sure we found them
($insertBeforeIndex) ||
    die "E $0: could not find insertBefore==$insertBefore in column headers of infile:\n$header\n";
($geneIndex) || die "E $0: could not find Gene field in column headers:\n$header\n";

# now make new header
my @newHeaders = @headers[0..$insertBeforeIndex-1];
# add new GTEX RATIO headers
foreach my $favTis (@favoriteTissues) {
    push(@newHeaders, "GTEX_$favTis"."_RATIO") ;
}
# add GTEX headers, @favoriteTissues first
foreach my $fti (@favTissIndex) {
    push(@newHeaders, $tissues[$fti]);
}
foreach my $i (0..$#tissues) {
    (grep(/^$i$/, @favTissIndex) == 0) && push(@newHeaders, $tissues[$i]);
}
push(@newHeaders, @headers[$insertBeforeIndex..$#headers]);

print join("\t",@newHeaders)."\n";


# data lines
while (my $line = <STDIN>) {
    chomp($line);
    my @fields = split(/\t/, $line, -1) ;
    my $gene = $fields[$geneIndex];
    ($gene =~ /,/) && die "E $0: line in inFile has several Genes, shouldn't happen:\n$line\n";

    # build line with expression values inserted where they should
    my @toPrint = @fields[0..$insertBeforeIndex-1];
    if (defined $gtex{$gene}) {
        push(@toPrint, @{$gtex{$gene}});
    }
    else {
        # no expression data for $gene, use empty strings
        push(@toPrint, ("") x (@favoriteTissues + @tissues)) ;
    }
    push(@toPrint, @fields[$insertBeforeIndex..$#fields]);
    print join("\t",@toPrint)."\n";
}

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";
