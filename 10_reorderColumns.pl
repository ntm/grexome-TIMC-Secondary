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


# 18/09/2019
# NTM

# Takes two arguments: $verbose (0 or 1), and a comma-separated list of favorite
# GTEX tissues.
# Parse on stdin a TSV file produced by extractCohorts.pl,
# preferably filtered by filterVariants.pl.
# Print to stdout a similar file but where the order of
# columns has been changed:
# - columns from @newOrder (if present) are printed first, in that order;
# - all other columns are printed in the order they are seen.


use strict;
use warnings;
use File::Basename qw(basename);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

(@ARGV == 2) || 
    die "E $0: needs two arguments: a boolean (verbose), and a comma-separated list of GTEX tissues\n";
my $verbose = $ARGV[0];
my $favoriteTissues = $ARGV[1];
# no sanity checking, these favorite tissues should have been passed to 08_addGTEX.pl 
# earlier and therefore already checked

# columns to be printed first, in that order (replace COHORT
# with the current cohort)
my @newOrder = qw(POSITION REF ALT SYMBOL KNOWN_CANDIDATE_GENE MANE_SELECT Gene Feature CANONICAL BIOTYPE COUNT_HR COUNT_COHORT_HV COUNT_COHORT_HET COUNT_COHORT_OTHERCAUSE_HV COUNT_COHORT_OTHERCAUSE_HET COUNT_COMPAT_HV COUNT_COMPAT_HET COUNT_NEGCTRL_HV COUNT_NEGCTRL_HET COUNT_OTHERGENO IMPACT Consequence CLIN_SIG HGVSc HGVSp cDNA_position CDS_position Protein_position gnomADe_AF gnomADg_AF ALFA_Total_AF RegeneronME_ALL_AF AllofUs_ALL_AF AllofUs_POPMAX_AF dbNSFP_POPMAX_AF dbNSFP_POPMAX_POP COHORT_HV COHORT_HET COHORT_OTHERCAUSE_HV COHORT_OTHERCAUSE_HET COMPAT_HV COMPAT_HET);
foreach my $favTissue (split(/,/, $favoriteTissues)) {
    push(@newOrder, "GTEX_$favTissue"."_RATIO");
}
foreach my $favTissue (split(/,/, $favoriteTissues)) {
    push(@newOrder, "GTEX_$favTissue");
}
# some additional tissues that we like to see early, hard-coded
push(@newOrder, qw(GTEX_blood GTEX_cerebellar_hemisphere GTEX_liver GTEX_lung));

# grab STDIN header
my $header = <STDIN>;
chomp($header);
my @titles = split(/\t/, $header);

# purge from @newOrder any fields that are missing (ignoring *COHORT* fields,
# these will be transformed and checked later)
{
    my %titles;
    foreach my $t (@titles) {
        $titles{$t} = 1;
    }
    my @newOrderMissing = ();
    my @newOrderFound = ();
    foreach my $n (@newOrder) {
        if (($titles{$n}) || ($n =~ /COHORT/)) {
            push(@newOrderFound, $n);
        }
        else {
            push(@newOrderMissing, $n);
        }
    }
    if ($verbose && (@newOrderMissing)) {
        warn "W $0: the following newOrder fields are missing: " . join(' ', @newOrderMissing) . "\n";
        warn "W $0: with non-human data missing fields are expected (eg fields from human-only plugins),\n";
        warn "W $0: but with human data this suggests that the code needs to be updated (open a github issue)\n";
    }
    @newOrder = @newOrderFound;
}

# build hash of @newOrder headers, value is the new column index
# for that header
my %title2index;
foreach my $i (0..$#newOrder) {
    $title2index{$newOrder[$i]} = $i;
}

# mapping from old to new column indexes:
# column that used to be at index $i will be printed at index $old2new[$i]
my @old2new;

# index where we want the next not-newOrder column to go
my $nextNotNewOrder = scalar(@newOrder);

foreach my $i (0..$#titles) {
    my $title = $titles[$i];
    if (($title =~ /^COUNT_(\w+)_OTHERCAUSE_/) || ($title =~ /^(\w+)_OTHERCAUSE_/)) {
        # replace $cohort with COHORT in COUNT_$cohort_OTHERCAUSE_* and
        # in $cohort_OTHERCAUSE_* (only in $title not in @titles)
        $title =~ s/$1/COHORT/;
    }
    elsif (($title !~ /NEGCTRL_/) && ($title !~ /COMPAT_/) &&
           (($title =~ /^COUNT_(\w+)_/) || ($title =~ /^(\w+)_HV$/) || ($title =~ /^(\w+)_HET$/))) {
        # also replace $cohort with COHORT in COUNT_$cohort_* and in $cohort_* ,
        # careful not to touch NEGCTRL or COMPAT or other random titles...
        $title =~ s/$1/COHORT/;
    }

    if (defined $title2index{$title}) {
        ($title2index{$title} == -1) &&
            die "E $0: title $titles[$i] was converted to $title but this has been seen already, fix the regexps!\n";
        $old2new[$i] = $title2index{$title};
        # set to -1 so we can make sure every @newOrder title was seen exactly once
        $title2index{$title} = -1;
    }
    else {
        $old2new[$i] = $nextNotNewOrder;
        $nextNotNewOrder++;
    }
}

# sanity: make sure every @newOrder column was seen exactly once
foreach my $t (keys(%title2index)) {
    ($title2index{$t} == -1) && (delete $title2index{$t});
}
if (my @missingTitles = keys(%title2index)) {
    die "E $0: some newOrder titles were not found: ".join(" ", keys(%title2index))."\n";
}

# new header
my @newHeader;
foreach my $i (0..$#old2new) {
    $newHeader[$old2new[$i]] = $titles[$i];
}
print join("\t", @newHeader)."\n";


# parse data
while(my $line = <STDIN>) {
    chomp($line);
    my @newLine;
    my @fields = split(/\t/, $line, -1);
    foreach my $i (0..$#fields) {
        $newLine[$old2new[$i]] = $fields[$i];
    }
    print join("\t", @newLine)."\n";
}
