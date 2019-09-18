#!/usr/bin/perl

# 18/09/2019
# NTM

# Parse on stdin a TSV file produced by extractCohorts.pl,
# preferably filtered by 7_filterVariants.pl.
# Print to stdout a similar file but where the order of
# columns has been changed:
# - columns from @newOrder are printed first, in that order;
# - all other columns are printed in the order they are seen.


use strict;
use warnings;



# columns to be printed first, in that order (replace COHORT
# with the current cohort)
my @newOrder = qw(POSITION REF ALT SYMBOL KNOWN_CANDIDATE_GENE COUNT_COHORT_HV COUNT_COHORT_HET COUNT_COHORT_HR COUNT_NEGCTRL_HV COUNT_NEGCTRL_HET COUNT_NEGCTRL_HR HV HET IMPACT Consequence HGVSc HGVSp Protein_position SIFT PolyPhen CADD_raw_rankscore MutationTaster_pred gnomAD_AF GTEX_testis_RATIO GTEX_ovary_RATIO GTEX_testis GTEX_ovary GTEX_blood GTEX_cerebellar_hemisphere GTEX_liver GTEX_lung);


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

my $header = <STDIN>;
chomp($header);
my @titles = split(/\t/, $header);
foreach my $i (0..$#titles) {
    my $title = $titles[$i];
    if ($title =~ /^COUNT_(\w+)_H[VRE]/) {
	# replace cohort name with COHORT in COUNT_*_HV HR HET 
	# (only in $title not in @titles)
	($1 ne "NEGCTRL") && ($title =~ s/$1/COHORT/);
    }

    if (defined $title2index{$title}) {
	$old2new[$i] = $title2index{$title};
	# delete so we can make sure every @newOrder title was seen
	delete($title2index{$title});
    }
    else {
	$old2new[$i] = $nextNotNewOrder;
	$nextNotNewOrder++;
    }

}

# sanity: make sure every @newOrder column exists
if (my @missingTitles = keys(%title2index)) {
    die "E: some newOrder titles were not found: ".join(" ", keys(%title2index))."\n";
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

