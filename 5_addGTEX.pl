#!/usr/bin/perl

# 25/03/2018
# NTM

# Take as args a list of "favorite tissues", these should be tissue names
# as they appear in the header of &gtexFile() (from grexome_se_config.pm)
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

#use lib '/home/nthierry/PierreRay/Grexome/SecondaryAnalyses/';
# export PERL5LIB=... before calling script, it's more flexible
use grexome_sec_config qw(gtexFile);

my @favoriteTissues = @ARGV;
my $gtexFile = &gtexFile();

# the GTEX columns will be placed just before the $insertBefore column,
# which must exist (this is checked)
my $insertBefore = "HV";

# GTEX data will be stored in hash:
# key is ENSG, value is a ref to an array holding the TPM values
my %gtex;

(-f $gtexFile) || die "E: gtex file $gtexFile doesn't exist\n";
open(GTEX, $gtexFile) || die "E: cannot open gtex file $gtexFile for reading\n";
# skip header
foreach my $i (1..4) {
    my $line = <GTEX>;
}

# next line should be column names, check and grab tissue names
# and index of $favoriteTissue in @tissues
my $line = <GTEX>;
chomp $line;
($line =~ s/^Gene ID\tGene Name\t//) || 
    die "E: line should be GTEX header but can't parse it:\n$line\n" ;
my @tissues = split(/\t/, $line);
my @favTissIndex = (-1) x @favoriteTissues;
# improve tissue name strings, for our header line
foreach my $i (0..$#tissues) {
    foreach my $fti (0..$#favoriteTissues) {
	if ($tissues[$i] eq $favoriteTissues[$fti]) {
	    ($favTissIndex[$fti] != -1) && die "found favorite tissue $tissues[$i] twice, WTF!\n";
	    $favTissIndex[$fti] = $i;
	    last;
	}
    }
    $tissues[$i] =~ s/ /_/g ;
    $tissues[$i] = "GTEX_".$tissues[$i];
}

foreach my $fti (0..$#favoriteTissues) {
    ($favTissIndex[$fti] == -1) && die "could not find favorite tissue $favoriteTissues[$fti] column in header.\n";
}

# now parse GTEX datalines
while ($line=<GTEX>) {
    chomp $line;
    # set LIMIT=-1 so we still produce (empty) array elements if we have trailing empty fields
    my @data = split(/\t/, $line, -1);
    (@data == @tissues+2) || die "E: wrong number of fields in:\n$line\n";
    # grab ENSG, ignore gene name
    my $ensg = shift(@data);
    shift(@data);
    $gtex{$ensg} = \@data ;
}

# OK, now parse infile
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
    die "E: could not find insertBefore==$insertBefore in column headers of infile:\n$header\n";
($geneIndex) || die "E: could not find Gene field in column headers:\n$header\n";

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
    ($gene =~ /,/) && die "line in inFile has several Genes, shouldn't happen:\n$line\n";

    # @thisGtex: array of strings holding expression values, one per tissue
    my @thisGtex =("") x @tissues ;
    # calculated gtex ratio for each favorite tissue
    my @favTissRatios = ("") x @favoriteTissues;

    if (defined $gtex{$gene}) {
	# for calculating GTEX_*_RATIO:
	# sum of all gtex values
	my $sumOfGtex = 0;
	# number of defined gtex values
	my $numberOfGtex = 0;
	my @gtexThisGene = @{$gtex{$gene}};
	foreach my $i (0..$#gtexThisGene) {
	    ($gtexThisGene[$i]) || next;
	    $thisGtex[$i] = $gtexThisGene[$i];
	    $sumOfGtex += $gtexThisGene[$i] ;
	    $numberOfGtex++;
	}

	# favExp / averageExp == favExp / (sumExp / numberGtex) == favExp * numberGtex / sumExp
	# so make sure we can divide by $sumOfGtex
	($sumOfGtex) || die "Sum of GTEX values is zero for gene $gene, impossible?\n$line\n";
	foreach my $ii (0..$#favTissIndex) {
	    ($thisGtex[$favTissIndex[$ii]]) && 
		($favTissRatios[$ii] = $thisGtex[$favTissIndex[$ii]] * $numberOfGtex / $sumOfGtex) ;
	}
    }
    # else: no expression data for $gene, use the default empty strings

    # OK, build line with expression values inserted where they should, 
    # and with GTEX_RATIOs first, then favorites, then others
    my @toPrint = @fields[0..$insertBeforeIndex-1];
    push(@toPrint, @favTissRatios);
    push(@toPrint, @thisGtex[@favTissIndex]);
    foreach my $i (0..$#tissues) {
	(grep(/^$i$/, @favTissIndex) == 0) && push(@toPrint, $thisGtex[$i]);
    }
    push(@toPrint, @fields[$insertBeforeIndex..$#fields]);

    print join("\t",@toPrint)."\n";
}