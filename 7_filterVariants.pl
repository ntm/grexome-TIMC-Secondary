#!/usr/bin/perl

# 26/03/2018
# NTM

# Parses on stdin a TSV file, produced by extractCohorts.pl
# or extractSamples.pl for example.
# Applies a bunch of filters (see args), and prints to stdout
# a similar file but where some lines have been filtered out.


use strict;
use warnings;
use Getopt::Long;

# default values for args
my $max_ctrl_hv = 3; # COUNT_NEGCTRL_HV <= $x
my $max_ctrl_het = 10; # COUNT_NEGCTRL_HET <= $x
my $min_cohort_hv = 0; # COUNT_$cohort_HV >= $x
my $min_hr = 100; # COUNT_$cohort_HR + COUNT_NEGCTRL_HR >= $x

my $no_mod = ''; # if enabled, filter out MODIFIER impacts
my $no_low = ''; # if enabled, filter out LOW impacts

my $pick = ''; # if enabled, only keep lines with PICK
# could add a filter on BIOTYPE values (eg protein_coding,
# processed_transcript, retained_intron, nonsense_mediated_decay)
# Not implementing now.

my $max_af_gnomad = 0.01; # gnomAD_AF <= $x
my $max_af_1kg = 0.03; # AF <= $x, this is 1KG phase 3
my $max_af_esp = 0.05; # AA_AF and EA_AF <= $x, this is ESP
GetOptions ("max_ctrl_hv=i" => \$max_ctrl_hv,
	    "max_ctrl_het=i" => \$max_ctrl_het,
	    "min_cohort_hv=i" => \$min_cohort_hv,
	    "min_hr=i" => \$min_hr,
	    "no_mod" => \$no_mod,
	    "no_low" => \$no_low,
	    "pick" => \$pick,
	    "max_af_gnomad=f" => \$max_af_gnomad,
	    "max_af_1kg=f" => \$max_af_1kg,
	    "max_af_esp=f" => \$max_af_esp)
    or die("Error in command line arguments\n");

# build string of all filter values, for logging
my $filterString = "max_ctrl_hv=$max_ctrl_hv max_ctrl_het=$max_ctrl_het";
($min_cohort_hv) && ($filterString .= " min_cohort_hv=$min_cohort_hv");
($min_hr) && ($filterString .= " min_hr=$min_hr");
($no_mod) && ($filterString .= " no_mod");
($no_low) && ($filterString .= " no_low");
($pick) && ($filterString .= " pick");
$filterString .= " max_af_gnomad=$max_af_gnomad max_af_1kg=$max_af_1kg max_af_esp=$max_af_esp";


# copy header, adding a final column with all filter values
my $header = <STDIN>;
chomp($header);
print "$header\t$filterString\n" ;
# build hash of header titles, value is the column number for that header
my %title2index;
my @titles = split(/\t/, $header);
foreach my $i (0..$#titles) {
    my $title = $titles[$i];
    if ($title =~ /^COUNT_(\w+)_HV/) {
	# replace cohort name with COHORT in COUNT_*_HV (only for the hash key, not in the outFile)
	($1 ne "NEGCTRL") && ($title = "COUNT_COHORT_HV");
    }
    elsif ($title =~ /^COUNT_(\w+)_HR/) {
	# same for HR, replace cohort name with COHORT in COUNT_*_HR hash key
	($1 ne "NEGCTRL") && ($title = "COUNT_COHORT_HR");
    }
    $title2index{$title} = $i;
}


# parse data
while(my $line = <STDIN>) {
    chomp($line);
    my @fields = split(/\t/, $line, -1);
    # apply all filters
    if (($pick) && (! $fields[$title2index{"PICK"}])) {
	next;
    }
    if ($fields[$title2index{"COUNT_NEGCTRL_HV"}] > $max_ctrl_hv) {
	next;
    }
    if ($fields[$title2index{"COUNT_COHORT_HV"}] < $min_cohort_hv) {
	next;
    }
    if ($fields[$title2index{"COUNT_COHORT_HR"}] + $fields[$title2index{"COUNT_NEGCTRL_HR"}]  < $min_hr) {
	next;
    }
    if ((defined $max_ctrl_het) && ($fields[$title2index{"COUNT_NEGCTRL_HET"}] > $max_ctrl_het)) {
	next;
    }
    if (($no_mod) && ($fields[$title2index{"IMPACT"}] eq "MODIFIER")) {
	next;
    }
    if (($no_low) && ($fields[$title2index{"IMPACT"}] eq "LOW")) {
	next;
    }
    if ($fields[$title2index{"gnomAD_AF"}]) {
	# sometimes we have several &-separated values, in this case
	# only filter if all values are high
	my $keep = 0;
	foreach my $gnomad (split(/&/, $fields[$title2index{"gnomAD_AF"}])) {
	    ($gnomad <= $max_af_gnomad) && ($keep = 1);
	}
	($keep) || next;
    }
    if ($fields[$title2index{"AF"}]) {
	#again several &-separated values
	my $keep = 0;
	foreach my $af (split(/&/, $fields[$title2index{"AF"}])) {
	    ($af <= $max_af_1kg) && ($keep = 1);
	}
	($keep) || next;
    }
    if ($fields[$title2index{"AA_AF"}]) {
	#again several &-separated values
	my $keep = 0;
	foreach my $esp (split(/&/, $fields[$title2index{"AA_AF"}])) {
	    ($esp <= $max_af_esp) && ($keep = 1);
	}
	($keep) || next;
    }
    if ($fields[$title2index{"EA_AF"}]) {
	#again several &-separated values
	my $keep = 0;
	foreach my $esp (split(/&/, $fields[$title2index{"EA_AF"}])) {
	    ($esp <= $max_af_esp) && ($keep = 1);
	}
	($keep) || next;
    }

    # passed all filters
    print "$line\n";
}
