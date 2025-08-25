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

# Parses on stdin a TSV file produced by extractCohorts.pl or extractSamples.pl,
# or even straight out of 05_vcf2tsv.pl (but then don't filter on COUNT*).
# Applies a bunch of filters (see args), and prints to stdout
# a similar file but where some lines have been filtered out.


use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename qw(basename);
use Getopt::Long;

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


# arguments for filtering: no default values, all filters disabled by default
my $max_ctrl_hv; # COUNT_NEGCTRL_HV <= $x
my $max_ctrl_het; # COUNT_NEGCTRL_HET <= $x
my $min_cohort_hv; # COUNT_$cohort_HV >= $x
my $min_hr; # COUNT_HR >= $x

my $no_mod = ''; # if enabled, filter out MODIFIER impacts
my $no_low = ''; # if enabled, filter out LOW impacts
my $no_pseudo = ''; # if enabled, filter out all *pseudogene BIOTYPEs
my $no_nmd = ''; # if enabled, filter out nonsense_mediated_decay BIOTYPE
my $canon = ''; # if enabled, only keep lines with CANONICAL==YES

# max global alt allele freqs in gnomADe, gnomADg, RegeneronME, AllofUs, ALFA (if available)
my $max_af_global;
# max per-population alt allele freqs in AllofUs and dbNSFP (if available)
my $max_af_perPop;

# if true, print timestamped start-done log messages to stderr
my $logTime = '';

GetOptions ("max_ctrl_hv=i" => \$max_ctrl_hv,
	    "max_ctrl_het=i" => \$max_ctrl_het,
	    "min_cohort_hv=i" => \$min_cohort_hv,
	    "min_hr=i" => \$min_hr,
	    "no_mod" => \$no_mod,
	    "no_low" => \$no_low,
	    "no_pseudo" => \$no_pseudo,
	    "no_nmd" => \$no_nmd,
	    "canonical" => \$canon,
	    "max_af_global=f" => \$max_af_global,
	    "max_af_perPop=f" => \$max_af_perPop,
	    "logtime" => \$logTime)
    or die("E $0: Error in command line arguments\n");

if ($logTime) {
    my $now = strftime("%F %T", localtime);
    warn "I $now: $0 - starting to run\n";
}

# build string of all filter values, for logging
my $filterString = "";
($max_ctrl_hv) && ($filterString .= "max_ctrl_hv=$max_ctrl_hv ");
($max_ctrl_het) && ($filterString .= "max_ctrl_het=$max_ctrl_het ");
($min_cohort_hv) && ($filterString .= "min_cohort_hv=$min_cohort_hv ");
($min_hr) && ($filterString .= "min_hr=$min_hr ");
($no_mod) && ($filterString .= "no_mod ");
($no_low) && ($filterString .= "no_low ");
($no_pseudo) && ($filterString .= "no_pseudo ");
($no_nmd) && ($filterString .= "no_nmd ");
($canon) && ($filterString .= "canonical ");
($max_af_global) && ($filterString .= "max_af_global=$max_af_global ");
($max_af_perPop) && ($filterString .= "max_af_perPop=$max_af_perPop ");
# remove trailing space and add leading tab if any filters are applied
if ($filterString) {
    $filterString = "\t$filterString";
    chop($filterString);
}

# copy header, adding a column with all filter values
my $header = <STDIN>;
chomp($header);
print "$header$filterString\n" ;
# build hash of header titles, value is the column number for that header
my %title2index;
my @titles = split(/\t/, $header);
foreach my $i (0..$#titles) {
    my $title = $titles[$i];
    if (($max_ctrl_hv) || ($max_ctrl_het) || ($min_cohort_hv) || ($min_hr)) {
	# if filtering on COUNTs we will need COUNT_$cohort_HV , we want a uniform hash
	# key COUNT_COHORT_HV in %title2index but we must ignore the other COUNT_*_HV columns
	if (($title !~ /_NEGCTRL_/) && ($title !~ /_COMPAT_/) && ($title !~ /_OTHERCAUSE_/) &&
	    ($title =~ /^COUNT_(\w+)_HV/)) {
	    # OK replace cohort name with COHORT as hash key
	    $title = "COUNT_COHORT_HV";
	}
    }
    # sanity
    (defined $title2index{$title}) &&
	die "E $0: title $title defined twice\n";
    $title2index{$title} = $i;
}

# make sure all titles we need are present
if (($max_ctrl_hv) || ($max_ctrl_het) || ($min_cohort_hv) || ($min_hr)) {
    foreach my $t ("COUNT_NEGCTRL_HV","COUNT_NEGCTRL_HET","COUNT_COHORT_HV","COUNT_HR") {
	(defined $title2index{$t}) ||
	    die "E $0: title $t required by script but missing, some VEP columns changed?\n";
    }
}
if (($no_mod) || ($no_low)) {
    (defined $title2index{"IMPACT"}) ||
	die "E $0: IMPACT required by script but missing, some VEP columns changed?\n";
}
if (($no_pseudo) || ($no_nmd)) {
    (defined $title2index{"BIOTYPE"}) ||
	die "E $0: BIOTYPE required by script but missing, some VEP columns changed?\n";
}
if ($canon) {
    (defined $title2index{"CANONICAL"}) ||
	die "E $0: CANONICAL required by script but missing, some VEP columns changed?\n";
}
if ($max_af_global) {
    foreach my $t ("gnomADe_AF","gnomADg_AF","RegeneronME_ALL_AF","AllofUs_ALL_AF","ALFA_Total_AF") {
	(defined $title2index{$t}) ||
	    die "E $0: $t required by script but missing, some VEP columns changed?\n";
    }
}
if ($max_af_perPop) {
    foreach my $t ("AllofUs_POPMAX_AF","dbNSFP_POPMAX_AF") {
	(defined $title2index{$t}) ||
	    die "E $0: $t required by script but missing, some VEP columns changed?\n";
    }
}

# parse data
while(my $line = <STDIN>) {
    chomp($line);
    my @fields = split(/\t/, $line, -1);
    # apply all filters
    if (($canon) && ($fields[$title2index{"CANONICAL"}] ne 'YES')) {
	next;
    }
    if (($no_mod) && ($fields[$title2index{"IMPACT"}] eq "MODIFIER")) {
	next;
    }
    if (($no_low) && ($fields[$title2index{"IMPACT"}] eq "LOW")) {
	next;
    }
    if (($no_pseudo) && ($fields[$title2index{"BIOTYPE"}] =~ /pseudogene$/)) {
	# there are a bunch of pseudogene biotypes but they all end with 'pseudogene'
	next;
    }
    if (($no_nmd) && ($fields[$title2index{"BIOTYPE"}] eq "nonsense_mediated_decay")) {
	next;
    }
    if ((defined $max_ctrl_hv) && ($fields[$title2index{"COUNT_NEGCTRL_HV"}] > $max_ctrl_hv)) {
	next;
    }
    if ((defined $max_ctrl_het) && ($fields[$title2index{"COUNT_NEGCTRL_HET"}] > $max_ctrl_het)) {
	next;
    }
    if ((defined $min_cohort_hv) && ($fields[$title2index{"COUNT_COHORT_HV"}] < $min_cohort_hv)) {
	next;
    }
    if  ((defined $min_hr) && ($fields[$title2index{"COUNT_HR"}]  < $min_hr)) {
	next;
    }
    if (defined $max_af_global) {
	my $keep = 1;
	foreach my $title ("gnomADe_AF","gnomADg_AF","RegeneronME_ALL_AF",
			   "AllofUs_ALL_AF","ALFA_Total_AF") {
	    if ($fields[$title2index{$title}]) {
		# sometimes we have several &-separated values, filter if any is high
		foreach my $af (split(/&/, $fields[$title2index{$title}])) {
		    if ($af > $max_af_global) {
			$keep = 0;
			last;
		    }
		}
		# if AF was too high, don't even look at other sources
		($keep) || last;
	    }
	}
	($keep) || next;
    }
    if (defined $max_af_perPop) {
	my $keep = 1;
	foreach my $title ("AllofUs_POPMAX_AF","dbNSFP_POPMAX_AF") {
	    if (($fields[$title2index{$title}]) && ($fields[$title2index{$title}] > $max_af_perPop)) {
		$keep = 0;
		last;
	    }
	}
	($keep) || next;
    }

    # passed all filters
    print "$line\n";
}

if ($logTime) {
    my $now = strftime("%F %T", localtime);
    warn "I $now: $0 - ALL DONE, completed successfully!\n";
}
