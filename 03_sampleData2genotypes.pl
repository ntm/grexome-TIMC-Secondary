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


# 24/03/2018
# NTM

# Read on stdin a (G)VCF file with one data column per sample.
# FORMAT MUST start with either GT:AF and have DP somewhere (for SNVs 
# and short indels), or it must start with GT:RR and have BF somewhere
# (for CNVs).
# Replace all the sample data columns by genotype columns 
# HV, HET, OTHER, HR.
# Each column has eg: 
# genotype1~sample1,sample2|genotype2~sample3|genotype3~sample4,sample5
# - HR has only one genotype: 0/0.
# - all genotypes in HET column will be REF/VAR, eg 0/1, 0/2...
# - all genotypes in HV column will by homovar, eg 1/1, 2/2 etc...
# - all other genotypes (ie VAR1/VAR2) will appear in OTHER column.
# For HV and HET genos at SNVs/short-indels, "sample" actually 
# becomes: $sample[$dp:$af].
# Similarly for CNV calls "sample" becomes  $sample[$BF:$RR].
#
# <NON_REF> is removed from ALTs if it was present, apart from that
# we don't touch the first 8 columns: we assume all ALTs are called
# in at least one genotype, and the REF+ALTs are already normalized.
# This is the case if the input GVCF went through filterBadCalls.pl .
#
# NOTE that FORMAT becomes "GENOS", so produced file is no longer
# a VCF stricto sensu (spec requires first FORMAT key is GT).


use strict;
use warnings;
use File::Basename qw(basename);
use POSIX qw(strftime);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


# parse header, just copy it except for FORMAT descriptions,
# which we replace with our FORMAT descriptions, and 
# the CHROM line, where we replace the sample data by genotype 
# columns, but store the list of sample names in @samples
my @samples;
while(my $line = <STDIN>) {
    if ($line =~ /^##/) {
	# ignore previous FORMAT descriptions, print other header lines
	($line =~ /^##FORMAT=/) || print $line;
    }
    elsif ($line =~ /^#CHROM/) {
	chomp($line);
	@samples = split(/\t/, $line);
	(@samples >= 10) || die "E $0: no sample ids in CHROM line?\n$line\n";
	# first 9 fields are just copied
	my $lineToPrint = shift(@samples);
	foreach my $i (2..9) {
	    $lineToPrint .= "\t".shift(@samples);
	}
	# now @samples has the sample ids, as we wanted.
	# add the new column names
	$lineToPrint .= "\tHV\tHET\tOTHER\tHR\n";

	# Now print the command line, the new FORMAT descriptions, and 
	# then the new CHROM line
	# add info with full command line run
	my $com = "$0 ".join(" ",@ARGV);
	chomp($com);
	$com .= " < ".`readlink -f /proc/$$/fd/0` ;
	chomp($com);
	$com .= " > ".`readlink -f /proc/$$/fd/1` ;
	chomp($com);
	$com .= " 2> ".`readlink -f /proc/$$/fd/2` ;
	chomp($com);
	print "##sampleData2genotypes=<commandLine=\"$com\">\n";
	print "##FORMAT=<ID=GENOS,Number=.,Type=String,Description=\"pipe-separated list of genotypes called at this position, with all corresponding sample identifiers for each genotype (with [DP:AF] for HET and HV SNVs/short-indels, and [BF:RR] for CNVs)\">\n";
	print $lineToPrint;
	last;
    }
    else {
	die "E $0: parsing header, found bad line:\n$line";
    }
}

# now parse data lines
while(my $line = <STDIN>) {
    chomp($line);

    my @data = split(/\t/, $line);
    (@data >= 10) || die "E $0: no sample data in line?\n$line\n";
    # first 8 fields are copied, except we remove <NON_REF> from ALTs if it was present
    $data[4] =~ s/,<NON_REF>//;
    my $lineToPrint = shift(@data);
    foreach my $i (1..7) {
	$lineToPrint .= "\t".shift(@data);
    }

    # from FORMAT we need GT:AF+DP (assumption: filterBadCalls.pl already ran
    # and created/updated AF and DP with correct values), or GT:RR+BF for CNVs.
    my $format = shift(@data);
    # GT and AF/RR should always be first, check it
    # Also disguise BF as DP since it's syntactically identical
    if ($format =~ /^GT:RR:/) {
	($format =~ s/:BF/:DP/) ||
	    die "E $0: FORMAT starts with GT:RR but doesn't contain BF in:\n$line\n";
    }
    elsif ($format !~ /^GT:AF:/) {
	die "E $0: FORMAT doesn't begin with 'GT:AF:' or 'GT:RR:' in:\n$line\n";
    }
    # find DP index
    my $dpCol = 0;
    my @format = split(/:/, $format);
    foreach my $i (2..$#format) {
	($format[$i] eq "DP") && ($dpCol = $i) && last;
    }
    ($dpCol==0) && die "E $0: DP/BF not found in format:\n$line\n";
    # print new FORMAT
    $lineToPrint .= "\tGENOS";

    # sanity check
    (@data == @samples) || 
	die "E $0: we don't have same number of samples and data columns in:\n$line\n";
    # construct comma-separated list of samples for each called genotype
    my %geno2samples = ();

    foreach my $i (0..$#data) {
	# add trailing ':' so we know GT is followed by : (eg for NOCALLS)
	$data[$i] .= ':';
	# ignore NOCALLs
	($data[$i] eq './.:') && next;
	my ($geno,$af);
	if ($data[$i] eq '0/0:') {
	    # HR at a CNV, disguise as normal HR ie set $af='.'
	    ($geno,$af) = ('0/0','.');
	}
	else {
	    ($data[$i] =~ m~^(\d+/\d+):([^:]+):~) || 
		die "E $0: cannot grab genotype and AF/RR for sample $i in $data[$i] in line:\n$line\n";
	    ($geno,$af) = ($1,$2);
	}

	if (defined $geno2samples{$geno}) { $geno2samples{$geno} .= ","; }
	else { $geno2samples{$geno} = ""; }
	$geno2samples{$geno} .= $samples[$i];
	if ($af ne '.') {
	    # if we have an AF this is a HET or HV call (possibly a CNV), find DP/BF
	    my @thisData = split(/:/, $data[$i]);
	    my $dp = 0;
	    ($dpCol) && ($thisData[$dpCol]) && ($thisData[$dpCol] ne '.') && ($dp = $thisData[$dpCol]);
	    ($dp) || die "E $0: AF/RR is $af but couldn't find DP/BF in $data[$i]\n$line\n";
	    # add [DP:AF] / [BF:RR] after sample ID
	    $geno2samples{$geno} .= "[$dp:$af]";
	}
    }

    # now print data for each genotype, removing data from %geno2samples as we go
    # format is: genotype1~sample1,sample2|genotype2~sample3|genotype3~sample4,sample5
    # HV
    $lineToPrint .= "\t";
    foreach my $geno (sort GENOSORT keys %geno2samples) {
	($geno eq '0/0') && next; # skip HR
	($geno =~ m~^(\d+)/\1$~) || next;
	$lineToPrint .= "$geno~".$geno2samples{$geno}."|";
	delete($geno2samples{$geno});
    }
    # remove last | if needed (not needed if there are no HVs)
    $lineToPrint =~ s/\|$// ;

    # HET == 0/*
    $lineToPrint .= "\t";
    foreach my $geno (sort GENOSORT keys %geno2samples) {
	($geno eq '0/0') && next; # skip HR
	($geno =~ m~^0/\d+$~) || next;
	$lineToPrint .= "$geno~".$geno2samples{$geno}."|";
	delete($geno2samples{$geno});
    }
    # remove last | if there was at least one HET geno
    $lineToPrint =~ s/\|$// ;

    # OTHER
    $lineToPrint .= "\t";
    foreach my $geno (sort GENOSORT keys %geno2samples) {
	($geno eq '0/0') && next; # skip HR
	($geno =~ m~^\d+/\d+$~) || die "E $0: cannot parse geno $geno in line:\n$line\n";
	$lineToPrint .= "$geno~".$geno2samples{$geno}."|";
	delete($geno2samples{$geno});
    }
    # remove last | if needed
    $lineToPrint =~ s/\|$// ;

    # HR: if no sample is HR, column is empty
    $lineToPrint .= "\t";
    if (defined $geno2samples{'0/0'}) {
	$lineToPrint .= "0/0~".$geno2samples{'0/0'};
	delete($geno2samples{'0/0'});
    }

    # sanity:
    (%geno2samples) &&
	die "E $0: after eating every geno, geno2samples not empty: %geno2samples, line:\n$line\n";

    print $lineToPrint."\n";
}

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";


#######################################
# sort function for genotypes x/y (x and y are ints)
sub GENOSORT {
    ($a =~ m~^(\d+)/(\d+)$~) || die "E $0 in GENOSORT: cannot decompose a $a\n";
    my ($a1,$a2) = ($1,$2);
    ($b =~ m~^(\d+)/(\d+)$~) || die "E $0 in GENOSORT: cannot decompose b $b\n";
    my ($b1,$b2) = ($1,$2);
    return(($a1 <=> $b1) || ($a2 <=> $b2));
}
