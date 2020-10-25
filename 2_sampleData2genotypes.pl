#!/usr/bin/perl

# 24/03/2018
# NTM

# Read on stdin a VCF file with one data column per sample.
# Replace all the sample data columns by genotype columns 
# HV, HET, OTHER, HR.
# Each column has eg: 
# genotype1~sample1,sample2|genotype2~sample3|genotype3~sample4,sample5
# - HR has only one genotype: 0/0.
# - all genotypes in HET column will be REF/VAR, eg 0/1, 0/2...
# - all genotypes in HV column will by homovar, eg 1/1, 2/2 etc...
# - all other genotypes (ie VAR1/VAR2) will appear in OTHER column.
#
# NOTE that FORMAT becomes "GENOS", so produced file is no longer
# a VCF stricto sensu (spec requires first FORMAT key is GT).
#
# UPDATE 09/08/2019: remove ALTs that don't appear in a genotype
# (and renumber remaining ALTs in the called genotypes). 
# UPDATE 14/08/2019: also re-normalize the variants. I don't check for
# collisions (ie if POS got increased by removing leading common bases, it
# might end up at same coord as a subsequent line), but a collision can
# only occur if two REFs overlapped in the input VCF... so the collision
# was already there (even if not explicit with identical POS).
# NOTE: If input variants were normalized, re-normalization should only occur
# when we remove unused ALTs.
# This is an initial cleanup of multiallelic sites: for example with 
# grexomes_0050-0520, 48.4% of sites with 2 ALTs will become monoallelic.
#
# UPDATE 20/08/2019: for HV and HET genos, "sample" is now: $sample[$dp:$af]

use strict;
use warnings;
use File::Basename qw(basename);
use POSIX qw(strftime);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

my $now = strftime("%F %T", localtime);
warn "I $0: $now - starting to run\n";


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
	my $com = qx/ps -o args= $$/;
	chomp($com);
	$com .= " < ".`readlink -f /proc/$$/fd/0` ;
	chomp($com);
	$com .= " > ".`readlink -f /proc/$$/fd/1` ;
	chomp($com);
	$com .= " 2> ".`readlink -f /proc/$$/fd/2` ;
	chomp($com);
	print "##sampleData2genotypes=<commandLine=\"$com\">\n";
	print "##FORMAT=<ID=GENOS,Number=.,Type=String,Description=\"pipe-separated list of genotypes called at this position, with all corresponding sample identifiers for each genotype (with [DP:AF] for HET and HV)\">\n";
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
    # first 5 fields are saved: we may need to modify POS, REF 
    # and ALT when removing unused ALTs and re-normalizing variants
    my ($chr,$pos,$varid,$ref,$alts) = splice(@data,0,5);
    my @alts = split(/,/,$alts);
    # next 3 are just copied
    my $lineToPrintEnd = "";
    foreach my $i (5..7) {
	$lineToPrintEnd .= "\t".shift(@data);
    }

    # from FORMAT we need GT, AF and DP/DPI/AD
    my $format = shift(@data);
    # GT and AF should always be first, check it
    ($format =~ /^GT:AF:/) || die "E $0: GT:AF: aren't the first FORMAT keys in:\n$line\n";
    # find DP, DPI and/or AD indexes
    my ($dpCol,$dpiCol, $adCol) = (0,0,0);
    my @format = split(/:/, $format);
    foreach my $i (2..$#format) {
	($format[$i] eq "DP") && ($dpCol = $i);
	($format[$i] eq "DPI") && ($dpiCol = $i);
	($format[$i] eq "AD") && ($adCol = $i);
    }
    ($dpCol==0) && ($dpiCol==0) && ($adCol==0) && 
	die "E $0: none of DP, DPI or AD in format:\n$line\n";
    # print new FORMAT
    $lineToPrintEnd .= "\tGENOS";

    # sanity check
    (@data == @samples) || 
	die "E $0: we don't have same number of samples and data columns in:\n$line\n";
    # construct comma-separated list of samples for each called genotype
    my %geno2samples = ();
    # remember seen alleles (for multiallelic cleanup):  $seenAlleles[$i]==1 if
    # allele $i was seen ($i==0 for REF, >=1 for ALTs)
    my @seenAlleles = (0)x(1 + @alts);

    foreach my $i (0..$#data) {
	# add trailing ':' so we know GT is followed by : (eg for NOCALLS)
	$data[$i] .= ':';
	# ignore NOCALLs
	($data[$i] =~ m~^\./\.:~) && next;
	($data[$i] =~ m~^(\d+)/(\d+):([^:]+):~) || 
	    die "E $0: cannot grab genotype and AF for sample $i in $data[$i] in line:\n$line\n";
	my ($geno1,$geno2,$af) = ($1,$2,$3);
	$seenAlleles[$geno1] = 1;
	$seenAlleles[$geno2] = 1;
	my $geno = "$geno1/$geno2";
	if (defined $geno2samples{$geno}) { $geno2samples{$geno} .= ","; }
	else { $geno2samples{$geno} = ""; }

	$geno2samples{$geno} .= $samples[$i];
	if ($af ne '.') {
	    # if we have an AF this is a HET or HV call
	    # find DP/DPI/sumOfADs (whichever is defined and not '.' and biggest)
	    my @thisData = split(/:/, $data[$i]);
	    my $dp = 0;
	    ($dpCol) && ($thisData[$dpCol]) && ($thisData[$dpCol] ne '.') && ($dp = $thisData[$dpCol]);
	    ($dpiCol) && ($thisData[$dpiCol]) && ($thisData[$dpiCol] ne '.') && ($dp < $thisData[$dpiCol]) &&
		($dp = $thisData[$dpiCol]);
	    if (($adCol) && ($thisData[$adCol]) && ($thisData[$adCol] ne '.')) {
		my $sumOfADs = 0;
		foreach my $ad (split(/,/,$thisData[$adCol])) {
		    $sumOfADs += $ad;
		}
		($dp < $sumOfADs) && ($dp = $sumOfADs);
	    }
	    ($dp) || die "E $0: AF is $af but couldn't find DP or DPI or sumOfADs in $data[$i]\n$line\n";
	    # add [DP:AF] after sample ID
	    $geno2samples{$geno} .= "[$dp:$af]";
	}
    }

    # now print data for each genotype, removing data from %geno2samples as we go
    # format is: genotype1~sample1,sample2|genotype2~sample3|genotype3~sample4,sample5

    # $old2new[$old] == $new if allele $old must become $new, 
    #                == -1 if it's never present (for sanity),
    #                undefined otherwise (no change needed)
    my @old2new;
    {
	my $skipped = 0;
	# start at 1 because we always keep REF
	foreach my $i (1..$#seenAlleles) {
	    if (! $seenAlleles[$i]) {
		$old2new[$i] = -1;
		$skipped++;
	    }
	    elsif ($skipped) {
		$old2new[$i] = $i - $skipped;
	    }
	}
    }
    # build new ALTs
    my @newAlts = ();
    foreach my $i (0..$#alts) {
	if ((!defined $old2new[$i+1]) || ($old2new[$i+1] != -1)) {
	    push(@newAlts,$alts[$i]);
	}
    }

    # normalize new ALTs:
    # 1. if length >= 2 for REF and all ALTS, and if REF and all ALTs have 
    #    common ending bases, remove them (keeping at least 1 base everywhere).
    while ($ref =~ /\w(\w)$/) {
	# ref has at least 2 chars
	my $lastRef = $1;
	my $removeLast = 1;
	foreach my $alt (@newAlts) {
	    if ($alt !~ /\w$lastRef$/) {
		# this alt is length one or doesn't end with $lastRef
		$removeLast = 0;
		last;
	    }
	}
	if ($removeLast) {
	    # OK remove last base from REF and all @alts
	    ($ref =~ s/$lastRef$//) || 
		die "E $0: WTF can't remove $lastRef from end of ref $ref\n";
	    foreach my $i (0..$#newAlts) {
		($newAlts[$i] =~ s/$lastRef$//) || 
		    die "E $0: WTF can't remove $lastRef from end of newAlt $i == $newAlts[$i]\n";
	    }
	}
	else {
	    # can't remove $lastRef, get out of while loop
	    last;
	}
    }
    # 2. if length >= 2 for REF and all ALTS, and if REF and all ALTs have 
    #    common starting bases, remove them (keeping at least 1 base everywhere)
    #    and adjust POS.
    while ($ref =~ /^(\w)\w/) {
	# ref has at least 2 chars
	my $firstRef = $1;
	my $removeFirst = 1;
	foreach my $alt (@newAlts) {
	    if ($alt !~ /^$firstRef\w/) {
		# this alt is length one or doesn't start with $firstRef
		$removeFirst = 0;
		last;
	    }
	}
	if ($removeFirst) {
	    # OK remove first base from REF and all @alts
	    ($ref =~ s/^$firstRef//) || 
		die "E $0: WTF can't remove $firstRef from start of ref $ref\n";
	    foreach my $i (0..$#newAlts) {
		($newAlts[$i] =~ s/^$firstRef//) || 
		    die "E $0: WTF can't remove $firstRef from start of alt $i == $newAlts[$i]\n";
	    }
	    # increment POS
	    $pos++;
	}
	else {
	    # can't remove $firstRef, get out of while loop
	    last;
	}
    }

    # build line to print with fixed POS,REF,ALTs
    my $lineToPrint = "$chr\t$pos\t$varid\t$ref\t".join(',',@newAlts)."$lineToPrintEnd";

    # HV
    $lineToPrint .= "\t";
    foreach my $geno (keys %geno2samples) {
	($geno eq '0/0') && next; # skip HR
	($geno =~ m~^(\d+)/\1$~) || next;
	my $allele = $1;
	my $newGeno = $geno;
	if (defined $old2new[$allele]) {
	    ($old2new[$allele] == -1) &&
		die "E $0: allele $allele is present in a geno $geno but it was never seen!\n$line\n";
	    $newGeno = $old2new[$allele]."/".$old2new[$allele];
	}
	$lineToPrint .= "$newGeno~".$geno2samples{$geno}."|";
	delete($geno2samples{$geno});
    }
    # remove last | if needed (not needed if there are no HVs)
    $lineToPrint =~ s/\|$// ;

    # HET == 0/*
    $lineToPrint .= "\t";
    foreach my $geno (keys %geno2samples) {
	($geno eq '0/0') && next; # skip HR
	($geno =~ m~^0/(\d+)$~) || next;
	my $allele = $1;
	my $newGeno = $geno;
	if (defined $old2new[$allele]) {
	    ($old2new[$allele] == -1) &&
		die "E $0: allele $allele is present in a geno $geno but it was never seen!\n$line\n";
	    $newGeno = "0/".$old2new[$allele];
	}
	$lineToPrint .= "$newGeno~".$geno2samples{$geno}."|";
	delete($geno2samples{$geno});
    }
    # remove last | if there was at least one HET geno
    $lineToPrint =~ s/\|$// ;

    # OTHER
    $lineToPrint .= "\t";
    foreach my $geno (sort keys %geno2samples) {
	($geno eq '0/0') && next; # skip HR
	($geno =~ m~^(\d+)/(\d+)$~) || die "E $0: cannot parse geno $geno in line:\n$line\n";
	my ($allele1,$allele2) = ($1,$2);
	my $newGeno = $geno;
	if (defined $old2new[$allele1]) {
	    ($old2new[$allele1] == -1) &&
		die "E $0: allele $allele1 is present in a geno $geno but it was never seen!\n$line\n";
	    $newGeno = $old2new[$allele1]."/$allele2";
	}
	if (defined $old2new[$allele2]) {
	    ($old2new[$allele2] == -1) &&
		die "E $0: allele $allele2 is present in a geno $geno but it was never seen!\n$line\n";
	    ($newGeno =~ s~/$allele2$~/$old2new[$allele2]~) ||
		die "E $0: cannot subst allele2 $allele2 in geno $geno, newGeno $newGeno, line:\n$line\n";
	}
	$lineToPrint .= "$newGeno~".$geno2samples{$geno}."|";
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
warn "I $0: $now - ALL DONE, completed successfully!\n";
