#!/usr/bin/perl

# 25/03/2018
# NTM

# Take in stdin a VCF file, print to stdout a TSV file with
# some of the VCF info.
# We need CSQ to be present, we're assuming the VCF went through VEP.
#
# When a variant impacts several transcripts and/or for multi-allelic lines
# (ie with several variants), we will print one line per transcript and per 
# variant. This results in some redundancy but makes it much easier to filter
# things and find the best candidates.
# For multi-allelic positions, the GENOS columns are modified: 
# HV or HET geotypes~samples for alleles other than ALLELE_NUM are moved to OTHER.

use strict;
use warnings;

# stuff we grab from INFO:
# CSQ==VEP, but we only want some VEP fields

# @goodVeps: VEP fields we want, in the order we want them printed.
# Current available annotations (11/07/2019) are:
# Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|
# HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|
# Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|
# HGNC_ID|CANONICAL|RefSeq|GENE_PHENO|SIFT|PolyPhen|
# HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|
# gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|
# gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|PUBMED|
# CADD_raw_rankscore|MutationTaster_pred|REVEL_rankscore|
# ada_score|rf_score
my @goodVeps = ("SYMBOL","Gene","IMPACT","Consequence","Feature","PICK","CANONICAL",
		"RefSeq","BIOTYPE","ALLELE_NUM","EXON","INTRON",
		"HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position",
		"SIFT","PolyPhen",
		"CADD_raw_rankscore","MutationTaster_pred","REVEL_rankscore",
		"ada_score","rf_score",
		"gnomAD_AF","AF","AA_AF","EA_AF",
		"Existing_variation","CLIN_SIG","SOMATIC","PHENO",
		"GENE_PHENO","PUBMED");

# VCF headers: ignore them but grab the:
# - VEP CSQ field names and store them in @vepNames;
# - data column names and store them as a string in $dataHeaders.
my @vepNames;
my $dataHeaders;
while (my $line = <STDIN>) {
    chomp($line);
    ($line =~ /^#/) ||
	die "E: parsing header but found non-header line:\n$line\n" ;

    # VEP CSQ header
    if ($line =~ /^##INFO=<ID=CSQ/) {
	($line =~ /Format: ([^"]+)">$/) ||
	    die "E: found CSQ header line but can't find Format:\n$line\n" ;
	@vepNames = split(/\|/,$1);
    }
    #CHROM == last header line
    elsif ($line =~ /^#CHROM/) {
	my @headers = split(/\t/,$line);
	(@headers == 13) || die "E: CHROM header line has wrong number of fields\n";
	$dataHeaders = join("\t",@headers[9..12]);
	last;
    }
}
(@vepNames) || 
    die "E: done parsing headers but vepNames still empty!\n" ;

# Make our own headers for the TSV
my $headerTsv = "POSITION\tREF\tALT";
# POSITION will be chrom:pos

# selected VEP fields will be printed in @goodVeps order
$headerTsv .= "\t".join("\t",@goodVeps);
# finally we print the sample IDs for each GENO, just copied from infile
$headerTsv .= "\t$dataHeaders";
# sanity: below we assume GENOS are HV,HET,OTHER,HR in that order
($dataHeaders eq "HV\tHET\tOTHER\tHR") ||
    die "E: GENOS are not in expected order but the code depends on it, dataHeaders is: $dataHeaders\n";
# done, print header
$headerTsv .= "\n" ;
print $headerTsv ;


# Deal with data
while (my $line =<STDIN>) {
    chomp($line);

    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@dataCols) = split(/\t/,$line);

    # VEP removed trailing empty fields, fill them if needed
    while (@dataCols < 4) {
	push(@dataCols, "");
    }

    # line to print
    my $toPrintStart = "$chr:$pos\t$ref\t$alt";

    # grab CSQs from $info
    my @csqs;
    foreach my $thisInfo (split(/;/,$info)) {
	# we only want the CSQ=
	($thisInfo =~ /^CSQ=(.+)$/) || next;
	my $csqs = $1;
	@csqs = split(/,/, $csqs) ;
    }

    foreach my $thisCsq (@csqs) {
	# third arg to split, -1, to generate elements even if last ones are empty
	my @csqTmp = split(/\|/, $thisCsq, -1) ;
	# store in hash: key is VEP key, value is value in this CSQ
	my %thisCsq;
	(@csqTmp == @vepNames) || 
	    die "wrong number of fields in CSQ (full line is below): $thisCsq\n$line\n" ;
	foreach my $i (0..$#vepNames) {
	    # replace all slashes by backslashes, for excel :(
	    $csqTmp[$i] =~ s~/~\\~g ;
	    $thisCsq{$vepNames[$i]} = $csqTmp[$i] ;
	}
	# HGVSp has %3D for = in synonymous variants, substitute
	($thisCsq{"HGVSp"}) && ($thisCsq{"HGVSp"} =~ s/%3D$/=/);
	# build vepToPrint string based on @goodVeps
	my $vepToPrint = "";
	foreach my $vepName (@goodVeps) {
	    if (! defined $thisCsq{$vepName}) {
		$vepToPrint .= "\t";
	    }
	    else {
		$vepToPrint .= "\t".$thisCsq{$vepName};
	    }
	}

	if ($alt !~ /,/) {
	    # mono-allelic, print immediately
	    print "$toPrintStart$vepToPrint\t".join("\t",@dataCols)."\n";
	}
	else {
	    # multi-allelic, have to fix the GENOS
	    my $allele = $thisCsq{"ALLELE_NUM"};
	    ($allele =~ /^\d+$/) || 
		die "E: ALLELE_NUM not a number or undefined? CSQ is $thisCsq\n";
	    my ($hvs,$hets,$others,$hrs) = @dataCols;
	    my @hvs = split(/\|/,$hvs);
	    my @hets = split(/\|/,$hets);
	    my @others = split(/\|/,$others);

	    my ($goodHv,$goodHet) = ("","");
	    # move HVs and HETs for other alleles to OTHER
	    foreach my $hv (@hvs) {
		if ($hv =~ /^$allele\/$allele~/) {
		    $goodHv = $hv;
		}
		else {
		    push(@others,$hv);
		}
	    }
	    foreach my $het (@hets) {
		if ($het =~ /^0\/$allele~/) {
		    $goodHet = $het;
		}
		else {
		    push(@others,$het);
		}
	    }
	    my $fixedData = "$goodHv\t$goodHet\t".join('|',@others)."\t$hrs";
	    print "$toPrintStart$vepToPrint\t$fixedData\n";
	}
    }
}
