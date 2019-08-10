#!/usr/bin/perl

# 25/03/2018
# NTM

# Take in stdin a VCF file, print to stdout a TSV file with
# some of the VCF info.
# We need some fields to be present (COUNT_*, CSQ), I am 
# assuming the VCF went through previous steps 3-5.
# We also need HV and HET count attribute names to be alphabetically
# smaller than the OTHER attribute name for the same cohort. This
# is crudely checked (just keep names as they are and we are fine).
#
# The COUNTs become numbers: the actual genotypes are stripped,
# but for multiallelic sites we only count in HV and HET the
# samples genotyped with the current allele_num; samples
# HV or HET with a different ALT allele are counted in OTHER.
# The actual genotypes will still be in the data columns.
# When a variant impacts several transcripts, we will print
# one line per transcript. Many fields will be identical between
# these lines, but many will be different... This should make
# it much easier to filter things and find the best candidates.

use strict;
use warnings;

# stuff we grab from INFO:
# COUNT_* and CSQ==VEP, but we only want some VEP fields

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


# store the goodInfo keys in hash, key==good key, value==1
# the COUNT_* keys will be grabbed from the VCF headers
my %goodInfo = ("CSQ"=>1);

# VCF headers: ignore them but grab the:
# - COUNT_* names and store them in %goodInfo;
# - VEP CSQ field names and store them in @vepNames;
# - data column names and store them as a string in $dataHeaders.
my @vepNames;
my $dataHeaders;
while (my $line = <STDIN>) {
    chomp($line);
    ($line =~ /^#/) ||
	die "in step 3, parsing header but found non-header line:\n$line\n" ;

    # COUNT_* definitions
    if ($line =~ /^##INFO=<ID=(COUNT_[^,]+),/) {
	$goodInfo{$1} = 1;
    }
    # VEP CSQ header
    elsif ($line =~ /^##INFO=<ID=CSQ/) {
	($line =~ /Format: ([^"]+)">$/) ||
	    die "found CSQ header line but can't find Format:\n$line\n" ;
	@vepNames = split(/\|/,$1);
    }
    #CHROM == last header line
    elsif ($line =~ /^#CHROM/) {
	my @headers = split(/\t/,$line);
	(@headers == 13) || die "CHROM header line has wrong number of fields\n";
	$dataHeaders = join("\t",@headers[9..12]);
	last;
    }
}
(@vepNames) || 
    die "in step 3, done parsing headers but vepNames still empty!\n" ;

# sanity check: we need the attribute names that store the counts for HVs 
# and HETs to be alphabetically smaller than the attribute of OTHER.
# I also hard-coded *HET, *HV and *OTHER below...
# I'm just testing for one cohort that the attribute names didn't change.
((defined $goodInfo{"COUNT_Astheno_HV"}) && (defined $goodInfo{"COUNT_Astheno_HET"}) &&
    (defined $goodInfo{"COUNT_Astheno_OTHER"})) ||
    die "E: the HV, HET and OTHER attr names are HARD-CODED, just don't change them please... or fix the code\n";

# Make our own headers for the TSV
my $headerTsv = "POSITION\tREF\tALT";
# POSITION will be chrom:pos

# COUNT_* fields will be printed in sort order
foreach my $key (sort keys %goodInfo) {
    ($key eq "CSQ") && next;
    $headerTsv .= "\t$key";
}
# and selected VEP fields will be printed in @goodVeps order
$headerTsv .= "\t".join("\t",@goodVeps);
# finally we print the sample IDs for each GENO, just copied from infile
$headerTsv .= "\t$dataHeaders";
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

    # parse %goodInfo fields from $info into %data
    my %data = () ;
    foreach my $thisInfo (split(/;/,$info)) {
	# some fields are just flags (no =), skip them
	($thisInfo =~ /^([^=]+)=(.+)$/) || next;
	my ($field,$value) = ($1,$2);
	(defined $goodInfo{$field}) && ($data{$field} = $value) ;
    }

    # prepare the stuff that will be printed at the end of each line, ie the GENOs
    my $toPrintEnd = "\t".join("\t",@dataCols);

    # split the VEP CSQ field
    my @csqs = split(/,/, $data{"CSQ"}) ;

    foreach my $thisCsq (@csqs) {
	# third arg to split, -1, to generate elements even if last ones are empty
	my @csqTmp = split(/\|/, $thisCsq, -1) ;
	# store in hash: key is VEP key, value is value in this CSQ
	my %thisCsq;
	(@csqTmp == @vepNames) || 
	    die "wrong number of fields in CSQ (full line is below): $thisCsq\n$line\n" ;
	foreach my $i (0..$#vepNames) {
	    # special case: for refseq genes we don't want Gene, it's the HGNC 
	    # id (a number) but we don't want it
	    ($vepNames[$i] eq "Gene") && ($csqTmp[$i] =~ /^\d+$/) && next;
	    # replace all slashes by backslashes, for excel :(
	    $csqTmp[$i] =~ s~/~\\~g ;
	    $thisCsq{$vepNames[$i]} = $csqTmp[$i] ;
	}
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

	# need alleleNum to sum the COUNTs correctly
	my $alleleNum = $thisCsq{"ALLELE_NUM"};
	($alleleNum) || die "E: no ALLELE_NUM value in CSQ $thisCsq\nline is\n$line\n";
	# COUNTs to print, printing in sort order (same as for $headerTsv)
	my $countsToPrint = "";
	# number of HVs and HETs that concern another alleleNum and must 
	# count as OTHER in this line
	my $other = 0;
	foreach my $key (sort keys %goodInfo) {
	    ($key eq "CSQ") && next;
	    my $total = 0;
	    if (($key =~ /_HV$/) || ($key =~ /_HET$/)) {
		if (defined $data{$key}) {
		    foreach my $count (split(/,/, $data{$key})) {
			($count =~ /^0\/(\d+):(\d+)$/) || 
			    ($count =~ /^(\d+)\/\1:(\d+)$/) || 
			    die "cannot decompose HV/HET count: $count\n";
			my ($allele,$count) = ($1,$2);
			if ($allele == $alleleNum) {
			    $total += $count;
			}
			else {
			    $other += $count;
			}
		    }
		}
	    }
	    elsif ($key =~ /_OTHER$/) {
		# use accumulator and reset
		$total += $other;
		$other = 0;
		if (defined $data{$key}) {
		    # add actual OTHER counts
		    foreach my $count (split(/,/, $data{$key})) {
			($count =~ /:(\d+)$/) || 
			    die "cannot grab count from OTHER: $count\n";
			$total += $1;
		    }
		}
	    }
	    elsif (defined $data{$key}) {
		# HR
		($data{$key} =~ /^0\/0:(\d+)$/) ||
		    die "cannot grab count for HR: $data{$key}\n";
		$total = $1;
	    }
	    $countsToPrint .= "\t$total";
	}
	
	# done with this csq, print full line
	print "$toPrintStart$countsToPrint$vepToPrint$toPrintEnd\n";
    }
}

