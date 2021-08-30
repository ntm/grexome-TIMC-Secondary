#!/usr/bin/perl

# 25/03/2018
# NTM

# Take in stdin a VCF file, print to stdout a TSV file with
# some of the VCF info.
# We need CSQ to be present, we're assuming the VCF went through VEP.
#
# We define a new MODHIGH impact, as follows:
# - missense variants are upgraded from MODER to MODHIGH if they are
#   considered deleterious by most methods (details in code, look for "missense");
# - similarly, splice_region_variants are upgraded from LOW to MODHIGH if they
#   are considered deleterious (look for "splice" in code).
#
# When a variant impacts several transcripts and/or for multi-allelic lines
# (ie with several variants), we will print one line per transcript and per 
# variant. This results in some redundancy but makes it much easier to filter
# things and find the best candidates.
# For multi-allelic positions, the GENOS columns are modified: 
# HV or HET genotypes~samples for alleles other than ALLELE_NUM are moved to OTHER.

use strict;
use warnings;
use File::Basename qw(basename);
use POSIX qw(strftime);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";

# stuff we grab from INFO:
# CSQ==VEP, but we only want some VEP fields

# @goodVeps: VEP fields we want, in the order we want them printed.
# Some of these are hard-coded in the "missense" and "splice"
# upgrade-to-MODHIGH code, make sure they get fixed there as
# well if they change names.
# Current available annotations (26/08/2021) are:
# Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|
# HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|
# Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|
# HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|RefSeq|SIFT|PolyPhen|HGVS_OFFSET|
# AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|
# gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|
# CLIN_SIG|SOMATIC|PHENO|PUBMED|
# CADD_raw_rankscore|MetaRNN_pred|MetaRNN_rankscore|MutationTaster_pred|REVEL_rankscore|
# ada_score|rf_score
my @goodVeps = ("SYMBOL","Gene","IMPACT","Consequence","Feature","CANONICAL",
		"BIOTYPE","VARIANT_CLASS","RefSeq","MANE_SELECT","MANE_PLUS_CLINICAL",
		"ALLELE_NUM","EXON","INTRON",
		"HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position",
		"SIFT","PolyPhen",
		"MetaRNN_pred","MetaRNN_rankscore","CADD_raw_rankscore",
		"MutationTaster_pred","REVEL_rankscore",
		"ada_score","rf_score",
		"gnomAD_AF","AF",
		"Existing_variation","CLIN_SIG","SOMATIC","PHENO","PUBMED");

# VCF headers: ignore them but grab the:
# - VEP CSQ field names and store them in @vepNames;
# - data column names and store them as a string in $dataHeaders.
my @vepNames;
my $dataHeaders;
while (my $line = <STDIN>) {
    chomp($line);
    ($line =~ /^#/) ||
	die "E $0: parsing header but found non-header line:\n$line\n" ;

    # VEP CSQ header
    if ($line =~ /^##INFO=<ID=CSQ/) {
	($line =~ /Format: ([^"]+)">$/) ||
	    die "E $0: found CSQ header line but can't find Format:\n$line\n" ;
	@vepNames = split(/\|/,$1);
    }
    #CHROM == last header line
    elsif ($line =~ /^#CHROM/) {
	my @headers = split(/\t/,$line);
	(@headers == 13) || die "E $0: CHROM header line has wrong number of fields\n";
	$dataHeaders = join("\t",@headers[9..12]);
	last;
    }
}
(@vepNames) || 
    die "E $0: done parsing headers but vepNames still empty!\n" ;

# sanity: make sure all my @goodVeps fields exist
{
    my %vepNames;
    foreach my $v (@vepNames) {
	$vepNames{$v} = 1;
    }
    foreach my $v (@goodVeps) {
	($vepNames{$v} == 1) ||
	    die "E $0: the VEP field $v doesn't seem to exist anymore, update goodVeps and maybe other things!\n";
    }
}

# Make our own headers for the TSV
my $headerTsv = "POSITION\tREF\tALT";
# POSITION will be chrom:pos

# selected VEP fields will be printed in @goodVeps order
$headerTsv .= "\t".join("\t",@goodVeps);
# finally we print the sample IDs for each GENO, just copied from infile
$headerTsv .= "\t$dataHeaders";
# sanity: below we assume GENOS are HV,HET,OTHER,HR in that order
($dataHeaders eq "HV\tHET\tOTHER\tHR") ||
    die "E $0: GENOS are not in expected order but the code depends on it, dataHeaders is: $dataHeaders\n";
# done, print header
$headerTsv .= "\n" ;
print $headerTsv ;


# Deal with data
while (my $line =<STDIN>) {
    chomp($line);

    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@dataCols) = split(/\t/,$line,-1);

    # VEP removed trailing empty fields (did it? nvm), fill them if needed
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
	    die "E $0: wrong number of fields in CSQ (full line is below): $thisCsq\n$line\n" ;
	foreach my $i (0..$#vepNames) {
	    # replace all slashes by backslashes, for excel :(
	    $csqTmp[$i] =~ s~/~\\~g ;
	    $thisCsq{$vepNames[$i]} = $csqTmp[$i] ;
	}

	# upgrade splice_region_variant from LOW to MODHIGH if:
	# (ada_score > 0.6) AND (rf_score > 0.6)
	# TODO: CADD v1.6 (AKA CADD-Splice) is supposedly good for splice-affecting variants,
	# I want to use CADD_raw_rankscore here. Unfortunately currently the dbNSFP VEP
	# plugin doesn't report CADD_raw_rankscore for splicing variants, despite the score
	# being present in the dbNSFP gz file
	if (($thisCsq{"IMPACT"} eq "LOW") && ($thisCsq{"Consequence"}) &&
	    ($thisCsq{"Consequence"} =~ /splice_region_variant/)) {
	    # NOTE: VEP consequence column can have several &-separated consequences,
	    # we just want splice_region_variant to be present somewhere
	    if (($thisCsq{"ada_score"}) && ($thisCsq{"ada_score"} > 0.6) &&
		($thisCsq{"rf_score"}) && ($thisCsq{"rf_score"} > 0.6)) {
		$thisCsq{"IMPACT"} = "MODHIGH";
	    }
	}

	# upgrade putatively deleterious missense variants:
	if (($thisCsq{"IMPACT"} eq "MODERATE") && ($thisCsq{"Consequence"}) &&
	    ($thisCsq{"Consequence"} =~ /missense_variant/)) {
	    # upgrade to MODHIGH if at least $minPassedFrac criteria are passed,
	    # among the following:
	    # - SIFT -> deleterious
	    # - Polyphen -> probably_damaging
	    # - CADD_raw_rankscore >= 0.7
	    # - MutationTaster_pred -> contains at least one A or D
	    # - REVEL_rankscore >= 0.7
	    # - MetaRNN_pred -> contains at least one D(amaging)
	    # 
	    # NOTE: sometimes we don't have any prediction for some predictors,
	    # due to dbNSFP using older transcripts and/or bugs in VEP or VEP plugins...
	    # $minPassedFrac allows to upgrade variants in those cases
	    my $minPassedFrac = 0.6;
	    my $passed = 0;
	    my $totalPreds = 0;
	    if ($thisCsq{"SIFT"}) {
		$totalPreds++;
		# deleterious\( so we don't get deleterious_low_confidence
		($thisCsq{"SIFT"} =~ /deleterious\(/) && ($passed++);
	    }
	    if ($thisCsq{"PolyPhen"}) {
		$totalPreds++;
		($thisCsq{"PolyPhen"} =~ /probably_damaging/) && ($passed++);
	    }
	    if ($thisCsq{"CADD_raw_rankscore"}) {
		$totalPreds++;
		($thisCsq{"CADD_raw_rankscore"} >= 0.7) && ($passed++);
	    }
	    if ($thisCsq{"MutationTaster_pred"}) {
		$totalPreds++;
		($thisCsq{"MutationTaster_pred"} =~ /[AD]/) && ($passed++);
	    }
	    if ($thisCsq{"REVEL_rankscore"}) {
		$totalPreds++;
		($thisCsq{"REVEL_rankscore"} >= 0.7) && ($passed++);
	    }
	    if ($thisCsq{"MetaRNN_pred"}) {
		$totalPreds++;
		($thisCsq{"MetaRNN_pred"} =~ /D/) && ($passed++);
	    }

	    if (($totalPreds) && ($passed / $totalPreds >= $minPassedFrac)) {
		# $totalPreds==0 can only happen due to bugs in VEP, but it does
		$thisCsq{"IMPACT"} = "MODHIGH";
	    }
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

	if ($alt !~ /,/) {
	    # mono-allelic, print immediately
	    print "$toPrintStart$vepToPrint\t".join("\t",@dataCols)."\n";
	}
	else {
	    # multi-allelic, have to fix the GENOS
	    my $allele = $thisCsq{"ALLELE_NUM"};
	    ($allele =~ /^\d+$/) || 
		die "E $0: ALLELE_NUM not a number or undefined? CSQ is $thisCsq\n";
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

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";
