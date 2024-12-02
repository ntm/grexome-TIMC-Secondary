#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2024
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


# 25/03/2018
# NTM

# Take in stdin a VCF file, print to stdout a TSV file with
# some of the VCF info.
# We need CSQ to be present, we're assuming the VCF went through VEP.
#
# We define a new MODHIGH impact, as follows:
# - missense variants are upgraded from MODER to MODHIGH if they are
#   considered deleterious by most methods (details in code, look for "missense");
# - variants predicted to affect splicing are upgraded to MODHIGH if they
#   are considered splice-affecting by most methods (currently by ada_score
#   and rf_score from dbscSNV, SpliceAI, and CADD-Splice, look for "splice").
#
# Also, in order to study miRNAs and other ncRNAs, some IMPACTs are
# upgraded from MODIFIER to LOW (look for "miRNA").
#
# Also, in an attempt to make use of VEP's --regulatory switch, variants
# affecting TFBS's and/or severely affecting regulatory regions are upgraded
# from MODIFIER to LOW or MODHIGH (look for "TF") [ACTUALLY NOT DOING THIS YET,
# --regulatory is too noisy for now].
#
# Overall for CNVs, VEP currently can produce 4 different consequences,
# but they are now (release 110, 16/10/2023) all HIGH.
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

# @goodVeps: VEP fields we want, in the order we want them printed (exception for
# SpliceAI_* fields, which we will replace by a single SpliceAI_DS containing
# the max delta score).
# Some of these are hard-coded in the "missense" and "splice" upgrade-to-MODHIGH
# code, make sure they get fixed there as well if they change names.
# Current available annotations produced by runVEP.pl (21/10/2024) are:
# Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|
# HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|
# Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|
# HGNC_ID|CANONICAL|MANE|MANE_SELECT|MANE_PLUS_CLINICAL|SIFT|PolyPhen|miRNA|HGVS_OFFSET|
# AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|
# gnomADe_FIN_AF|gnomADe_MID_AF|gnomADe_NFE_AF|gnomADe_REMAINING_AF|gnomADe_SAS_AF|
# gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|
# gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_REMAINING_AF|gnomADg_SAS_AF|
# CLIN_SIG|SOMATIC|PHENO|
# CADD_PHRED|CADD_RAW|ALFA_Total_AF|
# CADD_raw_rankscore|MetaRNN_pred|MetaRNN_rankscore|MutationTaster_pred|REVEL_rankscore|
# ada_score|rf_score|
# SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|
# SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL|
# SpliceAI_pred_SYMBOL|am_class|am_pathogenicity|pLI_gene_value
my @goodVeps = ("SYMBOL","Gene","IMPACT","Consequence","Feature","CANONICAL",
		"BIOTYPE","MANE_SELECT","pLI_gene_value","Existing_variation",
		"HGVSc","HGVSp","ALLELE_NUM","EXON","INTRON",
		"cDNA_position","CDS_position","Protein_position",
		"miRNA","VARIANT_CLASS",
		"gnomADe_AF","gnomADg_AF","AF","ALFA_Total_AF",
		"SIFT","PolyPhen","am_class","am_pathogenicity",
		"MetaRNN_pred","MetaRNN_rankscore","CADD_raw_rankscore",
		"MutationTaster_pred","REVEL_rankscore",
		"ada_score","rf_score","CADD_PHRED",
		"SpliceAI_pred_DS_AG","SpliceAI_pred_DS_AL","SpliceAI_pred_DS_DG",
		"SpliceAI_pred_DS_DL");

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

# With non-human data some @goodVeps don't exist, purge them out. Also replace
# the 4 SpliceAI* scores with a single SpliceAI_DS
{
    my %vepNames;
    foreach my $v (@vepNames) {
	$vepNames{$v} = 1;
    }
    my @missingVeps = ();
    my @vepsFound = ();
    foreach my $v (@goodVeps) {
	if ($vepNames{$v}) {
	    if ($v eq "SpliceAI_pred_DS_AG") {
		push(@vepsFound, "SpliceAI_DS");
	    }
	    elsif ($v =~ /^SpliceAI_pred/) {
		# ignore, we only want one SpliceAI_DS column
	    }
	    else {
		push(@vepsFound, $v);
	    }
	}
	else {
	    push(@missingVeps, $v);
	}
    }
    if (@missingVeps) {
	warn "W $0: the following VEP fields are missing: " . join(' ', @missingVeps) . "\n";
	warn "W $0: with non-human data missing fields are expected (eg fields from human-only plugins),\n";
	warn "W $0: but with human data this suggests that the code needs to be updated (open a github issue),\n";
	warn "W $0: or that your installation of VEP is outdated or missing some plugins / data\n";
    }
    @goodVeps = @vepsFound;
}

# Make our own headers for the TSV
my $headerTsv = "POSITION\tREF\tALT";
# POSITION will be chrom:pos

# selected VEP fields will be printed in @goodVeps order
$headerTsv .= "\t".join("\t", @goodVeps);
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

    # for CNVs, append the END= value to $alt
    if (($alt eq '<DUP>') || ($alt eq '<DEL>')) {
	foreach my $thisInfo (split(/;/,$info)) {
	    ($thisInfo =~ /^END=(\d+)$/) || next;
	    # remove trailing '>', we'll add it back after $end
	    chop($alt);
	    $alt .= ":$1>";
	    last;
	}
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

	# create SpliceAI_DS with max of the 4 SpliceAI* scores
	if ($thisCsq{"SpliceAI_pred_DS_AG"} ne "") {
	    $thisCsq{"SpliceAI_DS"} = $thisCsq{"SpliceAI_pred_DS_AG"};
	    foreach my $f ("SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL") {
		($thisCsq{$f} > $thisCsq{"SpliceAI_DS"}) && ($thisCsq{"SpliceAI_DS"} = $thisCsq{$f});
	    }
	}

	# upgrade putatively deleterious missense variants:
	if (($thisCsq{"IMPACT"} eq "MODERATE") && ($thisCsq{"Consequence"}) &&
	    ($thisCsq{"Consequence"} =~ /missense_variant/)) {
	    # upgrade to MODHIGH if more than $minPassedFracMissense criteria are
	    # passed among the following:
	    # - SIFT -> deleterious
	    # - Polyphen -> probably_damaging
	    # - CADD_raw_rankscore >= 0.7
	    # - MutationTaster_pred -> contains at least one A or D
	    # - REVEL_rankscore >= 0.7
	    # - MetaRNN_pred -> contains at least one D(amaging)
	    # - AlphaMissense class -> 'likely_pathogenic'
	    # 
	    # NOTE: sometimes we don't have any prediction for some predictors,
	    # e.g. for non-human data, or due to dbNSFP using older transcripts
	    # and/or bugs in VEP or VEP plugins...
	    # -> $minPassedFracMissense allows to still upgrade variants, but we also
	    # require that at least 2 predictors are available
	    my $minPassedFracMissense = 0.5;
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
	    if ($thisCsq{"am_class"}) {
		$totalPreds++;
		($thisCsq{"am_class"} =~ /likely_pathogenic/) && ($passed++);
	    }

	    # require at least 2 predictors
	    if (($totalPreds>1) && ($passed / $totalPreds > $minPassedFracMissense)) {
		$thisCsq{"IMPACT"} = "MODHIGH";
	    }
	}

	# upgrade any variant to MODHIGH if it putatively alters splicing:
	if (($thisCsq{"IMPACT"} ne "HIGH") && ($thisCsq{"IMPACT"} ne "MODHIGH")) {
	    # upgrade to MODHIGH if at least $minPassedFracSplicing criteria (and at 
	    # least 2) are passed, among the following:
	    # - SpliceAI max(DS_AG, DS_AL, DS_DG, DS_DL) > cutoff
	    # - CADD-PHRED > cutoff (presumably via CADD-Splice AKA CADD 1.6) 
	    # - ada_score > 0.6 (dbscSNV author-recommended cutoff)
	    # - rf_score > 0.6 (dbscSNV author-recommended cutoff)

	    # cutoffs:
	    # SpliceAI authors recommend 0.5 and say 0.8 is high-precision
	    my $spliceAI_cutoff = 0.5;
	    # for CADD_PHRED, 20-30 seems reasonable
	    my $cadd_cutoff = 20;

	    # non-human data and some variants for human data don't have scores -> use a frac
	    my $minPassedFracSplicing = 0.5;
	    my $passed = 0;
	    my $totalPreds = 0;
	    # SpliceAI
	    if (defined $thisCsq{"SpliceAI_DS"}) {
		$totalPreds++;
		($thisCsq{"SpliceAI_DS"} > $spliceAI_cutoff) && ($passed++);
	    }
	    # CADD-Splice
	    if ($thisCsq{"CADD_PHRED"} ne "") {
		$totalPreds++;
		($thisCsq{"CADD_PHRED"} > $cadd_cutoff) && ($passed++);
	    }
	    # dbscSNV
	    if ($thisCsq{"ada_score"} ne "") {
		$totalPreds++;
		($thisCsq{"ada_score"} > 0.6) && ($passed++);
	    }
	    if ($thisCsq{"rf_score"} ne "") {
		$totalPreds++;
		($thisCsq{"rf_score"} > 0.6) && ($passed++);
	    }

	    # require at least 2 predictors
	    if (($totalPreds > 1) && ($passed / $totalPreds > $minPassedFracSplicing)) {
		$thisCsq{"IMPACT"} = "MODHIGH";
	    }
	}

	if (($thisCsq{"IMPACT"} eq "MODIFIER") && ($thisCsq{"Consequence"})) {
	    # upgrade some miRNA and other ncRNA IMPACTs to LOW:
	    # non_coding_transcript_variant (except deep intronic ones), 
	    # non_coding_transcript_exon_variant, and mature_miRNA_variant
	    if ($thisCsq{"Consequence"} =~ /non_coding_transcript_exon_variant|mature_miRNA_variant/) {
		$thisCsq{"IMPACT"} = "LOW";
	    }
	    elsif (($thisCsq{"Consequence"} =~ /non_coding_transcript_variant/) &&
		   ($thisCsq{"Consequence"} ne 'intron_variant&non_coding_transcript_variant')) {
		# deepish intronic variants affecting ncRNAs should stay MODIFIER, but we do want
		# to upgrade eg 'splice_donor_region_variant&intron_variant&non_coding_transcript_variant'
		$thisCsq{"IMPACT"} = "LOW";
	    }

	    ###############
	    # tried upgrading some TFBS and regulatory region variants but it's
	    # very noisy and doesn't seem usable (yet) -> commenting out for now, 
	    # keeping the code because I expect --regulatory to become useful some day...
	    #
	    # upgrade severe TFBS and regulatory region variants to LOW or MODHIGH:
	    # TFBS_ablation and regulatory_region_ablation => MODHIGH 
	    # TF_binding_site_variant with [HIGH_INF_POS==Y] => MODHIGH
	    # TF_binding_site_variant without HIGH_INF_POS => LOW
	    # TFBS_amplification => LOW
	    # (regulatory_region_variant and regulatory_region_amplification stay MODIFIER)
	    # elsif ($thisCsq{"Consequence"} =~ /TFBS_ablation|regulatory_region_ablation/) {
	    # 	$thisCsq{"IMPACT"} = "MODHIGH";
	    # }
	    # elsif (($thisCsq{"Consequence"} =~ /TF_binding_site_variant/) && 
	    # 	   ($thisCsq{"HIGH_INF_POS"}) && ($thisCsq{"HIGH_INF_POS"} eq "Y")) {
	    # 	$thisCsq{"IMPACT"} = "MODHIGH";
	    # }
	    # elsif ($thisCsq{"Consequence"} =~ /TF_binding_site_variant/) {
	    # 	# no HIGH_INF_POS
	    # 	$thisCsq{"IMPACT"} = "LOW";
	    # }
	    # elsif ($thisCsq{"Consequence"} =~ /TFBS_amplification/) {
	    # 	$thisCsq{"IMPACT"} = "LOW";
	    # }
	    ###############
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
