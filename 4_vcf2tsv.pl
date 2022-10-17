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
# - variants predicted to affect splicing are upgraded to MODHIGH if they
#   are considered splice-affecting by most methods (currently by ada_score
#   and rf_score from dbscSNV, SpliceAI, and CADD-Splice, look for "splice");
# - feature_elongation variants are upgraded from MODIFIER to MODHIGH.
#
# In addition, feature_truncation variants are upgraded from MODIFIER to HIGH.
#
# Also, in order to study miRNAs and oter ncRNAs, some IMPACTs are
# upgraded from MODIFIER to LOW (look for "miRNA").
#
# Also, in an attempt to make use of VEP's --regulatory switch, variants
# affecting TFBS's and/or severely affecting regulatory regions are upgraded
# from MODIFIER to LOW or MODHIGH (look for "TF").
#
# Overall for CNVs, VEP currently can produce 4 different consequences:
# - transcript_ablation = DEL covering the complete transcript, HIGH;
# - feature_truncation = DEL with partial coverage of the transcript, 
#   MODIFER for VEP, we upgrade to HIGH;
# - transcript_amplification = DUP covering the complete transcript, HIGH;
# - feature_elongation = DUP with partial coverage, MODIFER upgraded to MODHIGH.
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
# Current available annotations produced by runVEP.pl (17/07/2022) are:
# Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|
# HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|
# Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|
# HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|RefSeq|SIFT|PolyPhen|HGVS_OFFSET|
# AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|
# gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|
# gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|
# gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|
# CLIN_SIG|SOMATIC|PHENO|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|
# CADD_PHRED|CADD_RAW|CADD_raw_rankscore|MetaRNN_pred|MetaRNN_rankscore|MutationTaster_pred|
# REVEL_rankscore|ada_score|rf_score|
# SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|SpliceAI_pred_DS_AG|
# SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL|SpliceAI_pred_SYMBOL

my @goodVeps = ("SYMBOL","Gene","IMPACT","Consequence","Feature","CANONICAL",
		"BIOTYPE","VARIANT_CLASS","RefSeq","MANE_SELECT","MANE_PLUS_CLINICAL",
		"ALLELE_NUM","EXON","INTRON",
		"HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position",
		"SIFT","PolyPhen",
		"MetaRNN_pred","MetaRNN_rankscore","CADD_raw_rankscore",
		"MutationTaster_pred","REVEL_rankscore",
		"ada_score","rf_score","CADD_PHRED",
		"SpliceAI_pred_DS_AG","SpliceAI_pred_DP_AG","SpliceAI_pred_DS_AL",
		"SpliceAI_pred_DP_AL","SpliceAI_pred_DS_DG","SpliceAI_pred_DP_DG",
		"SpliceAI_pred_DS_DL","SpliceAI_pred_DP_DL",
		"gnomADe_AF","gnomADg_AF","AF",
		"MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS",
		"Existing_variation");

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
	((defined $vepNames{$v}) && ($vepNames{$v} == 1)) ||
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

	# upgrade putatively deleterious missense variants:
	if (($thisCsq{"IMPACT"} eq "MODERATE") && ($thisCsq{"Consequence"}) &&
	    ($thisCsq{"Consequence"} =~ /missense_variant/)) {
	    # upgrade to MODHIGH if at least $minPassedFracMissense criteria are
	    # passed among the following:
	    # - SIFT -> deleterious
	    # - Polyphen -> probably_damaging
	    # - CADD_raw_rankscore >= 0.7
	    # - MutationTaster_pred -> contains at least one A or D
	    # - REVEL_rankscore >= 0.7
	    # - MetaRNN_pred -> contains at least one D(amaging)
	    # 
	    # NOTE: sometimes we don't have any prediction for some predictors,
	    # due to dbNSFP using older transcripts and/or bugs in VEP or VEP plugins...
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

	    # require at least 2 predictors
	    if (($totalPreds>1) && ($passed / $totalPreds > $minPassedFracMissense)) {
		$thisCsq{"IMPACT"} = "MODHIGH";
	    }
	}

	# upgrade CNVs with only partial coverage of the transcript:
	if (($thisCsq{"IMPACT"} ne "HIGH") && ($thisCsq{"Consequence"}) &&
	    ($thisCsq{"Consequence"} =~ /feature_truncation/)) {
	    # deletion of one or more exons, upgrade to HIGH
	    $thisCsq{"IMPACT"} = "HIGH";
	}
	elsif (($thisCsq{"IMPACT"} ne "HIGH") && ($thisCsq{"IMPACT"} ne "MODHIGH") &&
	       ($thisCsq{"Consequence"}) && ($thisCsq{"Consequence"} =~ /feature_elongation/)) {
	    # duplication of one or more exons, upgrade to MODHIGH
	    $thisCsq{"IMPACT"} = "MODHIGH";
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

	    # some variants don't have ada_score or rf_score -> use a frac
	    my $minPassedFracSplicing = 0.5;
	    my $passed = 0;
	    my $totalPreds = 0;
	    # SpliceAI
	    my @DS_headers = ("SpliceAI_pred_DS_AG","SpliceAI_pred_DS_AL",
			      "SpliceAI_pred_DS_DG","SpliceAI_pred_DS_DL");
	    foreach my $ds (@DS_headers) {
		($thisCsq{$ds}) && (++$totalPreds) && last;
	    }
	    foreach my $ds (@DS_headers) {
		($thisCsq{$ds}) && ($thisCsq{$ds} > $spliceAI_cutoff) && (++$passed) && last;
	    }
	    # CADD-Splice
	    if ($thisCsq{"CADD_PHRED"}) {
		$totalPreds++;
		($thisCsq{"CADD_PHRED"} > $cadd_cutoff) && ($passed++);
	    }
	    # dbscSNV
	    if ($thisCsq{"ada_score"}) {
		$totalPreds++;
		($thisCsq{"ada_score"} > 0.6) && ($passed++);
	    }
	    if ($thisCsq{"rf_score"}) {
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

	    # upgrade severe TFBS and regulatory region variants to LOW or MODHIGH:
	    # TFBS_ablation and regulatory_region_ablation => MODHIGH 
	    # TF_binding_site_variant with [HIGH_INF_POS==Y] => MODHIGH
	    # TF_binding_site_variant without HIGH_INF_POS => LOW
	    # TFBS_amplification => LOW
	    # (regulatory_region_variant and regulatory_region_amplification stay MODIFIER)
	    elsif ($thisCsq{"Consequence"} =~ /TFBS_ablation|regulatory_region_ablation/) {
		$thisCsq{"IMPACT"} = "MODHIGH";
	    }
	    elsif (($thisCsq{"Consequence"} =~ /TF_binding_site_variant/) && 
		   ($thisCsq{"HIGH_INF_POS"}) && ($thisCsq{"HIGH_INF_POS"} eq "Y")) {
		$thisCsq{"IMPACT"} = "MODHIGH";
	    }
	    elsif ($thisCsq{"Consequence"} =~ /TF_binding_site_variant/) {
		# no HIGH_INF_POS
		$thisCsq{"IMPACT"} = "LOW";
	    }
	    elsif ($thisCsq{"Consequence"} =~ /TFBS_amplification/) {
		$thisCsq{"IMPACT"} = "LOW";
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
