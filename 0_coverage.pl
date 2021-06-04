#!/usr/bin/perl

# 29/08/2019
# NTM

# Takes 5 args: a samples xlsx, a comma-separated list of $candidatesFiles xlsx,
# a gzipped tsv $transciptsFile as produced in Transcripts_Data/,
# a merged bgzipped $gvcf (must be  tabix-indexed), and an $outDir.
# If $outDir doesn't exist it is created; if it exists any pre-existing
# coverage*tsv file in it is not remade (considered good).
# For each sample found in $gvcf and lacking a coverage*tsv
# in $outDir, produce $outDir/coverage*.tsv with:
# for each coding transcript in $transcriptsFile that is transcribed
# from a gene listed in a $candidatesFiles or causal in $samplesFile,
# print one line per exon (limited to CDS parts of exons),
# for each exon we report the percentage of its bases 
# (going to +-10 into neighboring introns) that are covered 
# at least 10x, 20x or 50x.
# A global line for the whole CDS is also printed (Exon==ALL).
# In addition, every coding transcript corresponding to non-candidate
# genes gets printed as a single line (with Exon==ALL).
# Any gene present in $candidatesFile gets 1 in column 
# KNOWN_CANDIDATE_GENE (others get 0).
# Non-coding transcripts are currently skipped.
# 2 final lines are printed, with the global stats for ALL_CANDIDATES
# and ALL_GENES.

use strict;
use warnings;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use POSIX qw(strftime);

use lib "$RealBin";
use grexome_metaParse qw(parseCandidateGenes);

# use tabix module (installed in /usr/local/lib64/), this requires bioperl
use lib "/home/nthierry/Software/VariantEffectPredictor/ensembl-vep/";
use Bio::DB::HTS::Tabix;


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#########################################################

(@ARGV == 5) ||
    die "E $0: needs 5 args: a samples file, a comma-separated list of candidatesFiles, a tsv.gz, a GVCF and an outDir\n";
my ($samplesFile, $candidatesFiles, $transcriptsFile, $gvcf, $outDir) = @ARGV;

(-f $samplesFile) || die "E $0: the supplied samples file doesn't exist\n";
(-f $transcriptsFile) ||
    die "E $0: transcriptsFile $transcriptsFile doesn't exist or isn't a file\n";
(-f $gvcf) || die "E $0: GVCF $gvcf doesn't exist or isn't a file\n";
(-f "$gvcf.tbi") ||
    die "E $0: GVCF $gvcf exists but can't find its index file $gvcf.tbi, did you tabix-index the gvcf?\n";
(-d $outDir) || mkdir($outDir) || die "E $0: cannot mkdir outDir $outDir\n";

my $now = strftime("%F %T", localtime);
warn "I $0: $now - starting to run\n";

#########################################################
# parse known candidate genes file

# %candidateGenes: key==gene name, value==1 if gene is of interest
# but not seen yet, 2 if it's been seen
my %candidateGenes = ();


# $knownCandidateGenesR: hashref, key==$cohort, value is a hashref whose keys 
# are gene names and values are the "Confidence score" from a $candidatesFile,
# or 5 if the gene is "Causal" for a $cohort patient in $samplesFile.
my $knownCandidateGenesR = &parseCandidateGenes($candidatesFiles, $samplesFile);
foreach my $patho (keys(%$knownCandidateGenesR)) {
    foreach my $gene (keys($knownCandidateGenesR->{$patho})) {
	$candidateGenes{$gene} = 1;
    }
}

#########################################################
# grab sample ids from GVCF header
my $tabix = Bio::DB::HTS::Tabix->new( filename => "$gvcf" );

# ->header gives all headers in a single scalar, we just want the #CHROM line
my $header = $tabix->header;
($header =~ /\n(#CHROM\t.+)$/) ||
    die "E $0: cannot extract #CHROM line from header:\n$header\n";
my @fields = split(/\t/,$1);
my @samples = @fields[9..$#fields];

# samples that already have a coverage file are ignored (value 1)
my @samplesIgnored = (0) x scalar(@samples);

# array of filehandles open for writing, one for each non-ignored sample,
# same indexes as @samples
my @outFHs;
foreach my $i (0..$#samples) {
    my $sample = $samples[$i];
    my $outFile = "$outDir/coverage_$sample.csv";
    if (-e $outFile) {
	$samplesIgnored[$i] = 1;
    }
    else {
	open(my $FH, "> $outFile") || die "E $0: cannot open $outFile for writing\n";
	$outFHs[$i] = $FH;
	# print header
	print $FH "Gene\tKNOWN_CANDIDATE_GENE\tTranscript\tExon\tBases examined (+-10 around each exon)\tPercentage covered >= 50x\tPercentage covered >= 20x\tPercentage covered >= 10x\n";
    }
}

#########################################################

# counters for global stats:
# for each non-ignored sample, we will count the number of bases of
# candidate genes whose  coverage $cov is 
# >= 50x, 20x <= $cov < 50x, 10x <= $cov < 20x, or < 10x respectively
my @bases50Candidates = (0) x scalar(@samples);
my @bases20Candidates = (0) x scalar(@samples);
my @bases10Candidates = (0) x scalar(@samples);
my @bases0Candidates = (0) x scalar(@samples);
# same for all other sampled genes
my @bases50Sampled = (0) x scalar(@samples);
my @bases20Sampled = (0) x scalar(@samples);
my @bases10Sampled = (0) x scalar(@samples);
my @bases0Sampled = (0) x scalar(@samples);
# total lengths of candidates and of sampled transcripts
my ($lengthCandidates,$lengthSampled) = (0,0);


open(GENES, "gunzip -c $transcriptsFile |") || 
    die "E $0: cannot gunzip-open transcriptsFile $transcriptsFile\n";

# skip header
<GENES>;

while (my $line = <GENES>) {
    chomp($line);
    my @fields = split(/\t/,$line);
    (@fields == 8) || 
	die "E $0: wrong number of fields in line:\n$line\n";
    my ($transcript,$gene,$chr,$strand,$cdsStart,$cdsEnd,$starts,$ends) = @fields;

    # skip non-coding transcripts
    if ($cdsStart == $cdsEnd) {
	next;
    }

    # if gene is a candidate mark it as seen, or warn if it was seen earlier
    if (($candidateGenes{$gene}) && ($candidateGenes{$gene} == 1)) {
	$candidateGenes{$gene} = 2;
    }
    elsif ($candidateGenes{$gene}) {
	warn "W $0: found several transcripts for candidate gene $gene, is this expected?\n";
    }
    # else not a candidate gene, nothing to do

    my @starts = split(/,/,$starts);
    my $numExons = @starts;
    my @ends = split(/,/,$ends);
    ($numExons == @ends) || 
	die "E $0: mismatch between numbers of starts and ends in line:\n$line\n";

    # we will print one line per exon for candidate genes, but also a final 
    # line for the whole gene/transcript
    my $lengthGene = 0;
    # for each sample, count the bases whose cov is:
    # $cov >= 50x, 20x <= $cov < 50x, 10x <= $cov < 20x, or $cov < 10x respectively
    my @bases50Gene = (0) x scalar(@samples);
    my @bases20Gene = (0) x scalar(@samples);
    my @bases10Gene = (0) x scalar(@samples);
    my @bases0Gene = (0) x scalar(@samples);

    foreach my $i (0..$#starts) {
	# ignore 5'-UTR
	($ends[$i] < $cdsStart) && next;
	($starts[$i] < $cdsStart) && ($starts[$i] = $cdsStart);
	# ignore 3'-UTR
	($starts[$i] > $cdsEnd) && next;
	($ends[$i] > $cdsEnd) && ($ends[$i] = $cdsEnd);

	# apostrophe-space before gene name so excel doesn't corrupt file
	my $toPrint = "\' $gene\t";
	if ($candidateGenes{$gene}) {
	    $toPrint .= "1\t";
	}
	else {
	    $toPrint .= "0\t";
	}
	$toPrint .= "$transcript\t";
	# for Exon use eg 3\12 (so excel doesn't corrupt my file)
	$toPrint .= $i+1;
	$toPrint .= "\\$numExons\t";
	# end-start+1 is the exon length, add 10 bases on each side
	my $lengthExon = $ends[$i] - $starts[$i] + 21;
	$toPrint .= "$lengthExon\t";
	$lengthGene += $lengthExon;

	# count the bases of this exon covered at 50x / 20x / 10x / 0x for each sample
	my @bases50Exon = (0) x scalar(@samples);
	my @bases20Exon = (0) x scalar(@samples);
	my @bases10Exon = (0) x scalar(@samples);
	my @bases0Exon = (0) x scalar(@samples);

	# range of interest: start at -10 and end at +10
	my $range = "$chr:".($starts[$i]-10)."-".($ends[$i]+10);
	# grab GVCF lines in our range of interest
	my $tabixIter = $tabix->query($range);
	while (my $gvcfLine = $tabixIter->next) {
	    chomp($gvcfLine);
	    my @gvcfFields = split(/\t/,$gvcfLine);
	    my ($pos,$info,$format,@sampleData) = @gvcfFields[1,7..$#gvcfFields];
	    my $end = $pos;
	    # if this is a non-var block...
	    if ($info =~ /END=(\d+)/) {
		$end = $1;
	    }
	    # this line must overlap $range but it may go beyond, adjust if needed
	    if ($pos < $starts[$i]-10) {
		# tabix can return a line with a deletion that actually precedes 
		# $range, skip these lines
		($end >= $starts[$i]-10) || next;
		$pos = $starts[$i]-10;
	    }
	    if ($end > $ends[$i]+10) {
		#sanity
		($pos <= $ends[$i]+10) ||
		    die "E $0: tabix gave a line for $range that doesn't overlap it (after):\n$gvcfLine\n";
		$end = $ends[$i]+10;
	    }
	    # number of bases overlapping range in this line
	    my $numBases = $end - $pos + 1;

	    # %format: key is a FORMAT key (eg MIN_DP), value is the index of that key in $format
	    my %format;
	    { 
		my @format = split(/:/, $format);
		foreach my $i (0..$#format) {
		    $format{$format[$i]} = $i ;
		}
	    }

	    foreach my $i (0..$#samples) {
		# skip ignored samples
		($samplesIgnored[$i]) && next;
		# what's the coverage for this sample?
		my $coverage = 0;

		my @data = split(/:/,$sampleData[$i]);
		# if MIN_DP exists, use that
		if ((defined $format{"MIN_DP"}) && ($data[$format{"MIN_DP"}]) && ($data[$format{"MIN_DP"}] ne '.')) {
		    $coverage = $data[$format{"MIN_DP"}];
		}
		else {
		    # grab the depth (DP or DPI or sumOfADs, whichever is defined and higher)
		    if ((defined $format{"DP"}) && ($data[$format{"DP"}]) && ($data[$format{"DP"}] ne '.')) {
			$coverage = $data[$format{"DP"}];
		    }
		    if ((defined $format{"DPI"}) && ($data[$format{"DPI"}]) && ($data[$format{"DPI"}] ne '.')) {
			($coverage < $data[$format{"DPI"}]) && ($coverage = $data[$format{"DPI"}]);
		    }
		    if ((defined $format{"AD"}) && ($data[$format{"AD"}]) && ($data[$format{"AD"}] =~ /^[\d,]+$/)) {
			my $sumOfADs = 0;
			foreach my $ad (split(/,/,$data[$format{"AD"}])) {
			    $sumOfADs += $ad;
			}
			($coverage < $sumOfADs) && ($coverage = $sumOfADs);
		    }
		}
		if ($coverage >= 50) {
		    $bases50Exon[$i] += $numBases;
		}
		elsif ($coverage >= 20) {
		    $bases20Exon[$i] += $numBases;
		}
		elsif ($coverage >= 10) {
		    $bases10Exon[$i] += $numBases;
		}
		else {
		    $bases0Exon[$i] += $numBases;
		}
	    }		
	}

	foreach my $i (0..$#samples) {
	    ($samplesIgnored[$i]) && next;
	    # only print per-exon stats for candidate genes
	    if ($candidateGenes{$gene}) {
		my $frac = $bases50Exon[$i] / $lengthExon;
		my $toPrintEnd = sprintf("%.2f",$frac)."\t";
		$frac = ($bases50Exon[$i] + $bases20Exon[$i]) / $lengthExon;
		$toPrintEnd .= sprintf("%.2f",$frac)."\t";
		$frac = ($bases50Exon[$i] + $bases20Exon[$i] + $bases10Exon[$i]) / $lengthExon;
		$toPrintEnd .= sprintf("%.2f",$frac)."\n";
		print {$outFHs[$i]} $toPrint.$toPrintEnd;
	    }
	    # done with this exon but also record stats for the gene
	    $bases50Gene[$i] += $bases50Exon[$i];
	    $bases20Gene[$i] += $bases20Exon[$i];
	    $bases10Gene[$i] += $bases10Exon[$i];
	    $bases0Gene[$i] += $bases0Exon[$i];
	}
    }

    # done printing data for each exon of this gene (if it's a candidate),
    # now print global stats for $gene
    my $toPrint = "\' $gene\t";
    if ($candidateGenes{$gene}) { $toPrint .= "1\t"; }
    else { $toPrint .= "0\t"; }
    $toPrint .= "$transcript\tALL\t$lengthGene\t";

    foreach my $i (0..$#samples) {
	($samplesIgnored[$i]) && next;
	my $frac = $bases50Gene[$i] / $lengthGene;
	my $toPrintEnd = sprintf("%.2f",$frac)."\t";
	$frac = ($bases50Gene[$i] + $bases20Gene[$i]) / $lengthGene;
	$toPrintEnd .= sprintf("%.2f",$frac)."\t";
	$frac = ($bases50Gene[$i] + $bases20Gene[$i] + $bases10Gene[$i]) / $lengthGene;
	$toPrintEnd .= sprintf("%.2f",$frac)."\n";
	print {$outFHs[$i]} $toPrint.$toPrintEnd;

	# finally, update global stats
	if ($candidateGenes{$gene}) {
	    $bases50Candidates[$i] += $bases50Gene[$i];
	    $bases20Candidates[$i] += $bases20Gene[$i];
	    $bases10Candidates[$i] += $bases10Gene[$i];
	    $bases0Candidates[$i] += $bases0Gene[$i];
	}
	else {
	    $bases50Sampled[$i] += $bases50Gene[$i];
	    $bases20Sampled[$i] += $bases20Gene[$i];
	    $bases10Sampled[$i] += $bases10Gene[$i];
	    $bases0Sampled[$i] += $bases0Gene[$i];
	}
    }
    if ($candidateGenes{$gene}) {
	$lengthCandidates += $lengthGene;
    }
    else {
	$lengthSampled += $lengthGene;
    }
}

# print global stats
my $toPrint = "ALL_CANDIDATES\t1\tALL_CANDIDATES\tALL\t$lengthCandidates\t";

foreach my $i (0..$#samples) {
    ($samplesIgnored[$i]) && next;
    my $frac = $bases50Candidates[$i] / $lengthCandidates;
    my $toPrintEnd = sprintf("%.2f",$frac)."\t";
    $frac = ($bases50Candidates[$i] + $bases20Candidates[$i]) / $lengthCandidates;
    $toPrintEnd .= sprintf("%.2f",$frac)."\t";
    $frac = ($bases50Candidates[$i] + $bases20Candidates[$i] + $bases10Candidates[$i]) / $lengthCandidates;
    $toPrintEnd .= sprintf("%.2f",$frac)."\n";
    print {$outFHs[$i]} $toPrint.$toPrintEnd;
}
# same for all sampled transcripts
$toPrint = "ALL_GENES\t0\tALL_GENES\tALL\t$lengthSampled\t";
foreach my $i (0..$#samples) {
    ($samplesIgnored[$i]) && next;
    my $frac = $bases50Sampled[$i] / $lengthSampled;
    my $toPrintEnd = sprintf("%.2f",$frac)."\t";
    $frac = ($bases50Sampled[$i] + $bases20Sampled[$i]) / $lengthSampled;
    $toPrintEnd .= sprintf("%.2f",$frac)."\t";
    $frac = ($bases50Sampled[$i] + $bases20Sampled[$i] + $bases10Sampled[$i]) / $lengthSampled;
    $toPrintEnd .= sprintf("%.2f",$frac)."\n";
    print {$outFHs[$i]} $toPrint.$toPrintEnd;
}

foreach my $fh (@outFHs) {
    ($fh) && close($fh);
}

foreach my $gene (sort keys %candidateGenes) {
    if ($candidateGenes{$gene} == 1) {
	warn "W $0: all done but candidate gene $gene was never seen!\n";
    }
}

close(GENES);
$tabix->close;

$now = strftime("%F %T", localtime);
warn "I $0: $now - ALL DONE, completed successfully!\n";
