#!/usr/bin/perl

# 29/08/2019
# NTM

# Takes 3 args: a $candidatesFile xlsx, a gzipped tsv $transciptsFile as
# produced in Coverage_Data/, and a $gvcf (must be tabix-indexed).
# Print to stdout a tsv: 
# for each coding transcript in $transcriptsFile that is transcribed
# from a candidate gene listed in $candidatesFile, print
# one line per exon (limited to CDS parts of exons),
# for each exon we report the percentage of its bases 
# (going to +-10 into neighboring introns) that are covered 
# at least 10x or 20x.
# A global line for the whole CDS is also printed (Exon==ALL).
# In addition, any coding transcript corresponding to a non-candidate 
# gene gets printed as a single line (with Exon==ALL) with probability 
# $probaSample.
# Any gene present in $candidatesFile gets 1 in column 
# KNOWN_CANDIDATE_GENE (others get 0).
# Non-coding transcripts are currently skipped.
# 2 final lines are printed, with the global stats for ALL_CANDIDATES
# and ALL_SAMPLED_GENES.

use strict;
use warnings;
use Spreadsheet::XLSX;


#########################################################

# probability that a non-candidate coding transcript is analyzed
my $probaSample = 0.05;

#########################################################

(@ARGV == 3) || die "needs 3 args: a candidatesFile, a tsv.gz and a GVCF\n";
my ($candidatesFile, $transcriptsFile, $gvcf) = @ARGV;

(-f $candidatesFile) ||
    die "E: the supplied candidates file $candidatesFile doesn't exist\n";
(-f $transcriptsFile) ||
    die "transcriptsFile $transcriptsFile doesn't exist or isn't a file\n";
(-f $gvcf) ||
    die "GVCF $gvcf doesn't exist or isn't a file\n";
(-f "$gvcf.tbi") ||
    die "GVCF $gvcf exists but can't find its index file $gvcf.tbi, did you tabix-index the gvcf?\n";


#########################################################
# parse known candidate genes file

# %candidateGenes: key==gene name, value==1 if gene is of interest
# but not seen yet, 2 if it's been seen
my %candidateGenes = ();

{
    # code adapted from 6_extractCohorts.pl
    my $workbook = Spreadsheet::XLSX->new("$candidatesFile");
    (defined $workbook) ||
	die "E when parsing xlsx\n";
    ($workbook->worksheet_count() == 1) ||
	die "E parsing xlsx: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($geneCol) = (-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	($cell->value() eq "Candidate gene") &&
	    ($geneCol = $col);
     }
    ($geneCol >= 0) ||
	die "E parsing xlsx: no col title is Candidate gene\n";
    
    foreach my $row ($rowMin+1..$rowMax) {
	my $gene = $worksheet->get_cell($row, $geneCol)->unformatted();
	# clean up $gene
	$gene =~ s/^\s+//;
	$gene =~ s/\s+$//;

	$candidateGenes{$gene} = 1;
    }
}


#########################################################

# counters for global stats:
# we will count the number of bases of candidate genes whose 
# coverage $cov is >= 20x, 10x <= $cov < 20x, or < 10x respectively
my ($bases20Candidates,$bases10Candidates,$bases0Candidates) = (0,0,0);
# same for all other sampled genes
my ($bases20Sampled,$bases10Sampled,$bases0Sampled) = (0,0,0);


open(GENES, "gunzip -c $transcriptsFile |") || 
    die "cannot gunzip-open transcriptsFile $transcriptsFile\n";

# print header
print "Gene\tKNOWN_CANDIDATE_GENE\tTranscript\tExon\tBases examined (+-10 around each exon)\tPercentage covered >= 20x\tPercentage covered >= 10x\n";

while (my $line = <GENES>) {
    chomp($line);
    my @fields = split(/\t/,$line);
    (@fields == 7) || 
	die "wrong number of fields in line:\n$line\n";
    my ($transcript,$chr,$cdsStart,$cdsEnd,$starts,$ends,$gene) = @fields;

    # skip non-coding transcripts
    if ($cdsStart == $cdsEnd) {
	next;
    }

    # if gene is a candidate mark it as seen, or warn if it was seen earlier
    if (($candidateGenes{$gene}) && ($candidateGenes{$gene} == 1)) {
	$candidateGenes{$gene} = 2;
    }
    elsif ($candidateGenes{$gene}) {
	warn "W: found several transcripts for candidate gene $gene, is this expected?\n";
    }
    else {
	# not a candidate gene: only examine it with probability $probaSample
	(rand(1) > $probaSample) && next;
    }

    my @starts = split(/,/,$starts);
    my $numExons = @starts;
    my @ends = split(/,/,$ends);
    ($numExons == @ends) || 
	die "mismatch between numbers of starts and ends in line:\n$line\n";

    # we will print one line per exon for candidate genes, but also a final 
    # line for the whole gene/transcript
    my $lengthGene = 0;
    # we will count the number of bases in this gene whose coverage $cov is:
    # $cov >= 20x, 10x <= $cov < 20x, or $cov < 10x respectively
    my ($bases20Gene,$bases10Gene,$bases0Gene) = (0,0,0);

    foreach my $i (0..$#starts) {
	# ignore 5'-UTR
	($ends[$i] < $cdsStart) && next;
	($starts[$i] < $cdsStart) && ($starts[$i] = $cdsStart);
	# ignore 3'-UTR
	($starts[$i] > $cdsEnd) && next;
	($ends[$i] > $cdsEnd) && ($ends[$i] = $cdsEnd);

	my $toPrint = "$gene\t";
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
	my $length = $ends[$i] - $starts[$i] + 21;
	$toPrint .= "$length\t";
	$lengthGene += $length;

	# we will count the number of bases in this exon whose coverage $cov is:
	# $cov >= 20x, 10x <= $cov < 20x, or $cov < 10x respectively
	my ($basesCovered20,$basesCovered10,$basesCovered0) = (0,0,0);

	# range of interest: start at -10 and end at +10
	my $range = "$chr:".($starts[$i]-10)."-".($ends[$i]+10);
	# grab GVCF lines in our range of interest
	my $tabix = "tabix $gvcf $range";

	open(GVCF, "$tabix |") ||
	    die "cannot tabix-open the gvcf of interest with:\n$tabix\n";
	while (my $gvcfLine = <GVCF>) {
	    chomp($gvcfLine);
	    my @gvcfFields = split(/\t/,$gvcfLine);
	    my ($pos,$info,$format,$data) = @gvcfFields[1,7..9];
	    my $end = $pos;
	    # if this is a non-var block...
	    if ($info =~ /END=(\d+);/) {
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
		    die "tabix gave a line for $range that doesn't overlap it (after):\n$gvcfLine\n";
		$end = $ends[$i]+10;
	    }
	    # number of bases overlapping range in this line
	    my $numBases = $end - $pos + 1;

	    # what's the coverage?
	    my $coverage = 0;
	    # %format: key is a FORMAT key (eg MIN_DP), value is the index of that key in $format
	    my %format;
	    { 
		my @format = split(/:/, $format);
		foreach my $i (0..$#format) {
		    $format{$format[$i]} = $i ;
		}
	    }
	    my @data = split(/:/,$data);
	    # if MIN_DP exists, use that
	    if (defined $format{"MIN_DP"}) {
		$coverage = $data[$format{"MIN_DP"}];
	    }
	    else {
		# grab the depth (DP or DPI, whichever is defined and higher)
		if ((defined $format{"DP"}) && ($data[$format{"DP"}]) && ($data[$format{"DP"}] ne '.')) {
		    $coverage = $data[$format{"DP"}];
		}
		if ((defined $format{"DPI"}) && ($data[$format{"DPI"}]) && ($data[$format{"DPI"}] ne '.') &&
		    ($data[$format{"DPI"}] > $coverage)) {
		    $coverage = $data[$format{"DPI"}];
		}
	    }

	    if ($coverage >= 20) {
		$basesCovered20 += $numBases;
	    }
	    elsif ($coverage >= 10) {
		$basesCovered10 += $numBases;
	    }
	    else {
		$basesCovered0 += $numBases;
	    }

	}
	close(GVCF);

	# we didn't deal with the Strelka call-HR-before-indel bug, and we
	# didn't take into account the REF size (for DELs or MNPs), so the
	# result will be approximate... still, make sure we're close enough
	my $basesTotal = $basesCovered20 + $basesCovered10 + $basesCovered0;
	if ( (abs($length - $basesTotal) / $length) > 0.1) {
	    warn "W: more than 10% difference between bases counted ($basesTotal) and length for: $toPrint\n";
	}
	my $frac = $basesCovered20 / $basesTotal;
	$toPrint .= sprintf("%.2f",$frac)."\t";
	$frac = ($basesCovered20 + $basesCovered10) / $basesTotal;
	$toPrint .= sprintf("%.2f",$frac)."\n";

	# only print per-exon stats for candidate genes
	($candidateGenes{$gene}) && (print $toPrint);
	# done with this exon but also record stats for the gene
	$bases20Gene += $basesCovered20;
	$bases10Gene += $basesCovered10;
	$bases0Gene += $basesCovered0;
    }

    # done printing data for each exon of this gene (if it's a candidate),
    # now print global stats for $gene
    my $toPrint = "$gene\t";
    if ($candidateGenes{$gene}) { $toPrint .= "1\t"; }
    else { $toPrint .= "0\t"; }
    $toPrint .= "$transcript\tALL\t$lengthGene\t";
    my $basesTotalGene = $bases20Gene + $bases10Gene + $bases0Gene;
    if ( (abs($lengthGene - $basesTotalGene) / $lengthGene) > 0.1) {
	warn "W: GENE LEVEL, more than 10% difference between bases counted ($basesTotalGene) and length==$lengthGene for: $gene\n";
    }
    my $frac = $bases20Gene / $basesTotalGene;
    $toPrint .= sprintf("%.2f",$frac)."\t";
    $frac = ($bases20Gene + $bases10Gene) / $basesTotalGene;
    $toPrint .= sprintf("%.2f",$frac)."\n";
    print $toPrint;

    # finally, update global stats
    if ($candidateGenes{$gene}) {
	$bases20Candidates += $bases20Gene;
	$bases10Candidates += $bases10Gene;
	$bases0Candidates += $bases0Gene;
    }
    else {
	$bases20Sampled += $bases20Gene;
	$bases10Sampled += $bases10Gene;
	$bases0Sampled += $bases0Gene;
    }
}

# print global stats
my $toPrint = "ALL_CANDIDATES\t1\tALL_CANDIDATES\tALL\t";
my $basesTotal = $bases20Candidates + $bases10Candidates + $bases0Candidates;
$toPrint .= "$basesTotal\t";
my $frac = $bases20Candidates / $basesTotal;
$toPrint .= sprintf("%.2f",$frac)."\t";
$frac = ($bases20Candidates + $bases10Candidates) / $basesTotal;
$toPrint .= sprintf("%.2f",$frac)."\n";
print $toPrint;
# same for all sampled transcripts
$toPrint = "ALL_SAMPLED\t0\tALL_SAMPLED\tALL\t";
$basesTotal = $bases20Sampled + $bases10Sampled + $bases0Sampled;
$toPrint .= "$basesTotal\t";
$frac = $bases20Sampled / $basesTotal;
$toPrint .= sprintf("%.2f",$frac)."\t";
$frac = ($bases20Sampled + $bases10Sampled) / $basesTotal;
$toPrint .= sprintf("%.2f",$frac)."\n";
print $toPrint;


foreach my $gene (sort keys %candidateGenes) {
    if ($candidateGenes{$gene} == 1) {
	warn "W: all done but candidate gene $gene was never seen!\n";
    }
}

close(GENES);

