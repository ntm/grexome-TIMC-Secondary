#!/usr/bin/perl

# 29/08/2019
# NTM

# Takes 4 args: a $candidatesFile xlsx, a gzipped tsv $transciptsFile
# as produced in Coverage_Data/, a merged bgzipped $gvcf (must be 
# tabix-indexed), and an $outDir that doesn't exist.
# For each sample (grexome*) found in $gvcf, produce in $outDir a 
# coverage*.tsv file with:
# for each coding transcript in $transcriptsFile that is transcribed
# from a candidate gene listed in $candidatesFile, print
# one line per exon (limited to CDS parts of exons),
# for each exon we report the percentage of its bases 
# (going to +-10 into neighboring introns) that are covered 
# at least 10x, 20x or 50x.
# A global line for the whole CDS is also printed (Exon==ALL).
# In addition, one every $sampleEveryN coding transcript 
# corresponding to non-candidate genes gets printed as a single
# line (with Exon==ALL).
# Any gene present in $candidatesFile gets 1 in column 
# KNOWN_CANDIDATE_GENE (others get 0).
# Non-coding transcripts are currently skipped.
# 2 final lines are printed, with the global stats for ALL_CANDIDATES
# and ALL_SAMPLED_GENES.

use strict;
use warnings;
use Spreadsheet::XLSX;


#########################################################

# fraction of non-candidate coding transcripts that are analyzed:
# we analyze one coding transcript every $sampleEveryN
my $sampleEveryN = 20;

#########################################################

(@ARGV == 4) || die "needs 4 args: a candidatesFile, a tsv.gz, a GVCF and an outDir\n";
my ($candidatesFile, $transcriptsFile, $gvcf, $outDir) = @ARGV;

(-f $candidatesFile) ||
    die "E: the supplied candidates file $candidatesFile doesn't exist\n";
(-f $transcriptsFile) ||
    die "transcriptsFile $transcriptsFile doesn't exist or isn't a file\n";
(-f $gvcf) ||
    die "GVCF $gvcf doesn't exist or isn't a file\n";
(-f "$gvcf.tbi") ||
    die "GVCF $gvcf exists but can't find its index file $gvcf.tbi, did you tabix-index the gvcf?\n";
(-e $outDir) && 
    die "outDir $outDir already exists, remove it or choose another name.\n";
mkdir($outDir) || die "cannot mkdir outDir $outDir\n";


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
# grab sample ids from GVCF header
my @samples;

open(GVCF,"gunzip -c $gvcf |") || 
    die "cannot gunzip-open GVCF $gvcf for reading\n";
while(my $line = <GVCF>) {
    ($line =~ /^##/) && next;
    ($line =~ /^#CHROM/) || die "problem with GVCF $gvcf header, not CHROM?\n$line\n";
    chomp($line);
    my @fields = split(/\t/,$line);
    @samples = @fields[9..$#fields];
    last;
}
close(GVCF);

# array of filehandles open for writing, one for each grexome,
# same indexes as @samples
my @outFHs;
foreach my $grexome (@samples) {
    my $outFile = "$outDir/coverage_$grexome.csv";
    open(my $FH, "> $outFile") || die "cannot open $outFile for writing\n";
    push(@outFHs, $FH);
    # print header
    print $FH "Gene\tKNOWN_CANDIDATE_GENE\tTranscript\tExon\tBases examined (+-10 around each exon)\tPercentage covered >= 50x\tPercentage covered >= 20x\tPercentage covered >= 10x\n";
}

#########################################################

# counters for global stats:
# for each grexome, we will count the number of bases of
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


# remember how many non-candidate coding transcripts were seen, for sampling
my $transcriptsSeen = 0;

open(GENES, "gunzip -c $transcriptsFile |") || 
    die "cannot gunzip-open transcriptsFile $transcriptsFile\n";


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
	# not a candidate gene: only examine one every $sampleEveryN
	if ($transcriptsSeen+1 == $sampleEveryN) {
	    $transcriptsSeen = 0;
	}
	else {
	    $transcriptsSeen++;
	    next;
	}
    }

    my @starts = split(/,/,$starts);
    my $numExons = @starts;
    my @ends = split(/,/,$ends);
    ($numExons == @ends) || 
	die "mismatch between numbers of starts and ends in line:\n$line\n";

    # we will print one line per exon for candidate genes, but also a final 
    # line for the whole gene/transcript
    my $lengthGene = 0;
    # for each grexome, count the bases whose cov is:
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

	# print single quote before gene name so excel doesn't corrupt file
	my $toPrint = "\'$gene\t";
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

	# count the bases of this exon covered at 50x / 20x / 10x / 0x for each grexome
	my @bases50Exon = (0) x scalar(@samples);
	my @bases20Exon = (0) x scalar(@samples);
	my @bases10Exon = (0) x scalar(@samples);
	my @bases0Exon = (0) x scalar(@samples);

	# range of interest: start at -10 and end at +10
	my $range = "$chr:".($starts[$i]-10)."-".($ends[$i]+10);
	# grab GVCF lines in our range of interest
	my $tabix = "tabix $gvcf $range";

	open(GVCF, "$tabix |") ||
	    die "cannot tabix-open the gvcf of interest with:\n$tabix\n";
	while (my $gvcfLine = <GVCF>) {
	    chomp($gvcfLine);
	    my @gvcfFields = split(/\t/,$gvcfLine);
	    my ($pos,$info,$format,@sampleData) = @gvcfFields[1,7..$#gvcfFields];
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

	    # %format: key is a FORMAT key (eg MIN_DP), value is the index of that key in $format
	    my %format;
	    { 
		my @format = split(/:/, $format);
		foreach my $i (0..$#format) {
		    $format{$format[$i]} = $i ;
		}
	    }

	    foreach my $i (0..$#samples) {
		# what's the coverage for this sample?
		my $coverage = 0;

		my @data = split(/:/,$sampleData[$i]);
		# if MIN_DP exists, use that
		if ((defined $format{"MIN_DP"}) && ($data[$format{"MIN_DP"}]) && ($data[$format{"MIN_DP"}] ne '.')) {
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
	close(GVCF);

	foreach my $i (0..$#samples) {
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
    my $toPrint = "\'$gene\t";
    if ($candidateGenes{$gene}) { $toPrint .= "1\t"; }
    else { $toPrint .= "0\t"; }
    $toPrint .= "$transcript\tALL\t$lengthGene\t";

    foreach my $i (0..$#samples) {
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
	    $lengthCandidates += $lengthGene;
	}
	else {
	    $bases50Sampled[$i] += $bases50Gene[$i];
	    $bases20Sampled[$i] += $bases20Gene[$i];
	    $bases10Sampled[$i] += $bases10Gene[$i];
	    $bases0Sampled[$i] += $bases0Gene[$i];
	    $lengthSampled += $lengthGene;
	}
    }
}

# print global stats
my $toPrint = "ALL_CANDIDATES\t1\tALL_CANDIDATES\tALL\t$lengthCandidates\t";

foreach my $i (0..$#samples) {
    my $frac = $bases50Candidates[$i] / $lengthCandidates;
    my $toPrintEnd = sprintf("%.2f",$frac)."\t";
    $frac = ($bases50Candidates[$i] + $bases20Candidates[$i]) / $lengthCandidates;
    $toPrintEnd .= sprintf("%.2f",$frac)."\t";
    $frac = ($bases50Candidates[$i] + $bases20Candidates[$i] + $bases10Candidates[$i]) / $lengthCandidates;
    $toPrintEnd .= sprintf("%.2f",$frac)."\n";
    print {$outFHs[$i]} $toPrint.$toPrintEnd;
}
# same for all sampled transcripts
$toPrint = "ALL_SAMPLED\t0\tALL_SAMPLED\tALL\t$lengthSampled\t";
foreach my $i (0..$#samples) {
    my $frac = $bases50Sampled[$i] / $lengthSampled;
    my $toPrintEnd = sprintf("%.2f",$frac)."\t";
    $frac = ($bases50Sampled[$i] + $bases20Sampled[$i]) / $lengthSampled;
    $toPrintEnd .= sprintf("%.2f",$frac)."\t";
    $frac = ($bases50Sampled[$i] + $bases20Sampled[$i] + $bases10Sampled[$i]) / $lengthSampled;
    $toPrintEnd .= sprintf("%.2f",$frac)."\n";
    print {$outFHs[$i]} $toPrint.$toPrintEnd;
}

foreach my $fh (@outFHs) {
    close($fh);
}

foreach my $gene (sort keys %candidateGenes) {
    if ($candidateGenes{$gene} == 1) {
	warn "W: all done but candidate gene $gene was never seen!\n";
    }
}

close(GENES);

