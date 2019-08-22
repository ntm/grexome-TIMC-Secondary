#!/usr/bin/perl

# 19/08/2019, but starting from 6_extractCohorts.pl
# NTM

# Takes 2 arguments: $inDir $outDir
# - $inDir must contain gzipped cohort TSVs as produced by extractCohorts.pl;
# - $outDir doesn't exist, it will be created and filled with one TSV
#   per infile (never gzipped), adding .Transcripts to the name.
#
# Each cohort tsv is filtered with 8_finalFilters.pl to consider
# only rare variants (max_af_*) in picked transcripts, that aren't seen
# in too many CTRLs (max_ctrl_*), that are well genotyped in our dataset (min_hr),
# and that have HIGH or MODERATE impact (no_mod no_low).
# We then produce one TSV for each cohort.
# In each TSV we print one line for each transcript (=="Feature"), with:
# - COUNT_$cohort_HV_HIGH = number of distinct samples with at
#   least one HV HIGH variant
# - COUNT_$cohort_HV_MODER = number of distinct samples with at
#   least one HV MODERATE-or-HIGH variant
# - COUNT_$cohort_HET_HIGH = number of distinct samples with at
#   least TWO HET (or one HV) HIGH variants
# - COUNT_$cohort_HET_MODER = number of distinct samples with at
#   least TWO HET (or one HV) MODERATE-or-HIGH variants
# - 4 more columns COUNT_NEGCTRL_* with similar counts but counting
#   the control samples (using the extractCohorts criteria).
# - HV_HIGH, HV_MODER, HET_HIGH, HET_MODER: list of samples counted
#   in the corresponding COUNT columns (so the *MODER columns contain
#   the corresponding *HIGH samples, and the HET* columns contain
#   the corresponding HV* samples, but it's ok), we don't list the
#   NEGCTRL samples.
# In addition, several useful columns from the source file are
# retained, see @keptColumns.
# A transcript line is not printed if all 4 COUNT_$cohort columns
# are zero.

use strict;
use warnings;
use POSIX qw(strftime);
use Parallel::ForkManager;


# number of jobs
my $numJobs = 8;

# full path to finalFilter.pl, unfortunately hard-coded for now
# but with several possibbilities
my $filterBin;
my @possibleFilterBins = ("/home/nthierry/PierreRay/Grexome/SecondaryAnalyses/8_finalFilters.pl",
			  "/home/nthierry/VariantCalling/GrexomeFauve/SecondaryAnalyses/8_finalFilters.pl");
foreach my $f (@possibleFilterBins) {
    (-f $f) && ($filterBin = "perl $f") && last;
}
($filterBin) || 
    die "Sorry, can't find 8_finalFilters.pl, update \@possibleFilterBins\n";

# columns we want to keep, in this order:
my @keptColumns = qw(SYMBOL KNOWN_CANDIDATE_GENE Feature Gene RefSeq BIOTYPE);
# in addition we insert the new COUNT* columns right after the last @keptColumns
# and immediately followed by the HV_HIGH et al colums, and we then copy all 
# the GTEX_* columns (in the same order as in infile)

# for convenience, build a hash:
# key == a keptColumn, value== the column index where we want it in the output
my %keptCols;
foreach my $i (0..$#keptColumns) {
    $keptCols{$keptColumns[$i]} = $i;
}

# also for convenience: the types of samples to count
my @countTypes = ("HV_HIGH","HV_MODER","HET_HIGH","HET_MODER");

#########################################################

(@ARGV == 2) || die "needs 2 args: an inDir and a non-existant outDir\n";
my ($inDir, $outDir) = @ARGV;
(-d $inDir) ||
    die "inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "cannot opendir inDir $inDir\n";
(-e $outDir) && 
    die "found argument outDir $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "cannot mkdir outDir $outDir\n";


# create fork manager
my $pm = new Parallel::ForkManager($numJobs);


while (my $inFile = readdir(INDIR)) {
    ($inFile =~ /^\./) && next;
    $pm->start && next;
    my $cohort;
    if ($inFile =~ (/^(\w+)\.csv\.gz$/)) {
	$cohort = $1;
    }
    else {
	warn "W: cannot parse filename of inFile $inFile, skipping it\n";
	$pm->finish;
    }

    my $outFile = "$cohort.Transcripts.csv" ;
    open(OUT, "> $outDir/$outFile") ||
	die "E cannot open outfile $outDir/$outFile: $!\n";


    my $now = strftime("%F %T", localtime);
    warn "I: $now - starting $0 on $cohort\n";

    my $com = "$filterBin --max_ctrl_hv 10 --max_ctrl_het 50 --min_hr 100 --no_mod --no_low --pick";
    # using defaults for AFs 
    # $com .= " --max_af_gnomad 0.01 --max_af_1kg 0.03 --max_af_esp 0.05"
    open(FILTER, "gunzip -c $inDir/$inFile | $com | ") ||
	die "E: cannot gunzip-open and filter infile $inFile\n";

    # header line
    my $header = <FILTER>;
    chomp($header);
    my @headers = split(/\t/,$header);

    # $transCol == column of "Feature" (== transcript)
    my $transCol;
    # column of IMPACT
    my $impactCol;
    # columns of HV,HET,NEGCTRL_HV,NEGCTRL_HET
    my ($colHv,$colHet,$colNegHv,$colNegHet);
    # @destCols: for each column $i in infile: 
    # if $destCols[$i] >= 0 it is the column where that info goes in outfile
    # if $destCols[$i] == -1 source column is a GTEX
    my @destCols = ();

    # new headers
    my $newHeaders = join("\t",@keptColumns);
    # add COUNTs and SAMPLES
    foreach my $ct (@countTypes) {
	$newHeaders .= "\tCOUNT_$cohort"."_$ct";
    }
    foreach my $ct (@countTypes) {
	$newHeaders .= "\tCOUNT_NEGCTRL_$ct";
    }
    foreach my $ct (@countTypes) {
	$newHeaders .= "\t$ct";
    }
    # GTEX headers are added when parsing $headers and filling @destCols
    foreach my $hi (0..$#headers) {
	if (defined $keptCols{$headers[$hi]}) {
	    $destCols[$hi] = $keptCols{$headers[$hi]};
	    ($headers[$hi] eq "Feature") && ($transCol = $hi);
	}
	elsif ($headers[$hi] =~ /^GTEX_/) {
	    $destCols[$hi] = -1;
	    $newHeaders .= "\t$headers[$hi]";
	}
	elsif ($headers[$hi] eq "IMPACT") {
	    $impactCol = $hi;
	}
	elsif ($headers[$hi] eq "HV") {
	    $colHv = $hi;
	}
	elsif ($headers[$hi] eq "HET") {
	    $colHet = $hi;
	}
	elsif ($headers[$hi] eq "NEGCTRL_HV") {
	    $colNegHv = $hi;
	}
	elsif ($headers[$hi] eq "NEGCTRL_HET") {
	    $colNegHet = $hi;
	}
	# else ignore this column
    }

    # add filter params at the end and print
    print OUT "$newHeaders\t$headers[$#headers]\n";

    # body
    # lines desribing each transcript can be interspersed...
    # so we need %transcript2*, we will print when we change chroms
    # key==$transcript, value==start of line to print
    my %transcript2start;
    # key==$transcript, value==GTEX data to print (starting with \t)
    my %transcript2gtex;
    # key==$transcript, value is an arrayref with 8 hashrefs,
    # one each for @countTypes and then again for NEGCTRL_*,
    # each hash has key==$grexome, value==number of variants (of that
    # category), MODER lists the samples that are MODERATE-or-HIGH,
    # HET also lists the samples that are HV (but these count as 2 variants)
    my %transcript2samples;

    # chrom in previous line, when new line is different we print and empty hashes
    my $prevChr = "chr1";

    while(1) {
	# grab lines inside loop so we can print stuff for last chrom at EOF
	my $line = <FILTER>;
	my @fields;
	my $chr;
	if (defined $line) {
	    chomp($line);
	    @fields = split(/\t/,$line);
	    # check chrom, POS must be first field (not checked but...)
	    ($fields[0] =~ /^(chr[^:]+:)\d+$/) ||
		die "E: cannot grab chrom in line:\n$line\n";
	    $chr = $1;
	}
	else {
	    # EOF, use bogus chr
	    $chr = "chrBOGUS";
	}

	if ($chr ne $prevChr) {
	    foreach my $transcript (sort keys %transcript2start) {
		my $toPrint = $transcript2start{$transcript};
		# we will print except if all 4 COUNT_$cohort cols are zero
		my $printOK = 0;

		# for HET counts and lists we only want samples with at least 2 variants
		foreach my $t2si (2,3,6,7) {
		    foreach my $sample (keys %{$transcript2samples{$transcript}->[$t2si]}) {
			($transcript2samples{$transcript}->[$t2si]->{$sample} >= 2) ||
			    (delete $transcript2samples{$transcript}->[$t2si]->{$sample});
		    }
		}
		# the COUNT_* values are now simple:
		foreach my $t2si (0..7) {
		    $toPrint .= "\t".scalar(keys %{$transcript2samples{$transcript}->[$t2si]});
		}
		# and the samples lists also
		foreach my $t2si (0..3) {
		    $toPrint .= "\t".join(',',sort(keys(%{$transcript2samples{$transcript}->[$t2si]})));
		    (scalar(keys(%{$transcript2samples{$transcript}->[$t2si]})) != 0) && ($printOK = 1);
		}
		
		# transcript2gtex already starts with \t
		$toPrint .= $transcript2gtex{$transcript};

		($printOK) && (print OUT "$toPrint\n");

		# OK, clear hash entries for this transcript (not %transcript2start,
		# it will be cleared after the loop)
		delete($transcript2gtex{$transcript});
		delete($transcript2samples{$transcript});
	    }

	    # sanity, all hashes should be empty
	    (keys %transcript2gtex) && 
		die "E: finished printing chr $prevChr but still have keys in t2gtex %transcript2gtex\n";
	    (keys %transcript2samples) && 
		die "E: finished printing chr $prevChr but still have keys in t2samples %transcript2samples\n";
	    # clear %transcript2start and set $prevChr
	    %transcript2start = ();
	    if ($chr eq "chrBOGUS") {
		# all done with this infile
		last;
	    }
	    # else keep going
	    $prevChr = $chr;
	    # NOT next, still want to process and store this line
	}

	# process line
	my $transcript = $fields[$transCol];
	if (! $transcript2start{$transcript}) {
	    # first time we see $transcript, construct start and gtex strings
	    my @start = ();
	    my $gtex = "";
	    foreach my $fi (0..$#fields) {
		if (!defined($destCols[$fi])) {
		    next;
		}
		elsif ($destCols[$fi] >= 0) {
		    $start[$destCols[$fi]] = $fields[$fi];
		}
		elsif ($destCols[$fi] == -1) {
		    $gtex .= "\t$fields[$fi]";
		}
	    }
	    $transcript2start{$transcript} = join("\t",@start);
	    $transcript2gtex{$transcript} = $gtex;
	    # also initialize samples: arrayref with 8 hashrefs to empty hashes
	    $transcript2samples{$transcript} = [{},{},{},{},{},{},{},{}];
	}

	# in any case, update %transcript2samples
	my $impact = $fields[$impactCol];

	foreach my $col ($colHv,$colHet,$colNegHv,$colNegHet) {
	    # all columns are processed very similarly, we need some flags
	    # $negctrl is a flag, 1 if NEGCTRL_*, 0 otherwise
	    my $isNegctrl = 0;
	    if (($col == $colNegHv) || ($col == $colNegHet)) {
		$isNegctrl = 1;
	    }
	    # $isHV == 1 if HV, 0 if HET
	    my $isHV = 0;
	    if (($col == $colHv) || ($col == $colNegHv)) {
		$isHV = 1;
	    }

	    my $samples = $fields[$col];
	    ($samples) || next;
	    # remove genotype, we only want the grexome IDs
	    $samples =~ s/^\d+\/\d+~//;
	    
	    foreach my $sample (split(/,/,$samples)) {
		# ignore [DP;AF]
		($sample =~ /^(grexome\d\d\d\d)/) || 
		    die "E: cannot grab grexome from sample $sample\n";
		my $grexome = $1;
		if ($impact eq "HIGH") {
		    # always initialize to zero if needed before incrementing
		    if ($isHV) {
			# COUNT_HV_HIGH and COUNT_HV_MODER get +1
			foreach my $ai (4*$isNegctrl, 4*$isNegctrl+1) {
			    ($transcript2samples{$transcript}->[$ai]->{$grexome}) || 
				($transcript2samples{$transcript}->[$ai]->{$grexome} = 0);
			    $transcript2samples{$transcript}->[$ai]->{$grexome}++;
			}
		    }
		    # whether $isHV or not, COUNT_HET_HIGH and COUNT_HET_MODER get updated,
		    # but HV variants count as 2
		    foreach my $ai (4*$isNegctrl+2, 4*$isNegctrl+3) {
			($transcript2samples{$transcript}->[$ai]->{$grexome}) || 
			    ($transcript2samples{$transcript}->[$ai]->{$grexome} = 0);
			$transcript2samples{$transcript}->[$ai]->{$grexome} += (1+$isHV);
		    }
		}
		elsif ($impact eq "MODERATE") {
		    if ($isHV) {
			# COUNT_HV_MODER gets +1
			my $ai = 4*$isNegctrl+1;
			($transcript2samples{$transcript}->[$ai]->{$grexome}) || 
			    ($transcript2samples{$transcript}->[$ai]->{$grexome} = 0);
			$transcript2samples{$transcript}->[$ai]->{$grexome}++;
		    }
		    # whether $isHV or not, COUNT_HET_MODER gets +1 or +2
		    my $ai = 4*$isNegctrl+3;
		    ($transcript2samples{$transcript}->[$ai]->{$grexome}) || 
			($transcript2samples{$transcript}->[$ai]->{$grexome} = 0);
		    $transcript2samples{$transcript}->[$ai]->{$grexome} += (1+$isHV);
		}
		# else: ignore other IMPACTs but they were filtered anyways
	    }
	}
    }

    close(OUT);
    $now = strftime("%F %T", localtime);
    warn "I: $now - Finished with $cohort\n";
    $pm->finish;
}
closedir(INDIR);

$pm->wait_all_children;

