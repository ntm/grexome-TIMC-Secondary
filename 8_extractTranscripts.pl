#!/usr/bin/perl

# 19/08/2019, but starting from 6_extractCohorts.pl
# NTM

# Takes 2 arguments: $inDir $outDir
# - $inDir must contain cohort TSVs as produced by extractCohorts.pl,
#   filtered and reordered by filterVariants.pl and reorderColumns.pl;
# - $outDir doesn't exist, it will be created and filled with one TSV
#   per infile, adding .Transcripts to the name.
#
# We expect that cohort files were filtered to consider only rare variants 
# (max_af_*) in picked transcripts, that aren't seen in too many CTRLs 
# (max_ctrl_*) and that are well genotyped in our dataset (min_hr).
#
# We then produce one TSV for each cohort.
# In each TSV we print one line for each transcript (=="Feature"), with:
# - COUNT_HV_HIGH = number of distinct samples with at least
#   one HV HIGH variant
# - COUNT_HV_MODHIGH = number of distinct samples with at
#   least one HV MODHIGH-or-HIGH variant
# - COUNT_HV_MODER = number of distinct samples with at
#   least one HV MODERATE-or-MODHIGH-or-HIGH variant
# - COUNT_COMPHET_HIGH = number of distinct samples with at
#   least TWO HET (or one HV) HIGH variants
# - COUNT_COMPHET_MODHIGH = number of distinct samples with at
#   least TWO HET (or one HV) MODHIGH-or-HIGH variants
# - COUNT_COMPHET_MODER = number of distinct samples with at
#   least TWO HET (or one HV) MODERATE-or-MODHIGH-or-HIGH variants
# - 6 more columns COUNT_OTHERCAUSE_* with similar counts but counting
#   the samples with a "known causal variant" in another gene
# - another 6 columns COUNT_COMPAT_* with similar counts but counting
#   the samples belonging to compatible cohorts (as defined in extractCohorts)
# - 6 final columns COUNT_NEGCTRL_* with similar counts but counting
#   the control samples (as in extractCohorts again).
# - HV_HIGH, HV_MODHIGH, HV_MODER, COMPHET_HIGH, COMPHET_MODHIGH, COMPHET_MODER:
#   non-redundant list of samples counted in the corresponding COUNT columns
#   (so eg the *MODHIGH columns don't list the corresponding *HIGH samples,
#   even though COUNT*MODHIGH counts them);
# - OTHERCAUSE_COMPHET_MODHIGH, COMPAT_COMPHET_MODHIGH, NEGCTRL_COMPHET_MODHIGH:
#   all samples counted in corresponding COUNT columns (we don't list the
#   MODER samples).
# In addition, several useful columns from the source file are
# retained, see @keptColumns.
# A transcript line is not printed if all 6 COUNT_$cohort columns
# are zero.

use strict;
use warnings;
use POSIX qw(strftime);


#############################################
## hard-coded stuff that shouldn't change much

# @compatible: array of arrayrefs, each arrayref holds cohorts that
# should NOT be used as neg controls for each other.
# The cohort names must match the "pathology" column of the $metadata xlsx
# (this is checked).
# NOTE: @compatible IS DUPLICATED IN 6_extractCohorts.pl, IF IT IS CHANGED HERE IT 
# MUST ALSO BE CHANGED THERE 
my @compatible = (["Flag","Astheno","Headless"],
		  ["Azoo","Ovo","Macro","IOP"],
		  ["Globo","Macro","Terato"]);


# columns we want to keep, in this order:
my @keptColumns = qw(SYMBOL KNOWN_CANDIDATE_GENE Feature Gene RefSeq BIOTYPE);
# in addition we insert the new COUNT* columns right after the last @keptColumns
# and immediately followed by the HV_HIGH et al colums, and we then copy all 
# the GTEX_* columns (in the same order as in infile)

# among the @keptColumns some have cohort-specific data: list them
my @keptColumnsSpecific = qw(KNOWN_CANDIDATE_GENE);


# also for convenience: the types of samples to count
my @countTypes = ("HV_HIGH","HV_MODHIGH","HV_MODER","COMPHET_HIGH","COMPHET_MODHIGH","COMPHET_MODER");

#########################################################
# pre-process some of the hard-coded stuff

# build hashes from @keptColumns*:
# %keptCols for common data, %keptColsSpecific for cohort-specific data,
# keys == a keptColumn* title, value== the column index where we want it in the output
my %keptCols;
foreach my $i (0..$#keptColumns) {
    $keptCols{$keptColumns[$i]} = $i;
}
my %keptColsSpecific;
foreach my $col (@keptColumnsSpecific) {
    ($keptCols{$col}) || 
	die "E: column $col is in keptColsSpecific but not in keptColumns, add it there\n";
    $keptColsSpecific{$col} = $keptCols{$col};
    delete($keptCols{$col});
}

# %compatible: key is a cohort name, value is a hashref
# with keys == cohorts that shouldn't be used as negative 
# controls for this cohort, value==1
my %compatible = ();

foreach my $notConR (@compatible) {
    foreach my $cohort (@$notConR) {
	(defined $compatible{$cohort}) || ($compatible{$cohort} = {});
	foreach my $notC (@$notConR) {
	    ($notC eq $cohort) && next;
	    $compatible{$cohort}->{$notC} = 1;
	}
    }
}

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


# Accumulators for all the data we want to print:
# lines describing each transcript in a given infile can be interspersed,
# and we need to parse all infiles before printing anything (for the 
# NEGCTRL and COMPAT counters, otherwise we can't count variants
# that don't occur in $cohort).
# Some data concerning a given transcript must be identical in every infile, 
# some other is cohort-specific...

# COMMON DATA: 
# %transcript2start: key==$transcript, value==arrayref holding the start-of-line 
# common data for this transcript (one value per column), the cohort-specific
# fields are left undefined (only KNOWN_CANDIDATE_GENE currently)
my %transcript2start;
# %transcript2gtex: key==$transcript, value==GTEX data to print (starting with \t)
my %transcript2gtex;
# also remember CHR and POS (the first POS we see for a variant affecting the
# transcript), replacing X Y M with 23-25 for easy sorting
my %transcript2chr;
my %transcript2coord;

# COHORT-SPECIFIC DATA: 
# %transcript2cohort2start: key==$transcript, value==hashref with key==$cohort,
# value==arrayref holding the start-of-line cohort-specific data for this transcript
# (one value per column), the common fields are undefined
my %transcript2cohort2start;
# %transcript2cohort2samples: key==$transcript, value==hashref with key==cohort,
# value is an arrayref with 6 hashrefs, one each for @countTypes,
# each hash has key==$grexome, value==number of variants (of that
# category), MODHIGH includes HIGH samples and MODER includes MODHIGH
# and HIGH samples,
# COMPHET also lists the samples that are HV (but these count as 2 variants)
my %transcript2cohort2samples;

# for each $cohort found in $inDir, store the header line to print
my %cohort2header;


while (my $inFile = readdir(INDIR)) {
    ($inFile =~ /^\./) && next;
    my $cohort;
    if ($inFile =~ (/^(\w+)\.filtered\.pick\.csv$/)) {
	$cohort = $1;
    }
    else {
	warn "W: cannot parse filename of inFile $inFile, skipping it\n";
	next;
    }

    my $now = strftime("%F %T", localtime);
    warn "I: $now - starting $0 on $cohort\n";

    open(INFILE, "$inDir/$inFile") ||
	die "E: cannot open infile $inDir/$inFile\n";

    ###################################
    # header line
    my $header = <INFILE>;
    chomp($header);
    my @headers = split(/\t/,$header);

    # $transCol == column of "Feature" (== transcript)
    my $transCol;
    # column of POSITION
    my $posCol;
    # columns of HV,HET
    my ($colHv,$colHet);
    # @destCols and @destColsSpec: for each column $i in infile: 
    # if $destCols[$i] >= 0 it is the column where that info goes in transcript2start
    # elsif $destCols[$i] == -1 source column is a GTEX
    # elsif $destColsSpec[$i] it is the column where it goes in transcript2cohort2start
    my @destCols = ();
    my @destColsSpec = ();

    # new headers
    my $newHeaders = join("\t",@keptColumns);
    # add COUNTs and SAMPLES
    foreach my $ct (@countTypes) {
	$newHeaders .= "\tCOUNT_$cohort"."_$ct";
    }
    # separator column
    $newHeaders .= "\tCOMPAT";
    foreach my $ct (@countTypes) {
	$newHeaders .= "\tCOUNT_COMPAT_$ct";
    }
    # separator column
    $newHeaders .= "\tNEGCTRL";
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
	elsif (defined $keptColsSpecific{$headers[$hi]}) {
	    $destColsSpec[$hi] = $keptColsSpecific{$headers[$hi]};
	}
	elsif ($headers[$hi] eq "POSITION") {
	    $posCol = $hi;
	}
	elsif ($headers[$hi] eq "HV") {
	    $colHv = $hi;
	}
	elsif ($headers[$hi] eq "HET") {
	    $colHet = $hi;
	}
	# else ignore this column
    }

    # add filter params at the end and store
    $cohort2header{$cohort} = "$newHeaders\t$headers[$#headers]";

    ###################################
    # body lines
    while(my $line = <INFILE>) {
	chomp($line);
	my @fields= split(/\t/,$line,-1);

	# LOW and MODIFIER variants are ignored
	($impact eq "LOW") && next;
	($impact eq "MODIFIER") && next;

	my $transcript = $fields[$transCol];
	if (! $transcript2start{$transcript}) {
	    # first time we see $transcript, fill chr, coord, start and gtex
	    ($fields[$posCol] =~ /^chr([^:]+):(\d+)$/) ||
		die "E: cannot grab chrom:pos in line:\n$line\n";
	    my ($chr,$coord) = ($1,$2);
	    # for sorting we want just the chrom number, replace X Y M by 23-25
	    if ($chr eq "X") { $transcript2chr{$transcript} = "23" }
	    elsif ($chr eq "Y") { $transcript2chr{$transcript} = "24" }
	    elsif ($chr eq "M") { $transcript2chr{$transcript} = "25" }
	    else { $transcript2chr{$transcript} = $chr }
	    $transcript2coord{$transcript} = $coord;

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
	    $transcript2start{$transcript} = \@start;
	    $transcript2gtex{$transcript} = $gtex;
	    # also initialize transcript2cohort2samples and %transcript2cohort2start
	    # with empty hashrefs
	    $transcript2cohort2samples{$transcript} = {};
	    $transcript2cohort2start{$transcript} = {};
	}

	if (! $transcript2cohort2start{$transcript}->{$cohort}) {
	    # first time we see this transcript in this cohort
	    my @startSpec = ();
	    foreach my $fi (0..$#fields) {
		if ($destColsSpec[$fi]) {
		    $startSpec[$destColsSpec[$fi]] = $fields[$fi];
		}
	    }
	    $transcript2cohort2start{$transcript}->{$cohort} = \@startSpec;

	    # initialize samples: arrayref with 6 hashrefs to empty hashes
	    $transcript2cohort2samples{$transcript}->{$cohort} = [{},{},{},{},{},{}];
	}

	# in any case, update %transcript2cohort2samples
	# HV and HET columns are processed very similarly, use flag to know
	# where we are: $isHV == 0 if HET, 1 if HV, init to -1 because we always 
	# increment $isHV immediately in loop
	my $isHV = -1;
	# Similarly we pocess the various IMPACTs homogeneously with $impactStart:
	# $impactStart == 0 if impact is HIGH (need to update counters HIGH, MODHIGH and MODER),
	# 1 if it's MODHIGH (update MODHIGH and MODER), 2 if it's MODER (just update MODER)
	my $impactStart;
	if ($impact eq "HIGH") { $impactStart = 0 }
	elsif ($impact eq "MODHIGH") { $impactStart = 1 }
	elsif ($impact eq "MODERATE") { $impactStart = 2 }
	else {
	    die "E: impact is $impact, should have been skipped earlier. Line:\n$line\n";
	}

	foreach my $col ($colHet,$colHv) {
	    # increment immediately, even if no samples
	    $isHV++;
	    my $samples = $fields[$col];
	    ($samples) || next;
	    # remove genotype, we only want the grexome IDs
	    $samples =~ s/^\d+\/\d+~//;
	    
	    foreach my $sample (split(/,/,$samples)) {
		# ignore [DP:AF]
		($sample =~ /^(grexome\d\d\d\d)/) || 
		    die "E: cannot grab grexome from sample $sample\n";
		my $grexome = $1;

		# always initialize to zero if needed before incrementing
		if ($isHV) {
		    # COUNT_HV_HIGH, COUNT_HV_MODHIGH and COUNT_HV_MODER may get +1 
		    # (depending on $impactStart)
		    foreach my $ai ($impactStart..2) {
			($transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$grexome}) ||
			    ($transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$grexome} = 0);
			$transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$grexome}++;
		    }
		}
		# whether $isHV or not, COUNT_COMPHET_HIGH COUNT_COMPHET_MODHIGH and
		# COUNT_COMPHET_MODER may get updated, but HV variants count as 2
		foreach my $ai (3+$impactStart..5) {
		    ($transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$grexome}) ||
			($transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$grexome} = 0);
		    $transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$grexome} += (1+$isHV);
		}
	    }
	}
    }

    $now = strftime("%F %T", localtime);
    warn "I: $now - Finished parsing $cohort infile\n";
    close(INFILE);
}
closedir(INDIR);


#########################################################
# now print all results

# build the sorted list of transcripts once and for all:
# sort by chrom then by coord
my @transcripts = sort {($transcript2chr{$a} <=> $transcript2chr{$b}) || ($transcript2coord{$a} <=> $transcript2coord{$b})} keys(%transcript2chr);


# for COMPHET counts and lists we only want samples with at least 2 variants
foreach my $transcript (keys %transcript2cohort2samples) {
    foreach my $cohort (keys %{$transcript2cohort2samples{$transcript}}) {
	foreach my $t2si (3..5) {
	    foreach my $sample (keys %{$transcript2cohort2samples{$transcript}->{$cohort}->[$t2si]}) {
		($transcript2cohort2samples{$transcript}->{$cohort}->[$t2si]->{$sample} >= 2) ||
		    (delete $transcript2cohort2samples{$transcript}->{$cohort}->[$t2si]->{$sample});
	    }
	}
    }
}

# output filehandles, one per key==cohort
my %outFHs;
foreach my $cohort (sort keys(%cohort2header)) {
    my $outFile = "$cohort.Transcripts.csv" ;
    open(my $fh, "> $outDir/$outFile") ||
	die "E cannot open outfile $outDir/$outFile: $!\n";
    $outFHs{$cohort} = $fh;
    print $fh $cohort2header{$cohort}."\n";
}

# print data
foreach my $transcript (@transcripts) {
    # need to store lineStarts for each cohort, so the "don't print redundantly" feature works
    # key==$cohort, value==string to print
    my %toPrint;
    foreach my $cohort (keys(%{$transcript2cohort2start{$transcript}})) {
	# fill the cohort-specific fields in $transcript2start{$transcript}
	foreach my $i (0..$#{$transcript2cohort2start{$transcript}->{$cohort}}) {
	    (defined $transcript2cohort2start{$transcript}->{$cohort}->[$i]) || next;
	    $transcript2start{$transcript}->[$i] = $transcript2cohort2start{$transcript}->{$cohort}->[$i];
	}
	my $toPrint = join("\t",@{$transcript2start{$transcript}});

	# COUNT_* values: 6 for $cohort, "COMPAT", then 6 for COMPATs, "NEGCTRL", and finally 6 for NEGCTRLS
	my @counts = (0) x 20;
	$counts[6] = "COMPAT";
	$counts[13] = "NEGCTRL";
	foreach my $thisCohort (keys(%cohort2header)) {
	    # $indexInCounts is 0 ($cohort), 7 (a compatible cohort), or 14 (a negctrl cohort)
	    my $indexInCounts = 14;
	    if ($compatible{$cohort}->{$thisCohort}) {
		$indexInCounts = 7;
	    }
	    elsif ($cohort eq $thisCohort) {
		$indexInCounts = 0;
	    }
	    # else $thisCohort is a NEGCTRL, use the default == 12

	    foreach my $i (0..5) {
		$counts[$indexInCounts + $i] += scalar(keys(%{$transcript2cohort2samples{$transcript}->{$thisCohort}->[$i]}));
	    }
	}
	$toPrint .= "\t".join("\t",@counts);
	
	$toPrint{$cohort} = $toPrint;
    }

    foreach my $cohort (keys(%toPrint)) {
	# lists of samples from this cohort: don't print redundantly.
	# Doing this now so the counts still include all samples satisfying the condition,
	# I just want to avoid redundancy in the lists of samples
	foreach my $t (keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[0]})) {
	    # HV_HIGH -> remove from all others
	    foreach my $i (1..5) {
		delete($transcript2cohort2samples{$transcript}->{$cohort}->[$i]->{$t});
	    }
	}
	foreach my $t (keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[1]})) {
	    # HV_MODHIGH -> remove from HV_MODER, COMPHET_MODHIGH and COMPHET_MODER
	    foreach my $i (2,4,5) {
		delete($transcript2cohort2samples{$transcript}->{$cohort}->[$i]->{$t});
	    }
	}
	foreach my $t (keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[2]})) {
	    # HV_MODER -> remove from COMPHET_MODER
	    delete($transcript2cohort2samples{$transcript}->{$cohort}->[5]->{$t});
	}
	foreach my $t (keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[3]})) {
	    # COMPHET_HIGH -> remove from COMPHET_MODHIGH and COMPHET_MODER
	    foreach my $i (4,5) {
		delete($transcript2cohort2samples{$transcript}->{$cohort}->[$i]->{$t});
	    }
	}
	foreach my $t (keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[4]})) {
	    # COMPHET_MODHIGH -> remove from COMPHET_MODER
	    delete($transcript2cohort2samples{$transcript}->{$cohort}->[5]->{$t});
	}

	# we will actually print line except if all 6 sample lists are empty
	my $printOK = 0;
	my $toPrintEnd = "";
	foreach my $i (0..5) {
	    $toPrintEnd .= "\t".join(',',sort(keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[$i]})));
	    (scalar(keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[$i]})) != 0) && ($printOK = 1);
	}
	
	# transcript2gtex already starts with \t
	$toPrintEnd .= $transcript2gtex{$transcript};

	($printOK) && (print {$outFHs{$cohort}} $toPrint{$cohort}."$toPrintEnd\n");
    }
}

my $now = strftime("%F %T", localtime);
warn "I: $now - $0 all done!\n";
