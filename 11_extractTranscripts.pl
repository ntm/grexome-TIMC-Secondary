#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2025
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


# 19/08/2019, but starting from extractCohorts.pl
# NTM

# Takes as arguments: $inDir $outDir [$pathologies]
# - $inDir must contain cohort TSVs as produced by extractCohorts.pl,
#   possibly filtered/reordered by filterVariants.pl and reorderColumns.pl;
# - $outDir doesn't exist, it will be created and filled with one TSV
#   per infile, adding .Transcripts to the name;
# - $pathologies (optional) is the pathologies metadata file.
#
# We then produce one TSV for each cohort.
# In each TSV we print one line for each transcript (=="Feature"), with:
# - COUNTSAMPLES_HV_HIGH+ = number of distinct samples with at least
#   one HV HIGH variant
# - COUNTSAMPLES_HV_MODHIGH+ = number of distinct samples with at
#   least one HV MODHIGH-or-HIGH variant
# - COUNTSAMPLES_HV_MODER+ = number of distinct samples with at
#   least one HV MODERATE-or-MODHIGH-or-HIGH variant
# - COUNTSAMPLES_BIALLELIC_HIGH+ = number of distinct samples with at
#   least TWO HET (or one HV) HIGH variants
# - COUNTSAMPLES_BIALLELIC_MODHIGH+ = number of distinct samples with at
#   least TWO HET (or one HV) MODHIGH-or-HIGH variants
# - COUNTSAMPLES_BIALLELIC_MODER+ = number of distinct samples with at
#   least TWO HET (or one HV) MODERATE-or-MODHIGH-or-HIGH variants
# - 1 more column COUNTSAMPLES_OTHERCAUSE with similar counts but counting
#   the samples with a "known causal variant" in another gene, all concatenated
#   into a single :-separated string (with :: between HVs and HETs)
# - another column COUNTSAMPLES_COMPAT same as OTHERCAUSE but counting
#   the samples belonging to compatible cohorts (as defined in $pathologies)
# - 6 final columns COUNTSAMPLES_NEGCTRL_* similar to $cohort counts but counting
#   the control samples (as in extractCohorts again).
# - HV_HIGH, HV_MODHIGH, HV_MODER, BIALLELIC_HIGH, BIALLELIC_MODHIGH, BIALLELIC_MODER:
#   non-redundant list of samples counted in the corresponding COUNTSAMPLES columns
#   (so eg the *MODHIGH columns don't list the corresponding *HIGH samples,
#   even though COUNTSAMPLES*MODHIGH+ counts them);
# - OTHERCAUSE_BIALLELIC_MODHIGH+, COMPAT_BIALLELIC_MODHIGH+, NEGCTRL_BIALLELIC_MODHIGH+:
#   all samples counted in corresponding COUNTSAMPLES columns (we don't list the
#   MODER samples).
# In addition, several useful columns from the source file are
# retained, see @keptColumns.
# A transcript line is not printed if all 6 COUNTSAMPLES_$cohort columns
# are zero.

use strict;
use warnings;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use Getopt::Long;
use POSIX qw(strftime);

use lib "$RealBin";
use grexome_metaParse qw(parsePathologies);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## hard-coded stuff that shouldn't change much

# columns we want to keep, in this order:
my @keptColumns = qw(SYMBOL KNOWN_CANDIDATE_GENE Feature CANONICAL BIOTYPE Gene MANE_SELECT);
# in addition we insert the new COUNTSAMPLES* columns right after the last @keptColumns
# and immediately followed by the HV_HIGH et al colums, and we then copy all 
# the GTEX_* columns (in the same order as in infile)

# among the @keptColumns some have cohort-specific data: list them
my @keptColumnsSpecific = qw(KNOWN_CANDIDATE_GENE);


# also for convenience: the types of samples to count
my @countTypes = ("HV_HIGH","HV_MODHIGH","HV_MODER","BIALLELIC_HIGH","BIALLELIC_MODHIGH","BIALLELIC_MODER");


#############################################
## options / params from the command-line

# inDir contains cohort TSVs, $outDir must not pre-exist, no defaults
my ($inDir, $outDir);

# path+file of the pathologies XLSX file
my $pathologies = "";

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "\nParse cohort TSVs from inDir, create transcript TSVs in outDir.
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--indir string: subdir containing cohort TSVs as produced by extractCohorts.pl, 
                              possibly filtered/reordered by 10_filterAndReorderAll.pl;
--outdir string: subdir where resulting Transcripts TSV files will be created, must not pre-exist
--pathologies [optional] : pathologies metadata xlsx file, with path
--help : print this USAGE";

GetOptions ("indir=s" => \$inDir,
	    "outdir=s" => \$outDir,
	    "pathologies=s" => \$pathologies,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($inDir) ||
    die "E $0: you must provide an indir.\nUSAGE:$USAGE\n";
(-d $inDir) ||
    die "E $0: inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "E $0: cannot opendir inDir $inDir\n";

($outDir) ||
    die "E $0: you must provide an outdir\n";
(-e $outDir) && 
    die "E $0: provided outDir $outDir already exists, remove it or choose another name.\n";
mkdir($outDir) || die "E $0: cannot mkdir outDir $outDir\n";

# If $pathologies was provided we want to parse (and check) it now, it is used
# to populate $compatibleR
# $compatibleR: hashref, key is a cohort name, value is a hashref
# with keys == cohorts that are compatible with this cohort, value==1
my $compatibleR;
if ($pathologies) {
    (-f $pathologies) || die "E $0: the supplied pathologies file $pathologies doesn't exist\n";
    $compatibleR = &parsePathologies($pathologies);
}
# else $compatibleR stays undef

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


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
	die "E $0: column $col is in keptColsSpecific but not in keptColumns, add it there\n";
    $keptColsSpecific{$col} = $keptCols{$col};
    delete($keptCols{$col});
}

#########################################################

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
# also remember CHR and POS (the smallest POS we see in any cohort for a variant
#  affecting the transcript), replacing X Y M with 1023-1025 for easy sorting
my %transcript2chr;
my %transcript2coord;

# COHORT-SPECIFIC DATA: 
# %transcript2cohort2start: key==$transcript, value==hashref with key==$cohort,
# value==arrayref holding the start-of-line cohort-specific data for this transcript
# (one value per column), the common fields are undefined
my %transcript2cohort2start;
# %transcript2cohort2samples: key==$transcript, value==hashref with key==cohort,
# value is an arrayref with 6 hashrefs, one each for @countTypes,
# each hash has key==$sample, value==number of variants (of that
# category), MODHIGH includes HIGH samples and MODER includes MODHIGH
# and HIGH samples,
# BIALLELIC also lists the samples that are HV (but these count as 2 variants)
my %transcript2cohort2samples;
# %transcript2cohort2samplesOC: same as %transcript2cohort2samples but counting
# OTHERCAUSE variants
my %transcript2cohort2samplesOC;

# for each $cohort found in $inDir, store the header line to print
my %cohort2header;


while (my $inFile = readdir(INDIR)) {
    ($inFile =~ /^\./) && next;
    my $cohort;
    if ($inFile =~ /^(\w+)\..*csv.gz$/) {
	$cohort = $1;
    }
    else {
	warn "W $0: cannot parse filename of inFile $inFile, skipping it\n";
	next;
    }

    open(INFILE, "gunzip -c $inDir/$inFile |") ||
	die "E $0: cannot gunzip-open infile $inDir/$inFile\n";

    ###################################
    # header line
    my $header = <INFILE>;
    ($header) || die "E: $0 - input file $inDir/$inFile is empty\n";
    chomp($header);
    my @headers = split(/\t/,$header);

    # $transCol == column of "Feature" (== transcript)
    my $transCol;
    # column of POSITION
    my $posCol;
    # column of IMPACT
    my $impactCol;
    # columns of HV,HET, OTHERCAUSE_HV, OTHERCAUSE_HET
    my ($hvCol,$hetCol,$hvColOC,$hetColOC) = (-1,-1,-1,-1);
    # @destCols and @destColsSpec: for each column $i in infile: 
    # if $destCols[$i] >= 0 it is the column where that info goes in transcript2start
    # elsif $destCols[$i] == -1 source column is a GTEX
    # elsif $destColsSpec[$i] it is the column where it goes in transcript2cohort2start
    my @destCols = ();
    my @destColsSpec = ();

    # new headers
    my $newHeaders = join("\t",@keptColumns);
    # COUNTSAMPLES
    foreach my $ct (@countTypes) {
	$newHeaders .= "\tCOUNTSAMPLES_$ct+";
    }
    # single columns with OTHERCAUSE and COMPAT counts
    $newHeaders .= "\tCOUNTSAMPLES_OTHERCAUSE";
    $newHeaders .= "\tCOUNTSAMPLES_COMPAT";
    # all columns for NEGCTRL so we can filter/sort
    foreach my $ct (@countTypes) {
	$newHeaders .= "\tCOUNTSAMPLES_NEGCTRL_$ct+";
    }
    # SAMPLES from $cohort (excluding OCs)
    foreach my $ct (@countTypes) {
	$newHeaders .= "\t$ct";
    }
    # other SAMPLES: OC, COMPAT, NEGCTRL -> only list BIALLELIC_MODHIGH cumulatively
    $newHeaders .= "\tOTHERCAUSE_BIALLELIC_MODHIGH+";
    $newHeaders .= "\tCOMPAT_BIALLELIC_MODHIGH+";
    $newHeaders .= "\tNEGCTRL_BIALLELIC_MODHIGH+";
    
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
	elsif ($headers[$hi] eq "IMPACT") {
	    $impactCol = $hi;
	}
	elsif ($headers[$hi] eq $cohort."_HV") {
	    $hvCol = $hi;
	}
	elsif ($headers[$hi] eq $cohort."_OTHERCAUSE_HV") {
	    $hvColOC = $hi;
	}
	elsif ($headers[$hi] eq $cohort."_HET") {
	    $hetCol = $hi;
	}
	elsif ($headers[$hi] eq $cohort."_OTHERCAUSE_HET") {
	    $hetColOC = $hi;
	}
	# else ignore this column
    }
    # sanity check
    (($hvCol >= 0) && ($hvColOC >= 0) && ($hetCol >= 0) && ($hetColOC >= 0)) || 
	die "E $0: couldn't find one of HV/HET/OCHV/OCHET for $cohort: $hvCol $hetCol $hvColOC $hetColOC\n";

    # add filter params at the end and store
    $cohort2header{$cohort} = "$newHeaders\t$headers[$#headers]";

    ###################################
    # body lines
    while(my $line = <INFILE>) {
	chomp($line);
	my @fields= split(/\t/,$line,-1);

	my $impact = $fields[$impactCol];
	# LOW and MODIFIER variants are ignored
	($impact eq "LOW") && next;
	($impact eq "MODIFIER") && next;

	# grab chrom and coord
	($fields[$posCol] =~ /^chr([^:]+):(\d+)$/) ||
	    ($fields[$posCol] =~ /^([^:]+):(\d+)$/) ||
	    die "E $0: cannot grab chrom:pos in line:\n$line\n";
	my ($chr,$coord) = ($1,$2);

	my $transcript = $fields[$transCol];
	if (! $transcript2start{$transcript}) {
	    # first time we see $transcript, fill chr, coord, start and gtex
	    # for sorting we want just the chrom number, replace X Y M/MT by 1023-1025
	    if ($chr eq "X") { $transcript2chr{$transcript} = "1023" }
	    elsif ($chr eq "Y") { $transcript2chr{$transcript} = "1024" }
	    elsif (($chr eq "M") || ($chr eq "MT")) { $transcript2chr{$transcript} = "1025" }
	    elsif ($chr =~ /^\d+$/) { $transcript2chr{$transcript} = $chr }
	    else { die "E $0: found unexpected chromosome $chr\n" }
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
	    # also initialize %transcript2cohort2samples, %transcript2cohort2samplesOC,
	    # and %transcript2cohort2start with empty hashrefs
	    $transcript2cohort2samples{$transcript} = {};
	    $transcript2cohort2samplesOC{$transcript} = {};
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
	    # same for OTHERCAUSE samples
	    $transcript2cohort2samplesOC{$transcript}->{$cohort} = [{},{},{},{},{},{}];
	}

	# update transcript2coord if new coord is smaller
	($coord < $transcript2coord{$transcript}) && ($transcript2coord{$transcript} = $coord);
	# in any case, update %transcript2cohort2samples and transcript2cohort2samplesOC
	# we will read the sample columns in this order:
	my @cols = ($hetCol,$hvCol,$hetColOC,$hvColOC);
	# HV and HET columns are processed very similarly, index in @cols must
	# be even for HETs and odd for HVs;
	# similarly we pocess the various IMPACTs homogeneously with $impactStart:
	# $impactStart == 0 if impact is HIGH (need to update counters HIGH, MODHIGH and MODER),
	# 1 if it's MODHIGH (update MODHIGH and MODER), 2 if it's MODER (just update MODER)
	my $impactStart;
	if ($impact eq "HIGH") { $impactStart = 0 }
	elsif ($impact eq "MODHIGH") { $impactStart = 1 }
	elsif ($impact eq "MODERATE") { $impactStart = 2 }
	else {
	    die "E $0: impact is $impact, should have been skipped earlier. Line:\n$line\n";
	}

	foreach my $coli (0..$#cols) {
	    # $isHV = 0 for HETs and 1 for HVs
	    my $isHV = $coli % 2;
	    # $isOC is 0 for undiagnosed $cohort samples and 1 for OTHERCAUSE samples
	    my $isOC = int($coli / 2);
	    my $samples = $fields[$cols[$coli]];
	    ($samples) || next;

	    # remove genotype, we only want the sample IDs
	    $samples =~ s/^\d+\/\d+~//;
	    
	    foreach my $sample (split(/,/,$samples)) {
		# ignore [DP:AF] / [GQ:FR:BP] and possibly patientIDs
		($sample =~ /^([^(\[]+)/) ||
		    die "E $0: cannot grab sampleID from sample $sample\n";
		my $sampleID = $1;

		# always initialize to zero if needed before incrementing
		if ($isHV) {
		    # COUNTSAMPLES_HV_HIGH, COUNTSAMPLES_HV_MODHIGH and COUNTSAMPLES_HV_MODER or their
		    # OTHERCAUSE counterparts may get +1 (depending on $impactStart)
		    foreach my $ai ($impactStart..2) {
			($transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$sampleID}) ||
			    ($transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$sampleID} = 0);
			($transcript2cohort2samplesOC{$transcript}->{$cohort}->[$ai]->{$sampleID}) ||
			    ($transcript2cohort2samplesOC{$transcript}->{$cohort}->[$ai]->{$sampleID} = 0);
			if ($isOC) {
			    $transcript2cohort2samplesOC{$transcript}->{$cohort}->[$ai]->{$sampleID}++;
			}
			else {
			    $transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$sampleID}++;
			}
		    }
		}
		# whether $isHV or not, COUNTSAMPLES_BIALLELIC_HIGH COUNTSAMPLES_BIALLELIC_MODHIGH and
		# COUNTSAMPLES_BIALLELIC_MODER (or their OC counterparts) may get updated, but
		# HV variants count as 2
		foreach my $ai (3+$impactStart..5) {
		    ($transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$sampleID}) ||
			($transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$sampleID} = 0);
		    ($transcript2cohort2samplesOC{$transcript}->{$cohort}->[$ai]->{$sampleID}) ||
			($transcript2cohort2samplesOC{$transcript}->{$cohort}->[$ai]->{$sampleID} = 0);
		    if ($isOC) {
			$transcript2cohort2samplesOC{$transcript}->{$cohort}->[$ai]->{$sampleID} += (1+$isHV);
		    }
		    else {
			$transcript2cohort2samples{$transcript}->{$cohort}->[$ai]->{$sampleID} += (1+$isHV);
		    }
		}
	    }
	}
    }
    close(INFILE);
}
closedir(INDIR);

#########################################################
# now print all results

# build the sorted list of transcripts once and for all:
# sort by chrom, then by coord (of the first variant affecting
# the transcript), then by transcript name...
my @transcripts = sort {($transcript2chr{$a} <=> $transcript2chr{$b}) || ($transcript2coord{$a} <=> $transcript2coord{$b}) || ($a cmp $b)} keys(%transcript2chr);


# for HV counts and lists we want samples with at least one variant (some
# keys got created at zero and stayed there, delete them);
# and for BIALLELIC counts and lists we only want samples with at least 2 variants.
# the OC version is defined for exactly the same transcripts and cohorts as the
# normal version, but samples can be different
foreach my $transcript (keys %transcript2cohort2samples) {
    foreach my $cohort (keys %{$transcript2cohort2samples{$transcript}}) {
	# HVs
	foreach my $t2si (0..2) {
	    foreach my $sample (keys %{$transcript2cohort2samples{$transcript}->{$cohort}->[$t2si]}) {
		($transcript2cohort2samples{$transcript}->{$cohort}->[$t2si]->{$sample} >= 1) ||
		    (delete $transcript2cohort2samples{$transcript}->{$cohort}->[$t2si]->{$sample});
	    }
	    foreach my $sample (keys %{$transcript2cohort2samplesOC{$transcript}->{$cohort}->[$t2si]}) {
		($transcript2cohort2samplesOC{$transcript}->{$cohort}->[$t2si]->{$sample} >= 1) ||
		    (delete $transcript2cohort2samplesOC{$transcript}->{$cohort}->[$t2si]->{$sample});
	    }
	}
	# BIALLELICs
	foreach my $t2si (3..5) {
	    foreach my $sample (keys %{$transcript2cohort2samples{$transcript}->{$cohort}->[$t2si]}) {
		($transcript2cohort2samples{$transcript}->{$cohort}->[$t2si]->{$sample} >= 2) ||
		    (delete $transcript2cohort2samples{$transcript}->{$cohort}->[$t2si]->{$sample});
	    }
	    foreach my $sample (keys %{$transcript2cohort2samplesOC{$transcript}->{$cohort}->[$t2si]}) {
		($transcript2cohort2samplesOC{$transcript}->{$cohort}->[$t2si]->{$sample} >= 2) ||
		    (delete $transcript2cohort2samplesOC{$transcript}->{$cohort}->[$t2si]->{$sample});
	    }
	}
    }
}

# output filehandles, one per key==cohort
my %outFHs;
foreach my $cohort (sort keys(%cohort2header)) {
    my $outFile = "$cohort.Transcripts.csv" ;
    open(my $fh, "> $outDir/$outFile") ||
	die "E $0: cannot open outfile $outDir/$outFile: $!\n";
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
	    # we modify the cohort-independant transcript2start data, but it's OK
	    # because if a field is defined in one cohort it must be defined in
	    # every cohort (possibly with value "", which is what we want)
	    $transcript2start{$transcript}->[$i] = $transcript2cohort2start{$transcript}->{$cohort}->[$i];
	}
	my $toPrint = join("\t",@{$transcript2start{$transcript}});

	# COUNTSAMPLES_* values: 6 each for $cohort, OCs, COMPATs, NEGCTRLs
	my @counts = (0) x 24;
	foreach my $thisCohort (keys(%{$transcript2cohort2samples{$transcript}})) {
	    if ($cohort eq $thisCohort) {
		# update non-OC $cohort...
		foreach my $i (0..5) {
		    $counts[$i] += scalar(keys(%{$transcript2cohort2samples{$transcript}->{$thisCohort}->[$i]}));
		}
		# and also $cohort_OC, at indexes 6..11
		foreach my $i (0..5) {
		    $counts[6 + $i] += scalar(keys(%{$transcript2cohort2samplesOC{$transcript}->{$thisCohort}->[$i]}));
		}
	    }
	    else {
		# for COMPATs and NEGCTRLs we need to add non-OC and OC counts
		# $indexInCounts is 12 if $thisCohort is compatible with $cohort and 18 if it's a NEGCTRL
		my $indexInCounts = 18;
		if ($compatibleR->{$cohort}->{$thisCohort}) {
		    $indexInCounts = 12;
		}
		foreach my $i (0..5) {
		    $counts[$indexInCounts + $i] += scalar(keys(%{$transcript2cohort2samples{$transcript}->{$thisCohort}->[$i]}));
		    $counts[$indexInCounts + $i] += scalar(keys(%{$transcript2cohort2samplesOC{$transcript}->{$thisCohort}->[$i]}));
		}
	    }
	}
	# print $cohort counts individually
	foreach my $i (0..5) {
	    $toPrint .= "\t$counts[$i]";
	}
	# then OC and COMPAT as a ':'-separated single string with '::' between HV and HET
	foreach my $i (6,12) {
	    $toPrint .= "\t".join(':',@counts[$i..$i+2])."::".join(':',@counts[$i+3..$i+5]);
	}
	# anf finally NEGCTRLs individually
	foreach my $i (0..5) {
	    $toPrint .= "\t$counts[18+$i]";
	}

	# $toPrint{$cohort} stays undef for this cohort if it doesn't at least have a
	# BIALLELIC_MODER sample, otherwise save $toPrint
	if ($counts[5] > 0) {
	    $toPrint{$cohort} = $toPrint;
	}
    }
    
    # done counting, now remove redundancy (eg between HV_MODHIGH and HV_HIGH) in lists of
    # non-OTHERCAUSE samples from every cohort
    foreach my $cohort (keys(%{$transcript2cohort2start{$transcript}})) {
	foreach my $t (keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[0]})) {
	    # HV_HIGH -> remove from all others
	    foreach my $i (1..5) {
		delete($transcript2cohort2samples{$transcript}->{$cohort}->[$i]->{$t});
	    }
	}
	foreach my $t (keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[1]})) {
	    # HV_MODHIGH -> remove from HV_MODER, BIALLELIC_MODHIGH and BIALLELIC_MODER
	    foreach my $i (2,4,5) {
		delete($transcript2cohort2samples{$transcript}->{$cohort}->[$i]->{$t});
	    }
	}
	foreach my $t (keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[2]})) {
	    # HV_MODER -> remove from BIALLELIC_MODER
	    delete($transcript2cohort2samples{$transcript}->{$cohort}->[5]->{$t});
	}
	foreach my $t (keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[3]})) {
	    # BIALLELIC_HIGH -> remove from BIALLELIC_MODHIGH and BIALLELIC_MODER
	    foreach my $i (4,5) {
		delete($transcript2cohort2samples{$transcript}->{$cohort}->[$i]->{$t});
	    }
	}
	foreach my $t (keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[4]})) {
	    # BIALLELIC_MODHIGH -> remove from BIALLELIC_MODER
	    delete($transcript2cohort2samples{$transcript}->{$cohort}->[5]->{$t});
	}
    }

    # now finish building $toPrint, adding the lists of samples and GTEX data
    foreach my $cohort (keys(%toPrint)) {
	my $toPrint = $toPrint{$cohort};
	# the 6 cohort samplelists
	foreach my $i (0..5) {
	    $toPrint .= "\t".join(',',sort(keys(%{$transcript2cohort2samples{$transcript}->{$cohort}->[$i]})));
	}
	
	#  OTHERCAUSE_BIALLELIC_MODHIGH is complete, just copy it (BIALLELIC_MODHIGH -> index 4)
	$toPrint .= "\t".join(',',sort(keys(%{$transcript2cohort2samplesOC{$transcript}->{$cohort}->[4]})));
	
	# COMPAT_BIALLELIC_MODHIGH, NEGCTRL_BIALLELIC_MODHIGH: need to examine all other cohorts
	my @compat = ();
	my @negctrl = ();
	foreach my $thisCohort (keys(%cohort2header)) {
	    ($cohort eq $thisCohort) && next;
	    # otherwise we want all BIALLELIC_MODHIGH (or better) samples from $thisCohort,
	    # including the OTHERCAUSE ones
	    my @goodSamples;
	    # build full $cohort_BIALLELIC_MODHIGH+: need indexes 0,1,3,4
	    foreach my $i (0,1,3,4) {
		push(@goodSamples, keys(%{$transcript2cohort2samples{$transcript}->{$thisCohort}->[$i]}));
	    }
	    # add OCs (it's complete, just grab list at index 4)
	    push(@goodSamples, keys(%{$transcript2cohort2samplesOC{$transcript}->{$thisCohort}->[4]}));

	    if ($compatibleR->{$cohort}->{$thisCohort}) {
		push(@compat, @goodSamples); 
	    }
	    else {
		push(@negctrl, @goodSamples); 
	    }
	}
	$toPrint .= "\t".join(',',sort(@compat));
	$toPrint .= "\t".join(',',sort(@negctrl));

	# transcript2gtex already starts with \t
	$toPrint .= $transcript2gtex{$transcript};

	print {$outFHs{$cohort}} "$toPrint\n";
    }
}

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";
