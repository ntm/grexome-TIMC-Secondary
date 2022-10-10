#!/usr/bin/perl

# 25/03/2018
# NTM

# Parse COHORT TSVs and produce SAMPLE TSVs. See $USAGE.

use strict;
use warnings;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use Getopt::Long;
use POSIX qw(strftime);
use Parallel::ForkManager;

use lib "$RealBin";
use grexome_metaParse qw(parseSamples);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## options / params from the command-line

# samples metadata xlsx
my $samplesFile = "";

# $inDir containing cohort TSVs
my $inDir = "";

# $outDir will be created and filled with one TSV per sample
my $outDir = "";

# number of cohorts to process in parallel
my $jobs = 8;

# $covDir is optional, if provided it is a subdir containing per-sample 
# coverage files as produced by 0_coverage.pl
my $covDir = "";

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "\nParse COHORT TSV files in inDir, and produce SAMPLE TSVs in outDir.
Filenames in outDir will include patientID/specimenID.
If covDir is provided, the global coverage data (ALL_CANDIDATES and ALL_SAMPLED)
for each sample is grabbed from covDir and appended at the end of the header line.
For a sample, we only print lines from its cohort file and where it has an HV or 
HET genotype: this genotype is printed in new columns GENOTYPE and DP:AF/BF:RR,
inserted right after KNOWN_CANDIDATE_GENE.
Immediately after DP:AF/BF:RR we insert a new column BIALLELIC, value is one of:
  HIGH -> patient has >=2 HET or at >=1 HV HIGH variants;
  MODHIGH -> patient has >=2 HET or >=1 HV variants of impact HIGH or MODHIGH,
     but isn't in HIGH category;
  MODERATE -> same for impact >= MODERATE;
  LOW -> same for impact >= LOW;
  NO -> patient has at most one allele >= LOW.
Also the genoData columns (\$cohort_HV, NEGCTRL_HV, etc...) are not printed.

Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--samples : samples metadata xlsx file, with path
--indir : must contain cohort TSVs as produced by extractCohorts.pl,
          possibly filtered and reordered with 7_filterAndReorderAll.pl, and 
          possibly gzipped (but not with PatientIDs)
--outdir : subdir where SAMPLE TSVs will be created, must not pre-exist
--covdir : optional, if provided it must be a subdir containing per-sample 
           coverage files as produced by 0_coverage.pl
--jobs [$jobs] : number of cohorts to process in parallel
--help : print this USAGE";


GetOptions ("samples=s" => \$samplesFile,
	    "indir=s" => \$inDir,
	    "outdir=s" => \$outDir,
	    "covdir=s" => \$covDir,
	    "jobs=i" => \$jobs,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($samplesFile) || die "E $0: you must provide a samples file\n";
(-f $samplesFile) || die "E $0: the supplied samples file doesn't exist\n";

($inDir) || die "E $0: you must provide an inDir\n";
(-d $inDir) ||
    die "E: $0 - inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "E: $0 - cannot opendir inDir $inDir\n";

if ($covDir) {
    (-d $covDir) ||
	die "E: $0 - covDir $covDir was provided but it doesn't exist or isn't a directory\n";
}
else {
    warn "I: $0 - no covDir provided, coverage statistics won't be appended to the header lines\n";
}

($jobs > 0) || die "E $0: jobs = $jobs, really??\n";

($outDir) || die "E $0: you must provide an outDir\n";
(-e $outDir) && 
    die "E $0: outdir $outDir already exists, remove it or choose another name.\n";
mkdir($outDir) || die "E: $0 - cannot mkdir outDir $outDir\n";

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


#########################################################
# parse useful info from samplesFile

# key==pathologyID, value is an arrayref of all samples with this patho
my %cohort2samples = ();

# key==sample, value is patientID if it exists, specimenID otherwise
my $sample2patientR;

{
    my @parsed = &parseSamples($samplesFile);
    my $sample2pathoR = $parsed[0];
    $sample2patientR = $parsed[2];

    foreach my $sample (sort keys %$sample2pathoR) {
	my $patho = $sample2pathoR->{$sample};
	(defined $cohort2samples{$patho}) || ($cohort2samples{$patho} = []);
	push(@{$cohort2samples{$patho}}, $sample);
    }
}

#########################################################
# read infiles

my $pm = new Parallel::ForkManager($jobs);

while (my $inFile = readdir(INDIR)) {
    ($inFile =~ /^\./) && next;
    $pm->start && next;
    my ($cohort,$fileEnd,$gz);
    if ($inFile =~ (/^([^\.]+)\.(.*csv)$/)) {
	# $fileEnd allows for .canon etc...
	($cohort,$fileEnd) = ($1,$2);
    }
    elsif ($inFile =~ (/^([^\.]+)\.(.*csv)\.gz$/)) {
	($cohort,$fileEnd) = ($1,$2);
	$gz = 1;
    }
    else {
	warn "W $0: cannot parse filename of inFile $inDir/$inFile, skipping it\n";
	$pm->finish;
    }

    my $inFull = "$inDir/$inFile";
    ($gz) && ($inFull = "gunzip -c $inFull | ");
    open(IN, $inFull) ||
	die "E: $0 - cannot (gunzip-?)open cohort datafile $inDir/$inFile (as $inFull)\n";
    my $header = <IN>;
    ($header) || die "E: $0 - input file $inDir/$inFile is empty\n";
    chomp($header);
    my @header = split(/\t/,$header);

    # need $cohort_HV, $cohort_HET, $cohort_OTHERCAUSE_HV, $cohort_OTHERCAUSE_HET
    # and also other genoData columns (will all be removed)
    my ($hvCol,$hetCol,$hvColOC,$hetColOC) = (-1,-1,-1,-1);
    # $colsToRemove[$i] == 1 if column $i in infile must be removed
    my @colsToRemove;
    # reverse so we can splice out elements
    foreach my $i (reverse(0..$#header)) {
	# $toRemove: boolean, true iff column must be removed
	my $toRemove = 0;
	if ($header[$i] eq $cohort."_HV") {
	    $hvCol = $i;
	    $toRemove = 1;
	}
	elsif ($header[$i] eq $cohort."_OTHERCAUSE_HV") {
	    $hvColOC = $i;
	    $toRemove = 1;
	}
	elsif ($header[$i] eq $cohort."_HET") {
	    $hetCol = $i;
	    $toRemove = 1;
	}
	elsif ($header[$i] eq $cohort."_OTHERCAUSE_HET") {
	    $hetColOC = $i;
	    $toRemove = 1;
	}
	elsif (grep(/^$header[$i]$/, ("COMPAT_HV","COMPAT_HET","NEGCTRL_HV","NEGCTRL_HET","OTHERGENO"))) {
	    $toRemove = 1;
	}

	if ($toRemove) {
	    $colsToRemove[$i] = 1;
	    splice(@header,$i,1);
	}
    }
    (($hvCol >= 0) && ($hvColOC >= 0)) || 
	die "E: $0 - couldn't find ${cohort}_HV or ${cohort}_OTHERCAUSE_HV in header of infile $inFile\n";
    (($hetCol >= 0) && ($hetColOC >= 0)) || 
	die "E: $0 - couldn't find HET or OC_HET in header of infile $inFile\n";

    # KNOWN_CANDIDATE_GENE, Feature and IMPACT column indexes after
    # removing @colsToRemove columns
    my ($knownCandidateCol,$featureCol, $impactCol) = (-1,-1,-1);
 
    foreach my $i (0..$#header) {
	if ($header[$i] eq "KNOWN_CANDIDATE_GENE") {
	    $knownCandidateCol = $i;
	    $header[$i] .= "\tGENOTYPE";
	    $header[$i] .= "\tDP:AF/BF:RR";
	    $header[$i] .= "\tBIALLELIC";
	}
	elsif ($header[$i] eq "Feature") {
	    $featureCol = $i;
	}
	elsif ($header[$i] eq "IMPACT") {
	    $impactCol = $i;
	}
    }
    ($knownCandidateCol >= 0) || 
	die "E $0: couldn't find KNOWN_CANDIDATE_GENE in header of infile $inFile\n";
    ($featureCol >= 0) || 
	die "E $0: couldn't find Feature in header of infile $inFile\n";
    ($impactCol >= 0) || 
	die "E $0: couldn't find IMPACT in header of infile $inFile\n";

    $header = join("\t",@header);

    # hash of filehandles open for writing, one for each sample
    # from this cohort
    # will be gzipped if infiles were
    my %outFHs;

    ($cohort2samples{$cohort}) || 
	die "E: $0 - cohort $cohort parsed from filename of infile $inFile is not in $samplesFile\n";
    foreach my $sample (@{$cohort2samples{$cohort}}) {
	my $patient = $sample2patientR->{$sample};
	my $outFile = "$outDir/$cohort.$sample.$patient.$fileEnd";
	($gz) && ($outFile .= ".gz");
	my $outFull = " > $outFile";
	($gz) && ($outFull = " | gzip -c $outFull");
	open (my $FH, $outFull) || die "E: $0 - cannot (gzip-?)open $outFile for writing (as $outFull)\n";

	# grab global coverage data for $sample if a covDir was provided
	my $headerCov = "";
	if ($covDir) {
	    my $covFile = "$covDir/coverage_$sample.csv";
	    (-f $covFile) ||
		die "E $0: covDir provided but cannot find coverage data for $sample: $covFile\n";

	    # global coverage data is in last 2 lines
	    open(COV, "tail -n 2 $covFile |") ||
		die "E: $0 - cannot tail-grab cov data from covFile with: tail -n 2 $covFile\n";
	    # ALL_CANDIDATES
	    my $covLine = <COV>;
	    chomp $covLine;
	    my @covFields = split(/\t/,$covLine);
	    (@covFields == 8) || die "E $0: expecting 8 fields from candidates coverage line $covLine\n";
	    $headerCov = "\tCoverage_Candidates_50x=$covFields[5] Coverage_Candidates_20x=$covFields[6] Coverage_Candidates_10x=$covFields[7]";
	    # ALL_GENES
	    $covLine = <COV>;
	    chomp $covLine;
	    @covFields = split(/\t/,$covLine);
	    (@covFields == 8) || die "E $0: expecting 8 fields from sampled coverage line $covLine\n";
	    $headerCov .= "   Coverage_AllGenes_50x=$covFields[5] Coverage_AllGenes_20x=$covFields[6] Coverage_AllGenes_10x=$covFields[7]";
	    close(COV);
	}
	print $FH "$header$headerCov\n";
	$outFHs{$sample} = $FH ;
    }

    # now read the data
    # in order to print NB_* we need:
    # key == sample, values are arrayrefs of the beginnings and ends of lines
    # (respectively) that must be printed for this sample
    my %sample2lineStarts;
    my %sample2lineEnds;
    # key == sample, value is an arrayref, for each line to print it holds
    # the transcript that this line deals with
    my %sample2transcripts;
    # key == sample, value is a hashref whose keys are transcripts and values
    # are arrayrefs with 4 ints:
    # numbers of HIGH, MODHIGH, MODER, and LOW variant alleles found for this
    # transcript in this sample (so, an HV counts as 2 and a HET as 1)
    my %sample2trans2counters;

    while (my $line = <IN>) {
	chomp($line);
	my @fields = split(/\t/, $line, -1) ;
	# grab needed GENO data: rip out the actual geno (eg 1/1, we will just use HV
	# or HET for geno, the actual allele is in ALLELE_NUM), and store the
	# comma-separated lists of samples, HV and HVOC at index 0, HET and HETOC at index 1
	my @genoData = ("","");
	foreach my $fi ($hvCol,$hvColOC) {
	    my $gd = $fields[$fi];
	    if ($gd) {
		($gd =~ s/^[^~]+~//) || 
		    die "E: $0 - cannot rip out geno from $gd, infile $inFile\n";
		($genoData[0]) && ($genoData[0] .= ",");
		$genoData[0] .= $gd;
	    }
	}
	foreach my $fi ($hetCol,$hetColOC) {
	    my $gd = $fields[$fi];
	    if ($gd) {
		($gd =~ s/^[^~]+~//) || 
		    die "E: $0 - cannot rip out geno from $gd, infile $inFile\n";
		($genoData[1]) && ($genoData[1] .= ",");
		$genoData[1] .= $gd;
	    }
	}

	# splice all genoData columns out
	foreach my $i (reverse(0..$#colsToRemove)) {
	    ($colsToRemove[$i]) && splice(@fields,$i,1);
	}

	my $toPrintStart = join("\t",@fields[0..$knownCandidateCol])."\t";
	my $toPrintEnd = join("\t",@fields[($knownCandidateCol+1)..$#fields])."\n";
	my $transcript = $fields[$featureCol];
	my $impact = $fields[$impactCol];

	foreach my $i (0,1) {
	    if ($genoData[$i]) {
		my $geno;
		($i == 0) && ($geno = "HV");
		($i == 1) && ($geno = "HET");
		foreach my $sampleData (split(/,/,$genoData[$i])) {
		    # grab sample and [DP:AF] / [BF:RR], we know it must be there
		    # (allowing AF > 1 for Strelka bug)
		    ($sampleData =~ /^([^\[\s]+)\[(\d+:\d+\.\d\d?)\]$/) ||
			die  "E $0: inFile $inFile has a sampleData (in a genoData) that I can't parse: $sampleData\n";
		    my ($sample,$dpaf) = ($1,$2);

		    # initialize everything for this sample if needed
		    ($sample2lineStarts{$sample}) || ($sample2lineStarts{$sample} = []);
		    ($sample2lineEnds{$sample}) || ($sample2lineEnds{$sample} = []);
		    ($sample2transcripts{$sample}) || ($sample2transcripts{$sample} = []);
		    ($sample2trans2counters{$sample}) || ($sample2trans2counters{$sample} = {});
		    ($sample2trans2counters{$sample}->{$transcript}) || ($sample2trans2counters{$sample}->{$transcript} = [0,0,0,0]);

		    # now fill our data structures
		    push(@{$sample2lineStarts{$sample}}, "$toPrintStart$geno\t$dpaf");
		    push(@{$sample2lineEnds{$sample}}, "\t$toPrintEnd");
		    push(@{$sample2transcripts{$sample}}, $transcript);

		    if ($impact eq "HIGH") {
			# 2-$i is 1 for HET and 2 for HV...
			$sample2trans2counters{$sample}->{$transcript}->[0] += 2-$i;
		    }
		    elsif ($impact eq "MODHIGH") {
			$sample2trans2counters{$sample}->{$transcript}->[1] += 2-$i;
		    }
		    elsif ($impact eq "MODERATE") {
			$sample2trans2counters{$sample}->{$transcript}->[2] += 2-$i;
		    }
		    elsif ($impact eq "LOW") {
			$sample2trans2counters{$sample}->{$transcript}->[3] += 2-$i;
		    }
		    elsif ($impact eq "MODIFIER") {
			# NOOP: we don't want to flag transcripts affected by 2 HET MODIFIER 
			# variants, MODIFIERs of interest should be HV
		    }
		    else {
			die "E $0: unknown impact $impact in inFile $inFile, line:\n$line\n";
		    }
		}
	    }
	}
    }
    close(IN);

    # now print everything we accumulated
    foreach my $sample (keys %sample2lineStarts) {
	foreach my $i (0..$#{$sample2lineStarts{$sample}}) {
	    my $toPrint = $sample2lineStarts{$sample}->[$i];
	    my $transcript = $sample2transcripts{$sample}->[$i];
	    # $numAlleles: number of alleles of severity >= X, X is initially HIGH
	    # and will go down gradually
	    my $numAlleles = $sample2trans2counters{$sample}->{$transcript}->[0];
	    if ($numAlleles >= 2) {
		$toPrint .= "\tHIGH";
	    }
	    else {
		# add the number of MODHIGH alleles
		$numAlleles += $sample2trans2counters{$sample}->{$transcript}->[1];
		if ($numAlleles >= 2) {
		    $toPrint .= "\tMODHIGH";
		}
		else {
		    # add MODERATEs
		    $numAlleles += $sample2trans2counters{$sample}->{$transcript}->[2];
		    if ($numAlleles >= 2) {
			$toPrint .= "\tMODERATE";
		    }
		    else {
			# add LOWs
			$numAlleles += $sample2trans2counters{$sample}->{$transcript}->[3];
			if ($numAlleles >= 2) {
			    $toPrint .= "\tLOW";
			}
			else {
			    $toPrint .= "\tNO";
			}
		    }
		}
	    }
	    $toPrint .= $sample2lineEnds{$sample}->[$i];
	    print { $outFHs{$sample} } $toPrint;
	}
    }
    foreach my $fh (values %outFHs) {
	close($fh);
    }
    $pm->finish;
}
closedir(INDIR);

$pm->wait_all_children;

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";
