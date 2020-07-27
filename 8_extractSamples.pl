#!/usr/bin/perl

# 25/03/2018
# NTM

# Take 4 arguments: $metadata $inDir $covDir $outDir
# $metadata is the patient_summary xlsx file (with sampleID etc columns);
# $inDir must contain cohort TSVs as produced by extractCohorts.pl,
# possibly filtered and reordered with 7_filterAndReorderAll.pl, and 
# possibly gzipped (but not with PatientIDs);
# $covDir is a subdir containing per-sample coverage files, as produced
# by 0_coverage.pl;
# $outDir doesn't exist, it will be created and filled with one TSV
# per sample. 
# Filenames will include patientID/specimenID.
# The global coverage data (ALL_CANDIDATES and ALL_SAMPLED) for each
# sample is grabbed from $covDir and added at the end of the header line.
# For a sample, we only print lines from its cohort file and where 
# it has an HV or HET genotype: this genotype is printed in new columns
# GENOTYPE and DP:AF, inserted right after KNOWN_CANDIDATE_GENE.
# Immediately after DP:AF we insert a new column BOTH_ALLELES, value is one of:
#   HIGH -> patient has >=2 HET or at >=1 HV HIGH variants;
#   MODHIGH -> patient has >=2 HET or >=1 HV variants of impact HIGH or MODHIGH,
#     but isn't in HIGH category;
#   MODERATE -> same for impact >= MODERATE;
#   LOW -> same for impact >= LOW;
#   NO -> patient has at most one allele >= LOW.
# Also the genoData columns ($cohort_HV, NEGCTRL_HV, etc...) are not printed.

use strict;
use warnings;
use Spreadsheet::XLSX;


(@ARGV == 4) || 
    die "E: $0 - needs 4 args: a metadata XLSX, an inDir, a covDir and a non-existant outDir\n";
my ($metadata, $inDir, $covDir, $outDir) = @ARGV;
(-d $inDir) ||
    die "E: $0 - inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "E: $0 - cannot opendir inDir $inDir\n";
(-d $covDir) ||
    die "E: $0 - covDir $covDir doesn't exist or isn't a directory\n";
(-e $outDir) && 
    die "E: $0 - found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "E: $0 - cannot mkdir outDir $outDir\n";

#########################################################
# parse metadata file

# key==cohort name, value is an arrayref of all samples from this cohort
my %cohort2samples = ();

# key==sample, value is patientID if it exists, specimenID otherwise
my %sample2patient = ();

(-f $metadata) ||
    die "E: $0 - the supplied metadata file doesn't exist\n";
{
    my $workbook = Spreadsheet::XLSX->new("$metadata");
    (defined $workbook) ||
	die "E: $0 - E when parsing xlsx\n";
    ($workbook->worksheet_count() == 1) ||
	die "E: $0 - parsing xlsx: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($sampleCol, $cohortCol, $specimenCol, $patientCol) = (-1,-1,-1,-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	(defined $cell) || next;
	($cell->value() eq "sampleID") && ($sampleCol = $col);
	($cell->value() eq "pathology") && ($cohortCol = $col);
	($cell->value() eq "specimenID") && ($specimenCol = $col);
	($cell->value() eq "patientID") && ($patientCol = $col);
    }
    ($sampleCol >= 0) ||
	die "E: $0 - parsing xlsx: no column title is sampleID\n";
    ($cohortCol >= 0) ||
	die "E: $0 - parsing xlsx: no col title is pathology\n";
    ($specimenCol >= 0) ||
	die "E: $0 - parsing xlsx: no column title is specimenID\n";
    ($patientCol >= 0) ||
	die "E: $0 - parsing xlsx: no column title is patientID\n";
    
    foreach my $row ($rowMin+1..$rowMax) {
	(defined $worksheet->get_cell($row, $sampleCol)) ||
	    die "E: $0 - cell undefined for row $row, col $sampleCol\n";
	my $sample = $worksheet->get_cell($row, $sampleCol)->value;
	# skip "none" lines
	($sample eq "none") && next;
	my $cohort = $worksheet->get_cell($row, $cohortCol)->value;
	(defined $cohort2samples{$cohort}) || ($cohort2samples{$cohort} = []);
	push(@{$cohort2samples{$cohort}}, $sample);
	my $patient = $worksheet->get_cell($row, $specimenCol)->unformatted();
	if ($worksheet->get_cell($row, $patientCol)) {
	    my $tmp = $worksheet->get_cell($row, $patientCol)->unformatted();
	    $tmp =~ s/^\s+//;
	    $tmp =~ s/\s+$//;
	    ($tmp) && ($patient = $tmp);
	}
	$sample2patient{$sample} = $patient;
    }
}

#########################################################
# read infiles

while (my $inFile = readdir(INDIR)) {
    ($inFile =~ /^\./) && next;
    my ($cohort,$fileEnd,$gz);
    if ($inFile =~ (/^([^\.]+)\.(.*csv)$/)) {
	# $fileEnd allows for .filtered , .pick etc...
	($cohort,$fileEnd) = ($1,$2);
    }
    elsif ($inFile =~ (/^([^\.]+)\.(.*csv)\.gz$/)) {
	($cohort,$fileEnd) = ($1,$2);
	$gz = 1;
    }
    else {
	warn "W $0: cannot parse filename of inFile $inDir/$inFile, skipping it\n";
    }

    my $inFull = "$inDir/$inFile";
    ($gz) && ($inFull = "gunzip -c $inFull | ");
    open(IN, $inFull) ||
	die "E: $0 - cannot (gunzip-?)open cohort datafile $inDir/$inFile (as $inFull)\n";
    my $header = <IN>;
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
	die "E: $0 - couldn't find {$cohort}_HV or {$cohort}_OTHERCAUSE_HV in header of infile $inFile\n";
    (($hetCol >= 0) && ($hetColOC >= 0)) || 
	die "E: $0 - couldn't find HET or OC_HET in header of infile $inFile\n";

    # KNOWN_CANDIDATE_GENE, Feature and IMPACT column indexes after
    # removing @colsToRemove columns
    my ($knownCandidateCol,$featureCol, $impactCol) = (-1,-1,-1);
 
    foreach my $i (0..$#header) {
	if ($header[$i] eq "KNOWN_CANDIDATE_GENE") {
	    $knownCandidateCol = $i;
	    $header[$i] .= "\tGENOTYPE";
	    $header[$i] .= "\tDP:AF";
	    $header[$i] .= "\tBOTH_ALLELES";
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
	die "E: $0 - cohort $cohort parsed from filename of infile $inFile is not in $metadata\n";
    foreach my $sample (@{$cohort2samples{$cohort}}) {
	my $patient = $sample2patient{$sample};
	my $outFile = "$outDir/$cohort.$sample.$patient.$fileEnd";
	($gz) && ($outFile .= ".gz");
	my $outFull = " > $outFile";
	($gz) && ($outFull = " | gzip -c $outFull");
	open (my $FH, $outFull) || die "E: $0 - cannot (gzip-?)open $outFile for writing (as $outFull)\n";

	# grab global coverage data for $sample
	my $covFile = "$covDir/coverage_$sample.csv";
	(-f $covFile) || 
	    die "E $0: trying to grab coverage data for $sample but covFile doesn't exist: $covFile\n";
	# global coverage data is in last 2 lines
	open(COV, "tail -n 2 $covFile |") ||
	    die "E: $0 - cannot tail-grab cov data from covFile with: tail -n 2 $covFile\n";
	# ALL_CANDIDATES
	my $covLine = <COV>;
	chomp $covLine;
	my @covFields = split(/\t/,$covLine);
	(@covFields == 8) || die "E $0: expecting 8 fields from candidates coverage line $covLine\n";
	my $headerCov = "\tCoverage_Candidates_50x=$covFields[5] Coverage_Candidates_20x=$covFields[6] Coverage_Candidates_10x=$covFields[7]";
	# ALL_GENES
	$covLine = <COV>;
	chomp $covLine;
	@covFields = split(/\t/,$covLine);
	(@covFields == 8) || die "E $0: expecting 8 fields from sampled coverage line $covLine\n";
	$headerCov .= "   Coverage_AllGenes_50x=$covFields[5] Coverage_AllGenes_20x=$covFields[6] Coverage_AllGenes_10x=$covFields[7]";
	print $FH "$header$headerCov\n";
	$outFHs{$sample} = $FH ;
	close(COV);
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
		    # grab sample and [DP:AF], we know it must be there
		    # (allowing AF > 1 for Strelka bug)
		    ($sampleData =~ /^([^\[\s]+)\[(\d+:\d+\.\d\d)\]$/) ||
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
}
closedir(INDIR);
