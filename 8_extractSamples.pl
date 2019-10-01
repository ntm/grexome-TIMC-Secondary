#!/usr/bin/perl

# 25/03/2018
# NTM

# Take 4 arguments: $metadata $inDir $covDir $outDir
# $metadata is the patient_summary xlsx file (with grexomeID etc columns);
# $inDir must contain cohort TSVs as produced by extractCohorts.pl,
# possibly filtered and reordered with 7_filterAndReorderAll.pl, and 
# possibly gzipped;
# $covDir is a subdir containing per-grexome coverage files, as produced
# by 0_coverage.pl;
# $outDir doesn't exist, it will be created and filled with one TSV
# per sample. 
# Filenames will include patientID/specimenID.
# The global coverage data (ALL_CANDIDATES and ALL_SAMPLED) for each
# grexome is grabbed from $covDir and added at the end of the header line.
# For a sample, we only print lines from its cohort file and where 
# it has an HV or HET genotype: this genotype is printed in new columns
# GENOTYPE and DP:AF, inserted right after KNOWN_CANDIDATE_GENE.
# For samples with a causal gene, we also grab its genotype in the NEGCTRL_HV
# and NEGCTRL_HET columns (so the sample file doesn't have only variants
# impacting its causal gene).
# Between GENOTYPE and DP:AF we insert new columns NB_$geno_$impact_ThisSample_ThisTranscript
# with $geno == HV or HET and $impact == HIGH or MODER, counting the total
# number of $geno-$impact variants (passing all filters and) affecting
# this transcript in this grexome; and then NB_ALLVARIANTS_ThisSample_ThisTranscript
# with all variants affecting this transcript in this grexome.
# Also the HV, NEGCTRL_HV, HET etc... columns are not printed.

use strict;
use warnings;
use Spreadsheet::XLSX;


(@ARGV == 4) || 
    die "needs 4 args: a metadata XLSX, an inDir, a covDir and a non-existant outDir\n";
my ($metadata, $inDir, $covDir, $outDir) = @ARGV;
(-d $inDir) ||
    die "inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "cannot opendir inDir $inDir\n";
(-d $covDir) ||
    die "covDir $covDir doesn't exist or isn't a directory\n";
(-e $outDir) && 
    die "found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "cannot mkdir outDir $outDir\n";

#########################################################
# parse metadata file

# key==cohort name, value is an arrayref of all grexomes from this cohort
my %cohort2grexomes = ();
# also the other way around: key==grexome, value==its cohort
my %grexome2cohort = ();

# key==grexome, value is patientID if it exists, specimenID otherwise
my %grexome2patient = ();

(-f $metadata) ||
    die "E: the supplied metadata file doesn't exist\n";
{
    my $workbook = Spreadsheet::XLSX->new("$metadata");
    (defined $workbook) ||
	die "E when parsing xlsx\n";
    ($workbook->worksheet_count() == 1) ||
	die "E parsing xlsx: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($grexCol, $cohortCol, $specimenCol, $patientCol) = (-1,-1,-1,-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	(defined $cell) || next;
	($cell->value() eq "grexomeID") && ($grexCol = $col);
	($cell->value() eq "pathology") && ($cohortCol = $col);
	($cell->value() eq "specimenID") && ($specimenCol = $col);
	($cell->value() eq "patientID") && ($patientCol = $col);
    }
    ($grexCol >= 0) ||
	die "E parsing xlsx: no column title is grexomeID\n";
    ($cohortCol >= 0) ||
	die "E parsing xlsx: no col title is pathology\n";
    ($specimenCol >= 0) ||
	die "E parsing xlsx: no column title is specimenID\n";
    ($patientCol >= 0) ||
	die "E parsing xlsx: no column title is patientID\n";
    
    foreach my $row ($rowMin+1..$rowMax) {
	(defined $worksheet->get_cell($row, $grexCol)) ||
	    die "E: cell undefined for row $row, col $grexCol\n";
	my $grexome = $worksheet->get_cell($row, $grexCol)->value;
	# skip "none" lines
	($grexome eq "none") && next;
	my $cohort = $worksheet->get_cell($row, $cohortCol)->value;
	(defined $cohort2grexomes{$cohort}) || ($cohort2grexomes{$cohort} = []);
	push(@{$cohort2grexomes{$cohort}}, $grexome);
	$grexome2cohort{$grexome} = $cohort;
	my $patient = $worksheet->get_cell($row, $specimenCol)->unformatted();
	if ($worksheet->get_cell($row, $patientCol)) {
	    my $tmp = $worksheet->get_cell($row, $patientCol)->unformatted();
	    $tmp =~ s/^\s+//;
	    $tmp =~ s/\s+$//;
	    ($tmp) && ($patient = $tmp);
	}
	$grexome2patient{$grexome} = $patient;
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
	warn "W: cannot parse filename of inFile $inDir/$inFile, skipping it\n";
    }

    my $inFull = "$inDir/$inFile";
    ($gz) && ($inFull = "gunzip -c $inFull | ");
    open(IN, $inFull) ||
	die "cannot (gunzip-?)open cohort datafile $inDir/$inFile (as $inFull)\n";
    my $header = <IN>;
    chomp($header);
    my @header = split(/\t/,$header);

    # need HV HET and other Geno columns (will all be removed)
    my ($hvCol,$hetCol) = (-1,-1);
    my ($hvNegCol,$hetNegCol) = (-1,-1);
    # $colsToRemove[$i] == 1 if column $i in infile must be removed
    my @colsToRemove;
    # reverse so we can splice out elements
    foreach my $i (reverse(0..$#header)) {
	if ($header[$i] eq "HV") {
	    $hvCol = $i;
	    $colsToRemove[$i] = 1;
	    splice(@header,$i,1);
	}
	elsif ($header[$i] eq "HET") {
	    $hetCol = $i;
	    $colsToRemove[$i] = 1;
	    splice(@header,$i,1);
	}
	if ($header[$i] eq "NEGCTRL_HV") {
	    $hvNegCol = $i;
	    $colsToRemove[$i] = 1;
	    splice(@header,$i,1);
	}
	elsif ($header[$i] eq "NEGCTRL_HET") {
	    $hetNegCol = $i;
	    $colsToRemove[$i] = 1;
	    splice(@header,$i,1);
	}
	elsif (grep(/$header[$i]/, ("OTHER","NEGCTRL_OTHER"))) {
	    $colsToRemove[$i] = 1;
	    splice(@header,$i,1);
	}
    }
    ($hvCol >= 0) || ($hvNegCol >= 0) || 
	die "E: couldn't find HV or NEGCTRL_HV in header of infile $inFile\n";
    ($hetCol >= 0) || ($hetNegCol >= 0) || 
	die "E: couldn't find HET or NEGCTRL_HET in header of infile $inFile\n";

    # KNOWN_CANDIDATE_GENE, Feature and IMPACT column indexes after
    # removing @colsToRemove columns
    my ($knownCandidateCol,$featureCol, $impactCol) = (-1,-1,-1);
 
    foreach my $i (0..$#header) {
	if ($header[$i] eq "KNOWN_CANDIDATE_GENE") {
	    $knownCandidateCol = $i;
	    $header[$i] .= "\tGENOTYPE";
	    foreach my $g ("HV", "HET") {
		foreach my $im ("HIGH", "MODER") {
		    $header[$i] .= "\tNB_$g"."_$im"."_ThisSample_ThisTranscript";
		}
	    }
	    $header[$i] .= "\tNB_ALLVARIANTS_ThisSample_ThisTranscript";
	    $header[$i] .= "\tDP:AF";
	}
	elsif ($header[$i] eq "Feature") {
	    $featureCol = $i;
	}
	elsif ($header[$i] eq "IMPACT") {
	    $impactCol = $i;
	}
    }
    ($knownCandidateCol >= 0) || 
	die "E: couldn't find KNOWN_CANDIDATE_GENE in header of infile $inFile\n";
    ($featureCol >= 0) || 
	die "E: couldn't find Feature in header of infile $inFile\n";
    ($impactCol >= 0) || 
	die "E: couldn't find IMPACT in header of infile $inFile\n";

    $header = join("\t",@header);

    # hash of filehandles open for writing, one for each grexome
    # from this cohort
    # will be gzipped if infiles were
    my %outFHs;

    ($cohort2grexomes{$cohort}) || 
	die "cohort $cohort parsed from filename of infile $inFile is not in $metadata\n";
    foreach my $grexome (@{$cohort2grexomes{$cohort}}) {
	my $patient = $grexome2patient{$grexome};
	my $outFile = "$outDir/$cohort.$grexome.$patient.$fileEnd";
	($gz) && ($outFile .= ".gz");
	my $outFull = " > $outFile";
	($gz) && ($outFull = " | gzip -c $outFull");
	open (my $FH, $outFull) || die "cannot (gzip-?)open $outFile for writing (as $outFull)\n";

	# grab global coverage data for $grexome
	my $covFile = "$covDir/coverage_$grexome.csv";
	(-f $covFile) || 
	    die "E: trying to grab coverage data for $grexome but covFile doesn't exist: $covFile\n";
	# global coverage data is in last 2 lines
	open(COV, "tail -n 2 $covFile |") ||
	    die "cannot tail-grab cov data from covFile with: tail -n 2 $covFile\n";
	# ALL_CANDIDATES
	my $covLine = <COV>;
	chomp $covLine;
	my @covFields = split(/\t/,$covLine);
	(@covFields == 8) || die "E: expecting 8 fields from candidates coverage line $covLine\n";
	my $headerCov = "\tCoverage_Candidates_50x=$covFields[5] Coverage_Candidates_20x=$covFields[6] Coverage_Candidates_10x=$covFields[7]";
	# ALL_GENES
	$covLine = <COV>;
	chomp $covLine;
	@covFields = split(/\t/,$covLine);
	(@covFields == 8) || die "E: expecting 8 fields from sampled coverage line $covLine\n";
	$headerCov .= "   Coverage_AllGenes_50x=$covFields[5] Coverage_AllGenes_20x=$covFields[6] Coverage_AllGenes_10x=$covFields[7]";
	print $FH "$header$headerCov\n";
	$outFHs{$grexome} = $FH ;
	close(COV);
    }

    # now read the data
    # in order to print NB_* we need:
    # key == grexome, values are arrayrefs of the beginnings and ends of lines
    # (respectively) that must be printed for this grexome
    my %grex2lineStarts;
    my %grex2lineEnds;
    # key == grexome, value is an arrayref, for each line to print it holds
    # the transcript that this line deals with
    my %grex2transcripts;
    # key == grexome, value is a hashref whose keys are transcripts and values
    # are arrayrefs with 5 ints: numbers of HV_HIGH, HV_MODER, HET_HIGH, HET_MODER and ALL
    # found for this transcript in this grexome
    my %grex2trans2counters;

    while (my $line = <IN>) {
	chomp($line);
	my @fields = split(/\t/, $line, -1) ;
	# grab needed GENO data: HV at index 0, HET at index 1, 
	# NEGCTRL_HV at index 2, NEGCTRL_HET at index 3
	my @genoData;
	($fields[$hvCol]) && ($genoData[0] = $fields[$hvCol]);
	($fields[$hetCol]) && ($genoData[1] = $fields[$hetCol]);
	($fields[$hvNegCol]) && ($genoData[2] = $fields[$hvNegCol]);
	($fields[$hetNegCol]) && ($genoData[3] = $fields[$hetNegCol]);
	# splice all GENO data out
	foreach my $i (reverse(0..$#colsToRemove)) {
	    ($colsToRemove[$i]) && splice(@fields,$i,1);
	}

	my $toPrintStart = join("\t",@fields[0..$knownCandidateCol])."\t";
	my $toPrintEnd = join("\t",@fields[($knownCandidateCol+1)..$#fields])."\n";
	my $transcript = $fields[$featureCol];
	my $impact = $fields[$impactCol];

	foreach my $i (0..3) {
	    if ($genoData[$i]) {
		($genoData[$i] =~ /^([^~]+)~([^~\|]+)$/) || 
		    die "cannot parse HV/HET data $genoData[$i] from infile $inFile\n";
		my ($geno,$samples) = ($1,$2);
		# actually, just use HV or HET for geno, the actual allele is in ALLELE_NUM
		# we will still add DP:AF after HV/HET
		($i % 2 == 0) && ($geno = "HV");
		($i % 2 == 1) && ($geno = "HET");
		foreach my $sample (split(/,/,$samples)) {
		    # grab grexome and [DP:AF], we know it must be there in HV and HET columns
		    # (allowing AF > 1 for Strelka bug)
		    ($sample =~ /^(grexome\d\d\d\d)\[(\d+:\d+\.\d\d)\]$/) ||
			die  "E: inFile $inFile has a genotype call for a sample I can't parse: $sample\n";
		    my ($grexome,$dpaf) = ($1,$2);
		    # ignore grexomes from other cohorts in NEGCTRL* columns
		    if (($i >= 2) && ($grexome2cohort{$grexome} ne $cohort)) {
			next;
		    }
		    # initialize everything for this grexome if needed
		    ($grex2lineStarts{$grexome}) || ($grex2lineStarts{$grexome} = []);
		    ($grex2lineEnds{$grexome}) || ($grex2lineEnds{$grexome} = []);
		    ($grex2transcripts{$grexome}) || ($grex2transcripts{$grexome} = []);
		    ($grex2trans2counters{$grexome}) || ($grex2trans2counters{$grexome} = {});
		    ($grex2trans2counters{$grexome}->{$transcript}) || ($grex2trans2counters{$grexome}->{$transcript} = [0,0,0,0]);

		    # now fill our data structures
		    push(@{$grex2lineStarts{$grexome}}, "$toPrintStart$geno");
		    push(@{$grex2lineEnds{$grexome}}, "\t$dpaf\t$toPrintEnd");
		    push(@{$grex2transcripts{$grexome}}, $transcript);
		    my $indexToIncr = 0;
		    ($geno eq "HET") && ($indexToIncr += 2);
		    if ($impact eq "HIGH") {
			$grex2trans2counters{$grexome}->{$transcript}->[$indexToIncr]++;
		    }
		    elsif ($impact eq "MODERATE") {
			$grex2trans2counters{$grexome}->{$transcript}->[$indexToIncr+1]++;
		    }
		    # in any case increment the ALL counter
		    $grex2trans2counters{$grexome}->{$transcript}->[4]++;
		}
	    }
	}
    }
    close(IN);

    # now print everything we accumulated
    foreach my $grexome (keys %grex2lineStarts) {
	foreach my $i (0..$#{$grex2lineStarts{$grexome}}) {
	    my $toPrint = $grex2lineStarts{$grexome}->[$i];
	    my $transcript = $grex2transcripts{$grexome}->[$i];
	    $toPrint .= "\t".join("\t", @{$grex2trans2counters{$grexome}->{$transcript}});
	    $toPrint .= $grex2lineEnds{$grexome}->[$i];

	    print { $outFHs{$grexome} } $toPrint;
	}
    }
    foreach my $fh (values %outFHs) {
	close($fh);
    }
}
closedir(INDIR);
