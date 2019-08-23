#!/usr/bin/perl

# 25/03/2018
# NTM

# Take 3 arguments: $metadata $inDir $outDir
# $metadata is the patient_summary xlsx file (with grexomeID etc columns);
# $inDir must contain cohort TSVs as produced by extractCohorts.pl,
# possibly filtered with finalFilters.pl, and possibly gzipped;
# $outDir doesn't exist, it will be created and filled with one TSV
# per sample. 
# Filenames will include patientID/specimenID.
# For a sample, we only print lines from its cohort file and where 
# it has an HV or HET genotype: this genotype is printed in a new 
# column GENOTYPE, inserted right after KNOWN_CANDIDATE_GENE.
# Also the HV, NEGCTRL_HV, HET etc... columns are not printed.
#
# NOTE: the fact that HV et al are the last 6 columns (allowing for 
# one more column with filter logs) is hard-coded.
# If this isn't true the script dies.

use strict;
use warnings;
use Spreadsheet::XLSX;


(@ARGV == 3) || die "needs 3 args: a metadata XLSX, an inDir and a non-existant outDir\n";
my ($metadata, $inDir, $outDir) = @ARGV;
(-d $inDir) ||
    die "inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "cannot opendir inDir $inDir\n";
(-e $outDir) && 
    die "found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "cannot mkdir outDir $outDir\n";

#########################################################
# parse metadata file

# key==cohort name, value is an arrayref of all grexomes from this cohort
my %cohort2grexomes = ();

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
	my $grexome = $worksheet->get_cell($row, $grexCol)->value;
	# skip "none" lines
	($grexome eq "none") && next;
	my $cohort = $worksheet->get_cell($row, $cohortCol)->value;
	(defined $cohort2grexomes{$cohort}) || ($cohort2grexomes{$cohort} = []);
	push(@{$cohort2grexomes{$cohort}}, $grexome);
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

    # KNOWN_CANDIDATE_GENE column 
    my $knownCandidateCol = -1;

    my $inFull = "$inDir/$inFile";
    ($gz) && ($inFull = "gunzip -c $inFull | ");
    open(IN, $inFull) ||
	die "cannot (gunzip-?)open cohort datafile $inDir/$inFile (as $inFull)\n";
    my $header = <IN>;
    chomp($header);
    my @header = split(/\t/,$header);
    foreach my $i (0..$#header) {
	if ($header[$i] eq "KNOWN_CANDIDATE_GENE") {
	    $knownCandidateCol = $i;
	    $header[$i] .= "\tGENOTYPE";
	    last;
	}
    }
    ($knownCandidateCol >= 0) || 
	die "E: couldn't find KNOWN_CANDIDATE_GENE in header of infile $inFile\n";
    $header = join("\t",@header);
    ($header =~ s/\tHV\tNEGCTRL_HV\tHET\tNEGCTRL_HET\tOTHER\tNEGCTRL_OTHER$//) ||
	($header =~ s/\tHV\tNEGCTRL_HV\tHET\tNEGCTRL_HET\tOTHER\tNEGCTRL_OTHER(\tmax_ctrl_hv=[^\t]+)$/$1/) ||
	die "cannot remove HV HET OTHER from header of inFile $inFile\n$header\n";

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
	print $FH "$header\n";
	$outFHs{$grexome} = $FH ;
    }

    # now read the data
    while (my $line = <IN>) {
	chomp($line);
	my @fields = split(/\t/, $line, -1) ;
	my $toPrintStart = join("\t",@fields[0..$knownCandidateCol])."\t";
	my $toPrintEnd = join("\t",@fields[($knownCandidateCol+1)..($#fields-6)])."\n";

	foreach my $i ($#fields-5,$#fields-3) {
	    if ($fields[$i]) {
		($fields[$i] =~ /^([^~]+)~([^~\|]+)$/) || 
		    die "cannot parse HV/HET data $fields[$i] from infile $inFile\n";
		my ($geno,$samples) = ($1,$2);
		# actually, just use HV or HET for geno, the actual allele is in ALLELE_NUM
		# we will still add [DP:AF] after HV/HET
		($i == $#fields-5) && ($geno = "HV");
		($i == $#fields-3) && ($geno = "HET");
		foreach my $sample (split(/,/,$samples)) {
		    # grab grexome and [DP:AF], we know it must be there in HV and HET columns
		    # (allowing AF > 1 for Strelka bug)
		    ($sample =~ /^(grexome\d\d\d\d)(\[\d+:\d+\.\d\d\])$/) ||
			die  "E: inFile $inFile has a genotype call for a sample I can't parse: $sample\n";
		    my ($grexome,$dpaf) = ($1,$2);
		    print { $outFHs{$grexome} } "$toPrintStart$geno$dpaf\t$toPrintEnd" ;
		}
	    }
	}
    }
    close(IN);
    foreach my $fh (values %outFHs) {
	close($fh);
    }
}
closedir(INDIR);
