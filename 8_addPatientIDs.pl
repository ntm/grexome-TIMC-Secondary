#!/usr/bin/perl

# 12/08/2019
# NTM

# Take 3 arguments: $metadata $inDir $outDir
# $metadata is the patient_summary xlsx file (with sampleID etc columns);
# $inDir must contain cohort TSVs as produced by extractCohorts.pl,
# possibly filtered and reordered with 7_filterAndReorderAll.pl,
# and possibly gzipped;
# $outDir doesn't exist, it will be created and filled with 
# similar TSVs (gzipped if infiles were gzipped), but where every
# $sample identifier in the genoData columns becomes "$sample($patientID)",
# with $patientID taken from patientID column if it's not empty, 
# specimenID otherwise.
# Filenames get ".patientIDs" added before .csv.

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

# key==sample, value is patientID if it exists, specimenID otherwise
my %sample2patient = ();

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
    my ($sampleCol, $specimenCol, $patientCol) = (-1,-1,-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	(defined $cell) || next;
	($cell->value() eq "sampleID") && ($sampleCol = $col);
	($cell->value() eq "specimenID") && ($specimenCol = $col);
	($cell->value() eq "patientID") && ($patientCol = $col);
    }
    ($sampleCol >= 0) ||
	die "E parsing xlsx: no column title is sampleID\n";
    ($specimenCol >= 0) ||
	die "E parsing xlsx: no column title is specimenID\n";
      ($patientCol >= 0) ||
	  die "E parsing xlsx: no column title is patientID\n";
    
    foreach my $row ($rowMin+1..$rowMax) {
	my $sample = $worksheet->get_cell($row, $sampleCol)->value;
	# skip "none" lines
	($sample eq "none") && next;
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
    my ($fileStart,$gz);
    if ($inFile =~ /^(.+)\.csv$/) {
	$fileStart = $1;
    }
    elsif ($inFile =~ /^(.+)\.csv\.gz$/) {
	$fileStart = $1;
	$gz = 1;
    }
    else {
	warn "W: cannot parse filename of inFile $inDir/$inFile, skipping it\n";
    }

    my $inFull = "$inDir/$inFile";
    ($gz) && ($inFull = "gunzip -c $inFull | ");
    open(IN, $inFull) ||
	die "cannot (gunzip-?)open cohort datafile $inDir/$inFile (as $inFull)\n";

    my $outFile = "$outDir/$fileStart.patientIDs.csv";
    my $outFull = " > $outFile";
    ($gz) && ($outFull = " | gzip -c $outFull.gz");
    open (OUT, $outFull) || 
	die "cannot (gzip-?)open $outFile for writing (as $outFull)\n";

    while (my $line = <IN>) {
	chomp($line);
        # add trailing ',' so we know sampleID is always followed by some char
        $line .= ',';
	# chuck norris style: brutal but it works...
	foreach my $sample (keys %sample2patient) {
	    $line =~ s/$sample([\[,\s])/$sample($sample2patient{$sample})$1/g ;
	}
	# remove the trailing ,
	($line =~ s/,$//) || die "E: cannot remove trailing , in:\n$line\n";
	print OUT "$line\n";
    }
    close(IN);
    close(OUT);
}
closedir(INDIR);
