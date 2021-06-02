# NTM
# 01/04/2021

# Subs for parsing the grexome-TIMC metadata files.
# All subs die with an explicit error message if files are not as expected.
# Example metadata files are provided in the Documentation/ subdir.
# This module is shared by the grexome-TIMC-Primary and
# grexome-TIMC-Secondary pipelines.


package grexome_metaParse;

use strict;
use warnings;
use Spreadsheet::XLSX;
use Exporter;
our @ISA = ('Exporter');
our @EXPORT_OK = qw(parsePathologies);


#################################################################

# Parse pathologies metadata XLSX file: required columns are "pathologyID"
# and "compatibility groups" (can be in any order but they MUST exist).
#
# Return a hashref:
# - key is a pathologyID (used as cohort identifier)
# - value is a hashref, with keys==pathoIDs that are compatible with
#   the current patho and values==1
sub parsePathologies {
    my $subName = (caller(0))[3];
    (@_ == 1) || die "E: $subName needs one arg";
    my ($pathosFile) = @_;

    (-f $pathosFile) || die "E in $subName: provided file $pathosFile doesn't exist";

    # ref to %compatible will be returned, as defined above
    my %compatible = ();

    # we will populate temp hash %compatGroups as we go: key==compatGroup id,
    # value==arrayref of patho acronyms
    my %compatGroups = ();
    
    my $workbook = Spreadsheet::XLSX->new("$pathosFile");
    (defined $workbook) ||
	die "E in $subName: no workbook";
    ($workbook->worksheet_count() == 1) ||
	die "E in $subName: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($pathoCol,$compatCol) = (-1,-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	# if column has no header just ignore it
	(defined $cell) || next;
	my $val = $cell->unformatted();
	if ($val eq "pathologyID") { $pathoCol = $col; }
	elsif ($val eq "compatibility groups") { $compatCol = $col; }
    }
    ($pathoCol >= 0) ||
	die "E in $subName: missing required column title: 'pathologyID'";
    ($compatCol >= 0) ||
	die "E in $subName: missing required column title: 'compatibility groups'";

    foreach my $row ($rowMin+1..$rowMax) {
	my $patho = $worksheet->get_cell($row, $pathoCol);
	# skip lines without a pathoID
	($patho) || next;
	$patho = $patho->unformatted();
	# require alphanumeric strings
	($patho =~ /^\w+$/) ||
	    die "E in $subName: pathologyIDs must be alphanumeric strings, found \"$patho\" in row $row";
	(defined $compatible{$patho}) && 
	    die "E in $subName: there are 2 lines with same pathologyID $patho";
	# initialize with anonymous empty hash
	$compatible{$patho} = {};

	# 'Compatibility group' must be a comma-separated list of group identifiers (alphanum strings)
	my $compats = $worksheet->get_cell($row, $compatCol);
	(defined $compats) || next;
	$compats = $compats->unformatted();
	my @CGs = split(/,/, $compats);
	foreach my $cg (@CGs) {
	    ($cg =~ /^\w+$/) ||
		die "E in $subName: compat group identifiers must be alphanumeric strings, found $cg for $patho";
	    (defined $compatGroups{$cg}) || ($compatGroups{$cg} = []);
	    push(@{$compatGroups{$cg}}, $patho);
	}
    }

    # now populate %compatible from %compatGroups
    foreach my $cgs (values(%compatGroups)) {
	foreach my $patho (@$cgs) {
	    foreach my $compatPatho (@$cgs) {
		($compatPatho eq $patho) && next;
		$compatible{$patho}->{$compatPatho} = 1;
	    }
	}
    }
    
    return(\%compatible);
}

#################################################################

# Parse samples metadata XLSX file: required columns are "sampleID",
# "specimenID", "patientID", "pathologyID", and "Causal gene" (they
# can be in any order but they MUST exist).
# The optional second argument, if povided and non-empty, is the
# pathologies metadata XLSX file; it is used to make sure every
# pathologyID in samples.xlsx is defined in pathologies.xlsx (ie
# sanity-check for typoes).
#
# Return a list of 4 hashrefs (caller can use whichever it needs):
# - sample2patho: key is sampleID, value is pathologyID
# - sample2specimen: key is sampleID, value is specimenID
# - sample2patient: key is sampleID, value is patientID if it exists, specimenID otherwise
# - sample2causal: key is sampleID, value is casual gene if it exists (undef otherwise)
sub parseSamples {
    my $subName = (caller(0))[3];
    (@_ == 1) || (@_ == 2) || die "E: $subName needs one or two args";
    my ($samplesFile, $pathosFile) = @_;

    (-f $samplesFile) ||
	die "E in $subName: provided samples file $samplesFile doesn't exist";
    (! $pathosFile) || (-f $pathosFile) ||
	die "E in $subName: optional pathologies file $pathosFile was provided but doesn't exist";

    # sample2patho: key is sampleID, value is pathologyID
    my %sample2patho;
    # sample2specimen: key is sampleID, value is specimenID
    my %sample2specimen;
    # sample2patient: key is sampleID, value is patientID if it exists, specimenID otherwise
    my %sample2patient;
    # sample2causal: key is sampleID, value is casual gene if it exists (undef otherwise)
    my %sample2causal;

    ################
    # parse $pathosFile immediately if it was provided
    my $pathologiesR;
    if ($pathosFile) {
	$pathologiesR = &parsePathologies($pathosFile);
    }

    ################
    # parse headers
    my $workbook = Spreadsheet::XLSX->new("$samplesFile");
    (defined $workbook) ||
	die "E in $subName: no workbook";
    ($workbook->worksheet_count() == 1) ||
	die "E in $subName: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($sampleCol, $pathoCol, $specimenCol, $patientCol, $causalCol) = (-1,-1,-1,-1,-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	# if column has no header just ignore it
	(defined $cell) || next;
	my $val = $cell->unformatted();
	if ($val eq "sampleID") { $sampleCol = $col; }
	elsif ($val eq "pathologyID") { $pathoCol = $col; }
	elsif ($val eq "specimenID") { $specimenCol = $col; }
	elsif ($val eq "patientID") { $patientCol = $col; }
	elsif ($val eq "Causal gene") { $causalCol = $col; }
    }
    ($sampleCol >= 0) ||
	die "E in $subName: missing required column title: 'sampleID'";
    ($pathoCol >= 0) ||
	die "E in $subName: missing required column title: 'pathologyID'";
    ($specimenCol >= 0) ||
	die "E in $subName: missing required column title: 'specimenID'";
    ($patientCol >= 0) ||
	die "E in $subName: missing required column title: 'patientID'";
    ($causalCol >= 0) ||
	die "E in $subName: missing required column title: 'Causal gene'";

    ################
    # parse data rows
    foreach my $row ($rowMin+1..$rowMax) {
	my $sample = $worksheet->get_cell($row, $sampleCol);
	($sample) ||
	    die "E in $subName, row $row: every row MUST have a sampleID (use 0 for obsolete samples)";
	$sample = $sample->unformatted();
	# skip "0" lines == obsolete samples
	($sample eq "0") && next;
	(defined $sample2patho{$sample}) && 
	    die "E in $subName: found 2 lines with same sampleID $sample";
	
	################ sample2patho
	my $patho = $worksheet->get_cell($row, $pathoCol);
	($patho) || die "E in $subName, row $row: every row MUST have a pathologyID";
	$patho = $patho->unformatted();
	if ($pathosFile) {
	    (defined $pathologiesR->{$patho}) ||
		die "E in $subName, row $row: pathologyID $patho is not defined in the provided $pathosFile";
	}
	else {
	    # require alphanumeric strings
	    ($patho =~ /^\w+$/) ||
		die "E in $subName, row $row: pathologyIDs must be alphanumeric strings, found \"$patho\"";
	}
	$sample2patho{$sample} = $patho;

	################ sample2specimen
	my $specimen = $worksheet->get_cell($row, $specimenCol);
	($specimen) || die "E in $subName, row $row: every row MUST have a specimenID";
	$specimen = $specimen->unformatted();
	# require alphanumeric strings (really don't want spaces or shell metachars)
	($specimen =~ /^\w+$/) ||
	    die "E in $subName, row $row: specimenIDs must be alphanumeric strings, found \"$specimen\"";
	$sample2specimen{$sample} = $specimen;

	################ sample2patient
	my $patient = $worksheet->get_cell($row, $patientCol);
	# patientID value is optional, we use specimenID if it's missing
	if ($patient) {
	    $patient = $patient->unformatted();
	    # require alphanum (or underscores as always) or dashes, this could probably
	    # be relaxed
	    ($patient =~ /^[\w-]+$/) ||
		die "E in $subName, row $row: patientIDs must be alphanumeric (dashes allowed), found \"$patient\"";
	}
	else {
	    $patient = $specimen;
	}
	$sample2patient{$sample} = $patient;

	################ sample2causal
	my $causal = $worksheet->get_cell($row, $causalCol);
	if ($causal) {
	    $causal = $causal->unformatted();
	    # these are HUGO gene names -> must be alphanum+dashes
	    ($causal =~ /^[\w-]+$/) ||
		die "E in $subName, row $row: causalGene must be alphanumeric (dashes allowed), found \"$causal\"";
	    $sample2causal{$sample} = $causal;
	}
    }

    ################
    return(\%sample2patho, \%sample2specimen, \%sample2patient, \%sample2causal);
}




# module loaded ok
1;
