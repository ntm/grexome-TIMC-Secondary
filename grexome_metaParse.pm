# NTM
# 01/04/2021

# Subs for parsing the grexome-TIMC metadata files.
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

# Parse pathologies metadata XLSX file: required columns are "Acronym"
# and "Compatibility groups" (can be in any order but they MUST exist).
#
# Return a hashref:
# - key is a pathology acronym (used as cohort identifier)
# - value is a hashref, with keys==pathos that are compatible with
#   the current patho and values==1
sub parsePathologies {
    (@_ == 1) || die "E: parsePathologies needs one arg";
    my ($pathosFile) = @_;

    # ref to %compatible will be returned, as defined above
    my %compatible = ();

    # we will populate temp hash %compatGroups as we go: key==compatGroup id,
    # value==arrayref of patho acronyms
    my %compatGroups = ();
    
    my $workbook = Spreadsheet::XLSX->new("$pathosFile");
    (defined $workbook) ||
	die "E in parsePathologies: no workbook";
    ($workbook->worksheet_count() == 1) ||
	die "E in parsePathologies: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($pathoCol,$compatCol) = (-1,-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	# if column has no header just ignore it
	(defined $cell) || next;
	if ($cell->value() eq "Acronym") {
	    ($pathoCol = $col);
	}
	elsif ($cell->value() eq "Compatibility groups") {
	    ($compatCol = $col);
	}
    }
    ($pathoCol >= 0) ||
	die "E in parsePathologies: missing required column title: 'Acronym'";
    ($compatCol >= 0) ||
	die "E in parsePathologies: missing required column title: 'Compatibility groups'";

    foreach my $row ($rowMin+1..$rowMax) {
	my $patho = $worksheet->get_cell($row, $pathoCol);
	# skip lines without an acronym
	($patho) || next;
	$patho = $patho->unformatted();
	# require alphanumeric strings
	($patho =~ /^\w+$/) ||
	    die "E in parsePathologies: acronyms must be alphanumeric strings, found \"$patho\" in row $row";
	(defined $compatible{$patho}) && 
	    die "E in parsePathologies: there are 2 lines with same acronym $patho";
	# initialize with anonymous empty hash
	$compatible{$patho} = {};

	# 'Compatibility group' must be a comma-separated list of group identifiers (alphanum strings,
	# eg simply ints)
	my $compats = $worksheet->get_cell($row, $compatCol);
	(defined $compats) || next;
	$compats = $compats->unformatted();
	my @CGs = split(/,/, $compats);
	foreach my $cg (@CGs) {
	    ($cg =~ /^\w+$/) ||
		die "E in parsePathologies: compat group identifiers must be alphanumeric strings, found $cg for $patho";
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









# module loaded ok
1;
