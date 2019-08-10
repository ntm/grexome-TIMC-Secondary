#!/usr/bin/perl

# 24/03/2018
# NTM 

# Takes a single arg: $metadata, an xlsx file (with path) with
# "grexomeID" in some col and "pathology" in another (eg 
# patient_summary*.xlsx). 
# Read on stdin a VCF file with GENOS data columns, as produced
# by sampleData2genotypes.pl.
# Print to stdout a VCF file similar to the infile, but adding INFO
# fields: 
# COUNT_$cohort_HV, COUNT_$cohort_HET, COUNT_$cohort_OTHER, COUNT_$cohort_HR
# values are geno:count,geno:count.. eg COUNT_ovo_HV=1/1:12,2/2:3

use strict;
use warnings;
use Spreadsheet::XLSX;

(@ARGV == 1) || die "$0 needs one arg: the metadata xlsx\n";
my $metadata = shift(@ARGV);

# key==sample id, value is the $cohort this sample belongs to
my %sample2cohort = ();
# cohort names
my @cohorts = ();

# the genotype categories (don't change please, hard-coded in other scripts)
my @genoCategories = ("HV","HET","OTHER","HR");

#########################################################
# parse metadata file
(-f $metadata) ||
    die "E: the supplied metadata file doesn't exist\n";
{
    # for cohort names we use a temp hash to avoid redundancy
    my %cohorts;
    my $workbook = Spreadsheet::XLSX->new("$metadata");
    (defined $workbook) ||
	die "E when parsing xlsx\n";
    ($workbook->worksheet_count() == 1) ||
	die "E parsing xlsx: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($grexCol, $cohortCol) = (-1,-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	($cell->value() eq "grexomeID") &&
	    ($grexCol = $col);
	($cell->value() eq "pathology") &&
	    ($cohortCol = $col);
    }
    ($grexCol >= 0) ||
	die "E parsing xlsx: no column title is grexomeID\n";
    ($cohortCol >= 0) ||
	die "E parsing xlsx: no col title is pathology\n";
    
    foreach my $row ($rowMin+1..$rowMax) {
	my $grexome = $worksheet->get_cell($row, $grexCol)->value;
	# skip "none" lines
	($grexome eq "none") && next;
	(defined $sample2cohort{$grexome}) && 
	    die "E parsing xlsx: have 2 lines with grexome $grexome\n";
	my $cohort = $worksheet->get_cell($row, $cohortCol)->value;
	$sample2cohort{$grexome} = $cohort;
	$cohorts{$cohort} = 1;
    }
    @cohorts = sort(keys(%cohorts));
}


#########################################################
# parse stdin

# header: copy and add descriptions for new INFO fields
while(my $line = <STDIN>) {
    if ($line =~ /^##/) {
	print $line;
    }
    elsif ($line =~ /^#CHROM/) {
	# add new INFO field descriptions
	foreach my $cohort (@cohorts) {
	    foreach my $geno (@genoCategories) {
		print "##INFO=<ID=COUNT_$cohort"."_$geno,Number=.,Type=String,Description=\"Number of samples from the $cohort cohort having each genotype that fall in the $geno category, comma-separated\">\n";
	    }
	}

	# add info with full command line run
	my $com = qx/ps -o args= $$/;
	chomp($com);
	$com .= " < ".`readlink -f /proc/$$/fd/0` ;
	chomp($com);
	$com .= " > ".`readlink -f /proc/$$/fd/1` ;
	chomp($com);
	$com .= " 2> ".`readlink -f /proc/$$/fd/2` ;
	chomp($com);
	print "##countPerCohort=<commandLine=\"$com\">\n";

	# make sure this VCF has the correct GENO columns
	my $goodEnd = join("\t", ("FORMAT", @genoCategories))."\n";
	($line =~ /$goodEnd$/) ||
	    die "infile #CHROM line doesn't have correct GENO columns\n";
	# print #CHROM line
	print $line;
	last;
     }
    else {
	die "E: parsing header, found bad line:\n$line";
    }
}

# now parse data lines
while (my $line = <STDIN>) {
    chomp($line);

    my @data = split(/\t/, $line, -1);
    (@data == 13) || die "wrong number (".scalar(@data).") of fields in line:\n$line\n";
    # first 7 fields are just copied
    my $lineToPrint = shift(@data);
    foreach my $i (2..7) {
	$lineToPrint .= "\t".shift(@data);
    }
    # grab current INFO
    my $info = shift(@data);

    # remaining fields will be copied as well later

    # count the samples per geno and per cohort
    # key is eg COUNT_ovo_HR, value is what will be printed for this INFO field
    my %counts;
    foreach my $cohort (@cohorts) {
	foreach my $geno (@genoCategories) {
	    $counts{"COUNT_$cohort"."_$geno"} = "";
	}
    }

    foreach my $i (0..$#genoCategories) {
	foreach my $realGenoData (split(/\|/,$data[1+$i])) {
	    ($realGenoData =~ s/^([^~]+)~//) || die "cannot grab realGeno from $realGenoData\n";
	    my $realGeno = $1;
	    # count the number of samples from each cohort for this realGeno
	    # key is cohort
	    my %counts4realGeno = ();
	    foreach my $cohort (@cohorts) {
		$counts4realGeno{$cohort} = 0;
	    }
	    foreach my $sample (split(/,/,$realGenoData)) {
		(defined $sample2cohort{$sample}) ||
		    die "sample $sample is not present in the $metadata xls file\n";
		my $cohort = $sample2cohort{$sample} ;
		$counts4realGeno{$cohort}++ ; 
	    }

	    # OK, add counts for realGeno to %counts
	    foreach my $cohort (@cohorts) {
		# if no counts for this cohort, skip
		($counts4realGeno{$cohort}) || next;
		my $key = "COUNT_$cohort"."_".$genoCategories[$i] ;
		($counts{$key}) && ($counts{$key} .= ",");
		$counts{$key} .= "$realGeno:".$counts4realGeno{$cohort};
	    }
	}
    }

    # OK, add the data to INFO
    foreach my $cohort (@cohorts) {
	foreach my $geno (@genoCategories) {
	    ($counts{"COUNT_$cohort"."_$geno"}) && 
		($info .= ";COUNT_$cohort"."_$geno=".$counts{"COUNT_$cohort"."_$geno"});
	}
    }

    # clean up if INFO was just '.'
    $info =~ s/^\.;//;

    # done, print line
    $lineToPrint .= "\t$info\t".join("\t", @data)."\n";
    print $lineToPrint;
}
