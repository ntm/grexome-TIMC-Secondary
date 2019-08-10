#!/usr/bin/perl

# 25/03/2018
# NTM

# Take as argument a $metadata xlsx file, and an $outDir that doesn't
# exist, and read on stdin a fully annotated TSV file;
# make $outDir and create in it one TSV file per cohort.
# The cohorts are defined in $metadata.
# For each sample, any identified causal (mutation in a) gene is grabbed
# from $metadata.
# @notControls defined at the top of this script says which cohorts
# should NOT be used as negative controls for each other.
#
# Changes from infile:
# - for each cohort, the COUNTs for all cohorts except those defined
#   in @notControls are summed and printed in new columns
#   COUNT_NEGCTRL_* (one per GENO).
# - COUNTs are printed in the order: HV HET OTHER HR.
# - the HR GENO column is removed.
# - each outFile lists only the samples that belong to the relevant
#   cohort in the remaining GENO columns, and samples from neg control
#   cohorts are now listed in new NEGCTRL_* columns placed immediately
#   after the HV/HET/OTHER columns.
# - skip lines where no samples from the cohort are HV|HET|OTHER with
#   the allele ALLELE_NUM (which this line is dealing with).
# - if a sample has an identified causal gene, this sample is not COUNTed
#   as HV or HET in its own cohort except for mutations hitting said gene
#   (ie the HV/HET COUNT is decremented if sample is HV/HET, but the sample
#   remains in the HV/HET GENO column); the sample is still counted as NEGCTRL
#   for other cohorts.


use strict;
use warnings;
use Spreadsheet::XLSX;


# @notControls: array of arrayrefs, each arrayref holds cohorts that
# should NOT be used as neg controls for each other.
# The cohort names must match the "pathology" column of the $metadata xlsx
# (this is checked).
my @notControls = (["Flag","Astheno"],
		   ["Azoo","Ovo"],
		   ["Globo","Macro","Terato"]);


(@ARGV == 2) || die "needs two args: a metadata xlsx and a non-existing outDir\n";
my ($metadata, $outDir) = @ARGV;

(-e $outDir) && 
    die "found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "cannot mkdir outDir $outDir\n";

#########################################################
# parse metadata file

# key==sample id, value is the $cohort this sample belongs to
my %sample2cohort = ();
# cohort names
my @cohorts = ();
# causal gene, key==sample id, value == HGNC gene name
my %sample2causal = ();
# for sanity: if $sample2causal{$sample} is defined, $sample2causalSeen{$sample}
# will be set to 1 when we see the gene. At the end if a causal gene was never seen
# there is probably a typo in the xlsx.
my %sample2causalSeen = ();


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
    my ($grexCol, $cohortCol,$causalCol) = (-1,-1,-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	($cell->value() eq "grexomeID") &&
	    ($grexCol = $col);
	($cell->value() eq "pathology") &&
	    ($cohortCol = $col);
	($cell->value() eq "Causal gene") &&
	    ($causalCol = $col);
     }
    ($grexCol >= 0) ||
	die "E parsing xlsx: no column title is grexomeID\n";
    ($cohortCol >= 0) ||
	die "E parsing xlsx: no col title is pathology\n";
    ($causalCol >= 0) ||
	die "E parsing xlsx: no col title is Causal gene\n";
    
    foreach my $row ($rowMin+1..$rowMax) {
	my $grexome = $worksheet->get_cell($row, $grexCol)->value;
	# skip "none" lines
	($grexome eq "none") && next;
	(defined $sample2cohort{$grexome}) && 
	    die "E parsing xlsx: have 2 lines with grexome $grexome\n";
	my $cohort = $worksheet->get_cell($row, $cohortCol)->value;
	$sample2cohort{$grexome} = $cohort;
	$cohorts{$cohort} = 1;
	my $causal = $worksheet->get_cell($row, $causalCol)->value;
	if ($causal) {
	    $sample2causal{$grexome} = $causal;
	    $sample2causalSeen{$grexome} = 0;
	}
    }
    @cohorts = sort(keys(%cohorts));
}


#########################################################
# check @notControls cohort names and store in %notControls hash

# %notControls: key is a cohort name, value is a hashref
# with keys == cohorts that shouldn't be used as negative 
# controls for this cohort, value==1
my %notControls = ();

foreach my $notConR (@notControls) {
    foreach my $cohort (@$notConR) {
	(grep($cohort eq $_, @cohorts)) ||
	    die "E in extractCohorts: cohort $cohort from notControls is not in cohorts @cohorts\n";
	(defined $notControls{$cohort}) &&
	    die "E in extractCohorts: cohort $cohort from notControls present twice?!\n";
	my %notCtmp;
	foreach my $notC (@$notConR) {
	    ($notC eq $cohort) && next;
	    $notCtmp{$notC} = 1;
	}
	$notControls{$cohort} = \%notCtmp;
    }
}

#########################################################

# hash of filehandles open for writing, one for each cohort
my %outFHs;
foreach my $cohort (@cohorts) {
    my $outFile = "$outDir/$cohort.csv";
    open (my $FH, "> $outFile") || die "cannot open $outFile for writing";
    $outFHs{$cohort} = $FH ;
}


# key == cohort, value == ref to array of ints (one int per column),
# value 0 at index $i means sum it for COUNT_NEGCTRL_HV,
# value 1 means sum it for COUNT_NEGCTRL_HET,
# value 2 means sum it for COUNT_NEGCTRL_OTHER,
# value 3 means sum it for COUNT_NEGCTRL_HR,
# value 4 means drop column $i in this cohort's outfile,
# value 5 means keep it as-is and print before COUNTs,
# value 6 means keep it as-is but print after COUNTs,
# value 7 means col has HV samples, need to separate samples from this cohort and NEGCTRLs,
#   and decrement COUNT_$cohort_HV for samples from this cohort with causal genes
# value 8 means col has HET samples, need to separate samples from this cohort and NEGCTRLs,
#   and decrement COUNT_$cohort_HET for samples from this cohort with causal genes
# value 9 means col has OTHER samples, need to separate samples from this cohort and NEGCTRLs
# value 10 means col has COUNT_$cohort_HV
# value 11 means col has COUNT_$cohort_HET
# value 12 means col has COUNT_$cohort_OTHER
# value 13 means col has COUNT_$cohort_HR
my %keepIndexes;


# header: we want an array with the header titles
my $header = <STDIN>;
chomp($header);
my @headers = split(/\t/, $header);


# indexes of the ALLELE_NUM and SYMBOL columns 
my ($alleleNumIndex, $symbolIndex);
foreach my $i (0..$#headers) {
    if ($headers[$i] eq "ALLELE_NUM") {
	$alleleNumIndex = $i;
    }
    elsif ($headers[$i] eq "SYMBOL") {
	$symbolIndex = $i;
    }
}
($alleleNumIndex) || die "could not find ALLELE_NUM in headers\n";
($symbolIndex) || die "could not find SYMBOL in headers\n";

# fill %keepIndexes and print new headers
foreach my $cohort (@cohorts) {
    my @keepIndexes;
    my $toPrint = "";
    # we assume all COUNT columns are adjacent.
    # Assumption is checked with $countsSeen: value 0 means not seen yet, 
    # 1 means we are in COUNTs, -1 means we already left and printed COUNTs
    my $countsSeen = 0;
    # the COUNT headers for $cohort are stored in @countsToPrint at index (10-$keepIndex)
    my @countsToPrint = ();
    foreach my $title (@headers) {
	if (! $toPrint) {
	    # we always want to keep the first column and 
	    # this simplifies things (no \t)
	    $toPrint = $title;
	    push(@keepIndexes,5);
	}
	elsif ($title =~ /^COUNT_(\w+)_([^_]+)$/) {
	    my ($curCoh,$curGeno) = ($1,$2);
	    ($countsSeen==-1) && 
		die "Parsing headers, found COUNT but we already left and printed them\n";
	    $countsSeen = 1;
	    if ($curCoh eq $cohort) {
		if ($curGeno eq "HV") {
		    $countsToPrint[0] = $title ;
		    push(@keepIndexes,10);
		}
		elsif ($curGeno eq "HET") {
		    $countsToPrint[1] = $title ;
		    push(@keepIndexes,11);
		}
		elsif ($curGeno eq "OTHER") {
		    $countsToPrint[2] = $title ;
		    push(@keepIndexes,12);
		}
		elsif ($curGeno eq "HR") {
		    $countsToPrint[3] = $title ;
		    push(@keepIndexes,13);
		}
		else {
		    die "reading headers, found COUNT $title but invalid geno $curGeno\n";
		}
	    }
	    # for other cohorts don't print anything, but fill keepIndexes
	    elsif (defined ${$notControls{$cohort}}{$curCoh}) {
		# don't use this as neg control, just skip
		push(@keepIndexes,4);
	    }
	    # otherwise this is a NEGCTRL, look at geno
	    elsif ($curGeno eq "HV") {
		push(@keepIndexes,0);
	    }
	    elsif ($curGeno eq "HET") {
		push(@keepIndexes,1);
	    }
	    elsif ($curGeno eq "OTHER") {
		push(@keepIndexes,2);
	    }
	    elsif ($curGeno eq "HR") {
		push(@keepIndexes,3);
	    }
	    else {
		die "WTF can't happen";
	    }
	}
	elsif ($title =~ /^COUNT_/) {
	    die "title is a COUNT but I didnt parse it: $title\n";
	}
	elsif ($title eq "HR") {
	    # skip
	    push(@keepIndexes,4);
	}
	elsif ($title eq "HV") {
	    ($countsSeen == 1) && 
		die "HV header comes right after COUNTs, code doesnt deal with this.\n";
	    $toPrint .= "\t$title\tNEGCTRL_$title";
	    push(@keepIndexes,7);
	}
	elsif ($title eq "HET") {
	    ($countsSeen == 1) && 
		die "HET header comes right after COUNTs, code doesnt deal with this.\n";
	    $toPrint .= "\t$title\tNEGCTRL_$title";
	    push(@keepIndexes,8);
	}
	elsif ($title eq "OTHER") {
	    ($countsSeen == 1) && 
		die "OTHER header comes right after COUNTs, code doesnt deal with this.\n";
	    $toPrint .= "\t$title\tNEGCTRL_$title";
	    push(@keepIndexes,9);
	}
	else {
	    # keep all other fields but use $countsSeen to decide if it's before or after COUNTs
	    if ($countsSeen==0) {
		$toPrint .= "\t$title";
		push(@keepIndexes,5);
	    }
	    else {
		if ($countsSeen==1) {
		    # first column after COUNTs
		    $countsSeen = -1;
		    # print the COUNTs headers in new order
		    $toPrint .= "\t".join("\t",@countsToPrint);
		    # print the new NEGCTRL headers
		    $toPrint .= "\tCOUNT_NEGCTRL_HV\tCOUNT_NEGCTRL_HET\tCOUNT_NEGCTRL_OTHER\tCOUNT_NEGCTRL_HR";
		}
		$toPrint .= "\t$title";
		push(@keepIndexes,6);
	    }
	}
    }
    # sanity:
    (@keepIndexes == @headers) || 
	die "keepIndexes doesn't have same number of elements as headers:\n@keepIndexes\n@headers\n";
    print { $outFHs{$cohort} } "$toPrint\n";
    $keepIndexes{$cohort} = \@keepIndexes ;
    # make sure the COUNT columns were printed
    ($countsSeen == -1) ||
	die "the COUNT_* columns were not seen or are last in the infile, this isn't supported.\n";
}


# Parse data lines
while (my $line = <STDIN>) {
    chomp($line);
    my @fields = split(/\t/, $line, -1) ;
    # grab alleleNum and symbol now, they're the same for all cohorts
    my $alleleNum = $fields[$alleleNumIndex];
    my $symbol = $fields[$symbolIndex];

    foreach my $cohort (@cohorts) {
	# we will prepare $toPrintBeforeCounts and $toPrintAfterCounts, and
	# print them with @countsToPrint and @negCtrlSums in between

	# initialize $toPrintBeforeCounts with first field
	my $toPrintBeforeCounts = "$fields[0]";
	my $toPrintAfterCounts = "";

	# negCtrlSums will sum the NEGCTRL samples for each geno;
	# sums for HV, HET, OTHER, HR are stored at indexes 0-3 (same value 
	# used in @keepIndexes). Initialize with 4 zeroes
	my @negCtrlSums = (0) x 4;

	# the COUNTs for $cohort will be printed in order HV, HET, OTHER, HR,
	# store values in array at indexes 0-3 (@keepIndex value - 10)
	my @countsToPrint = (0) x 4;

	# we only keep lines where at least one sample from $cohort is HV/HET/OTHER
	# with allele $alleleNum (after adjusting for diagnosed samples)
	my $keepLine = 0;
	foreach my $i (1..$#fields) {
	    if ($keepIndexes{$cohort}->[$i] <= 3) {
		$negCtrlSums[$keepIndexes{$cohort}->[$i]] += $fields[$i] ;
	    }

	    elsif ($keepIndexes{$cohort}->[$i] == 4) {
		next;
	    }
	    elsif ($keepIndexes{$cohort}->[$i] == 5) {
		$toPrintBeforeCounts .= "\t".$fields[$i];
	    }
	    elsif ($keepIndexes{$cohort}->[$i] == 6) {
		$toPrintAfterCounts .= "\t".$fields[$i];
	    }

	    elsif ($keepIndexes{$cohort}->[$i] == 7) {
		# realGenos (with samples) to print for this GENO category
		my @toPrintGenos = ();
		# same but for negative controls
		my @toPrintGenosBad = ();
		foreach my $realGenoData (split(/\|/,$fields[$i])) {
		    ($realGenoData =~ s/^([^~]+)~//) || die "cannot grab realGeno from $realGenoData\n";
		    my $realGeno = $1;
		    ($realGeno =~ /^$alleleNum\//) || ($realGeno =~ /\/$alleleNum$/) || next;
		    my @goodSamples = ();
		    my @badSamples = (); # neg ctrls
		    foreach my $sample (split(/,/,$realGenoData)) {
			if ($sample2cohort{$sample} eq $cohort) {
			    push(@goodSamples,$sample);
			}
			elsif (! defined ${$notControls{$cohort}}{$sample2cohort{$sample}}) {
			    push(@badSamples,$sample);
			}
		    }
		    if (@goodSamples) {
			push(@toPrintGenos, $realGeno."~".join(',',@goodSamples));
		    }
		    if (@badSamples) {
			push(@toPrintGenosBad, $realGeno."~".join(',',@badSamples));
		    }
		}
		(@toPrintGenos) && ($keepLine=1);
		$toPrint .= "\t".join('|',@toPrintGenos);
		$toPrint .= "\t".join('|',@toPrintGenosBad);
	    }
	    elsif (($keepIndexes{$cohort}->[$i] >= 10) && ($keepIndexes{$cohort}->[$i] <= 13)) {
		$countsToPrint[$keepIndexes{$cohort}->[$i] - 10] = $fields[$i] ;
	    }
	    else {
		die "bad keepIndex value $keepIndexes{$cohort}->[$i], impossible\n";
	    }
	}

	my $toPrint = "$toPrintBeforeCounts\t".join("\t", @countsToPrint)."\t".join("\t", @negCtrlSums)."$toPrintAfterCounts\n";

	($keepLine) && (print { $outFHs{$cohort} } "$toPrint\n");
    }
}

foreach my $fh (values %outFHs) {
    close($fh);
}
