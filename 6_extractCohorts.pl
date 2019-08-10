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
# For each $cohort, the GENO columns HV/HET/OTHER/HR are modified as follows:
# - the HR GENO column is removed (but samples are COUNTed, see below).
# - we make new NEGCTRL_* columns placed immediately after the HV/HET/OTHER columns.
# - NEGCTRL_* columns list all samples falling in that GENO category from other
#   cohorts, except those defined in @notControls.
# - if a sample from $cohort has an identified $causalGene, this sample is
#   removed from the HV/HET columns, except in lines where SYMBOL==$causalGene.
#
# New COUNT_$cohort_$geno and COUNT_NEGCTRL_$geno columns are created
# for each  GENO (HV, HET, OTHER, HR) in that order.
# These columns contain the total number of samples listed in the
# corresponding GENO column (except for HR, which has no GENO column).
# For HR we count all samples (ie don't care about @notControls or $causalGene).
# The COUNTs are inserted right after the SYMBOL column.
#
# Lines where no samples from the cohort are HV|HET (for this alt allele)
# are skipped. We rely on the fact that vcf2tsv.pl moved HV/HET genotypes
# concerning other alleles to OTHER (but we check it).


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
# for sanity: $causalSeen{$gene} is set to 0 for every causal gene,
# value changed to 1 when we see that gene. At the end if a causal 
# gene was never seen there is probably a typo in the xlsx.
my %causalSeen = ();


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
	if ($worksheet->get_cell($row, $causalCol)) {
	    my $causal = $worksheet->get_cell($row, $causalCol)->value;
	    # clean up a bit, remove leading or trailing whitespace
	    $causal =~ s/^\s+//;
	    $causal =~ s/\s+$//;
	    $sample2causal{$grexome} = $causal;
	    $causalSeen{$causal} = 0;
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

#########################################################
# headers

my $header = <STDIN>;
chomp($header);
my @headers = split(/\t/, $header);

# the genotype categories (don't change this)
my @genoCategories = ("HV","HET","OTHER","HR");

# useful columns: SYMBOL and genos in @genoCategories order
my $symbolCol;
my @genoCols;
foreach my $i (0..$#headers) {
    ($headers[$i] eq "SYMBOL") && ($symbolCol = $i);
    foreach my $gi (0..$#genoCategories) {
	($headers[$i] eq $genoCategories[$gi]) && 
	    ($genoCols[$gi] = $i);
    }
}
($symbolCol) || die "could not find SYMBOL in headers\n";
foreach my $gi (0..$#genoCategories) {
    ($genoCols[$gi]) || die "cound not find $genoCategories[$gi] in headers\n";
    ($genoCols[$gi] == ($genoCols[0]+$gi)) || 
	die "E: GENO columns are not subsequent and in genoCategories order, the code relies on this\n";
    # actually maybe the code just requires that they be 
    # consecutive columns and start with HV (not sure)
}

# print new headers
foreach my $cohort (@cohorts) {
    # we always want to keep the first column and 
    # this simplifies things (no \t)
    my $toPrint = "$headers[0]";
    foreach my $i (1..$#headers) {
	if ($i == $symbolCol) {
	    $toPrint .= "\t$headers[$i]";
	    # COUNTs go right after SYMBOL
	    foreach my $geno (@genoCategories) {
		$toPrint .= "\tCOUNT_$cohort"."_$geno";
	    }
	    # COUNT_NEGCTRL_* right after
	    foreach my $geno (@genoCategories) {
		$toPrint .= "\tCOUNT_NEGCTRL_$geno";
	    }
	}
	elsif (($i == $genoCols[0]) || ($i == $genoCols[1]) || ($i == $genoCols[2])) {
	    # HV/HET/OTHER
	    $toPrint .= "\t$headers[$i]\tNEGCTRL_$headers[$i]";
	}
	elsif ($i == $genoCols[3]) {
	    # HR is not printed
	}
	else {
	    # all other columns are kept as-is
	    $toPrint .= "\t$headers[$i]";
	}
    }
    print { $outFHs{$cohort} } "$toPrint\n";
}

#########################################################
# Parse data lines

while (my $line = <STDIN>) {
    chomp($line);
    my @fields = split(/\t/, $line, -1) ;

    # $symbol doesn't depend on cohorts
    my $symbol = $fields[$symbolCol];
    (defined $causalSeen{$symbol}) && ($causalSeen{$symbol}=1);

    
  COHORT:
    foreach my $cohort (@cohorts) {
	# build array of 8 counts for $cohort: HV,HET,OTHER,HR and again for NEGCTRLs
	my @counts = (0) x 8;
	# also build array of 6 GENO columns: HV,NEGCTRL_HV,HET,NEGCTRL_HET,OTHER,NEGCTRL_OTHER
	my @genos = ("") x 6;

	# parse data
	foreach my $gi (0..3) {
	    my @genoData = split(/\|/,$fields[$genoCols[$gi]]);
	    # sanity: at most one genotype except for OTHER column
	    (@genoData <= 1) || ($gi==2) || 
		die "E: more than one genoData for genotype $genoCols[$gi], impossible. Line:\n$line\n";
	    foreach my $genoData (@genoData) {
		($genoData =~ /^(\d+\/\d+)~([^~\|]+)$/) ||
		    die "E: cannot parse GENOS data $genoData in line:\n$line\n";
		# $geno is the genotype (eg 1/1 or 0/2)
		my $geno = $1;
		my @samples = split(/,/,$2);
		# @goodSamples will hold samples that should be counted for $cohort
		# @badSamples will hold samples that should be counted as NEGCTRLs for $cohort
		my @goodSamples = ();
		my @badSamples = ();
		foreach my $sample (@samples) {
		    if ($sample2cohort{$sample} eq $cohort) {
			# $sample belongs to cohort
			if (($gi >= 2) || (! defined $sample2causal{$sample}) || ($sample2causal{$sample} eq $symbol)) {
			    # we are not in HV|HET or sample has no causal gene or it's the current gene
			    push(@goodSamples,$sample);
			}
			# else ignore this sample == NOOP
		    }
		    elsif (($gi==3) || (! defined ${$notControls{$cohort}}{$sample2cohort{$sample}})) {
			# sample is from another cohort that can be used as control 
			# for $cohort, or we are in HR (where all samples are counted)
			push(@badSamples,$sample);
		    }
		    # else sample is from a different cohort but it's in @notControls: NOOP
		}
		
		# OK, store counts and GENOs (careful with indexes)
		if (@goodSamples) {
		    $counts[$gi] += scalar(@goodSamples);
		    if ($gi < 3) {
			# don't store the GENO for HR ie gi==3
			($genos[$gi * 2]) && ($genos[$gi * 2] .= '|');
			$genos[$gi * 2] .= "$geno~".join(',',@goodSamples);
		    }
		}
		if (@badSamples) {
		    $counts[$gi + 4] += scalar(@badSamples);
		    if ($gi < 3) {
			($genos[1 + $gi * 2]) && ($genos[1 + $gi * 2] .= '|');
			$genos[1 + $gi * 2] .= "$geno~".join(',',@badSamples);
		    }
		}
	    }

	    # if we just finished with HV and HET but there's no sample in either,
	    # we can skip this line in this cohort
	    if (($gi == 1) && ($counts[0] == 0) && ($counts[1] == 0)) {
		next COHORT;
	    }
	    # otherwise: parse OTHER and HR data
	}

	# OK we have some data to print for $cohort
	# always keep first field 
	my $toPrint = "$fields[0]";
	foreach my $i (1..$#fields) {
	    if ($i == $symbolCol) {
		$toPrint .= "\t$fields[$i]";
		# print all COUNTs
		$toPrint .= "\t".join("\t",@counts);
	    }
	    elsif ($i == $genoCols[0]) {
		# HV -> print all 6 GENO columns
		# NOTE: we rely on the fact that the GENOs are consecutive and start with HV
		$toPrint .= "\t".join("\t",@genos);
	    }
	    elsif (! grep(/^$i$/, @genoCols)) {
		$toPrint .= "\t$fields[$i]";
	    }
	}
	print { $outFHs{$cohort} } "$toPrint\n";
    }
}

# sanity
foreach my $gene (keys(%causalSeen)) {
    ($causalSeen{$gene}) ||
	warn "W: causal gene $gene mentioned in $metadata was never seen, most likely a typo in $metadata\n";
}

foreach my $fh (values %outFHs) {
    close($fh);
}
