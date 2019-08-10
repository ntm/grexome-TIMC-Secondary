#!/usr/bin/perl

# NTM
# 09/07/2019 (but starting from 1_filterBadCalls.pl which is older)


# Parses on stdin a Strelka GVCF file with one or more 
# sample data columns;
# takes as args:
# --minDP int (must have DP or DPI >= $minDP),
# --minGQX int (must have GQX >= $minGQ), 
# --strandDisc [NOT IMPLEMENTED YET], 
# --minFracVarReads float (must have fraction of variant reads for the 
#                          called variant AD[i]/DP >= $minFracVarReads, filter is not
#                          applied if called genotype is HR or VAR1/VAR2 with 2 non-ref alleles.
# Prints to stdout a VCF file where:
# - non-variant lines are removed;
# - the variant calls in data columns are replaced by ./. if a condition
#   is not met, or if previous call was '.' (the Strelka NOCALL);
# - lines where every sample is now ./. or 0/0 are skipped.


use strict;
use warnings;
use Getopt::Long;

# default values for args
my $minDP = 0;
my $minGQX = 0;
my $minFracVarReads = 0;

GetOptions ("minDP=i" => \$minDP,
            "minGQX=i"   => \$minGQX,
            "minFracVarReads=f"  => \$minFracVarReads)
    or die("Error in command line arguments\n");

# parse header, just copy it
while(my $line = <STDIN>) {
    if ($line =~ /^##/) {
	print $line;
    }
    elsif ($line =~ /^#CHROM/) {
	# add info with full command line run
	my $com = qx/ps -o args= $$/;
	chomp($com);
	$com .= " < ".`readlink -f /proc/$$/fd/0` ;
	chomp($com);
	$com .= " > ".`readlink -f /proc/$$/fd/1` ;
	chomp($com);
	$com .= " 2> ".`readlink -f /proc/$$/fd/2` ;
	chomp($com);
	print "##filterBadCalls=<commandLine=\"$com\">\n";
	print $line;
	last;
    }
    else {
	die "E: parsing header, found bad line:\n$line";
    }
}

# now parse data lines
while(my $line = <STDIN>) {
    chomp($line);
    # $keepLine: boolean, true if at least one sample is not ./. and 0/0 after filtering
    my $keepLine = 0;
    my @data = split(/\t/, $line);
    (@data >= 10) || die "no sample data in line?\n$line\n";
    # if no ALT in line, skip immediately
    ($data[4] eq '.') && next;
    # first 9 fields are just copied, but grab FORMAT string
    my $lineToPrint = shift(@data);
    foreach my $i (2..8) {
	$lineToPrint .= "\t".shift(@data);
    }
    my $format = shift(@data);
    $lineToPrint .= "\t$format";
    # %format: key is a FORMAT key (eg GQX), value is the index of that key in $format
    my %format;
    { 
	my @format = split(/:/, $format);
	foreach my $i (0..$#format) {
	    $format{$format[$i]} = $i ;
	}
    }
    # sanity: make sure the fields we need are there
    # we use either DP or DPI as long as one is present and not '.'
    (defined $format{"GQX"}) || die "no GQX key in FORMAT string for line:\n$line\n";
    (defined $format{"DP"}) || (defined $format{"DPI"}) ||
	die "no DP or DPI key in FORMAT string for line:\n$line\n";
    (defined $format{"GT"}) || die "no GT key in FORMAT string for line:\n$line\n";
    (defined $format{"AD"}) || die "no AD key in FORMAT string for line:\n$line\n";

    # now deal with actual data fields
    while(my $data = shift(@data)) {
	# if genotype is already '.' or './.' == NOCALL, just use ./.
	if (($data =~ m~^\.$~) || ($data =~ m~^\.:~) || ($data =~ m~^\./\.~)) {
	    $lineToPrint .= "\t./." ;
	    next;
	}
	# otherwise examine content and apply filters
	my @thisData = split(/:/, $data) ;

	if ((! $thisData[$format{"GQX"}]) || ($thisData[$format{"GQX"}] eq '.') ||
	    ($thisData[$format{"GQX"}] < $minGQX)) {
	    # GQX undefined or too low, change to NOCALL
	    $lineToPrint .= "\t./.";
	    next;
	}

	# grab the depth (DP or DPI, whichever is defined and higher)
	my $thisDP = -1;
	if ((defined $format{"DP"}) && ($thisData[$format{"DP"}]) && ($thisData[$format{"DP"}] ne '.')) {
	    $thisDP = $thisData[$format{"DP"}];
	}
	if ((defined $format{"DPI"}) && ($thisData[$format{"DPI"}]) && ($thisData[$format{"DPI"}] ne '.') &&
	    ($thisData[$format{"DPI"}] > $thisDP)) {
	    $thisDP = $thisData[$format{"DPI"}];
	}

	# if depth too low or undefined for this sample, change to NOCALL
	if ($thisDP < $minDP) {
	    $lineToPrint .= "\t./.";
	    next;
	}
	# with GATK I had some issues with DP=0 calls, causing illegal divisions
	# by zero when calculating fracVarReads, but that is now skipped above

	# Strelka makes some phased calls sometimes, homogenize
	$thisData[$format{"GT"}] =~ s~\|~/~ ;
	# Strelka also makes some hemizygous calls (eg when the position
	# is in a HET deletion), makes sense but still, homogenize
	$thisData[$format{"GT"}] =~ s~^(\d+)$~$1/$1~;

	# minFracVarReads doesn't apply to HR
	if ($thisData[$format{"GT"}] ne '0/0') {
	    # AD should always be there when something other than HR was called
	    if ((! $thisData[$format{"AD"}]) || ($thisData[$format{"AD"}] =~ /^[\.,]+$/)) {
		die "E: GT is not 0/0 but we don't have AD or AD data is blank in:\n$line\nright after:\n$lineToPrint\n";
	    }

	    # for minFracVarReads we need GT for the called variant index in AD
	    my ($geno1,$geno2) = split(/\//, $thisData[$format{"GT"}]);
	    ((defined $geno1) && (defined $geno2)) ||
		die "E: a sample's genotype cannot be split: ".$thisData[$format{"GT"}]."in:\n$line\n";
	    if (($geno1 == 0) || ($geno2 == 0) || ($geno1 == $geno2)) {
		# ok we can deal with this
		($geno2 == 0) && ($geno2 = $geno1);
		# in all cases $geno2 is now the index of the VAR
		my @ads = split(/,/, $thisData[$format{"AD"}]);
		my $fracVarReads = $ads[$geno2] / $thisDP ;
		if ($fracVarReads < $minFracVarReads) {
		    # fracVarReads too low, change to NOCALL
		    $lineToPrint .= "\t./.";
		    next;
		}
	    }
	    # else this GT is VAR1/VAR2, don't try to filter on fracVarReads
	}

	# other filters (eg strandDisc) would go here

	# OK data passed all filters, print and set $keepLine if not HOMOREF
	$lineToPrint .= "\t$data";
	($thisData[$format{"GT"}] ne '0/0') && ($keepLine = 1);
    }
    # done with $line, print if at least one sample is not NOCALL|HOMOREF
    ($keepLine) && (print "$lineToPrint\n");
}


