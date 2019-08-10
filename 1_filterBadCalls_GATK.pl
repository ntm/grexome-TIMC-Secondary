#!/usr/bin/perl


# Parses on stdin a VCF file with one or more sample data columns;
# takes as args:
# --minDP int (must have DP >= $minDP),
# --minGQ int (must have GQ >= $minGQ), 
# --strandDisc [NOT IMPLEMENTED YET], 
# --minFracVarReads float (must have fraction of variant reads for the 
#                          called variant >= $minFracVarReads, filter is not
#                          applied if called genotype is VAR1/VAR2 with 2 non-ref alleles.
# Prints to stdout a similar VCF file where the variant calls in data columns
# are replaced by ./. if a condition is not met,
# and lines where every sample is now ./. are skipped.

use strict;
use warnings;
use Getopt::Long;

# default values for args
my $minDP = 0;
my $minGQ = 0;
my $minFracVarReads = 0;

GetOptions ("minDP=i" => \$minDP,
            "minGQ=i"   => \$minGQ,
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
    # $keepLine: boolean, true if at least one sample is not ./. after filtering
    my $keepLine = 0;
    my @data = split(/\t/, $line);
    (@data >= 10) || die "no sample data in line?\n$line\n";
    # first 9 fields are just copied, but grab FORMAT string
    my $lineToPrint = shift(@data);
    foreach my $i (2..8) {
	$lineToPrint .= "\t".shift(@data);
    }
    my $format = shift(@data);
    $lineToPrint .= "\t$format";
    # %format: key is a FORMAT key (eg GQ), value is the index of that key in $format
    my %format;
    { 
	my @format = split(/:/, $format);
	foreach my $i (0..$#format) {
	    $format{$format[$i]} = $i ;
	}
    }
    # sanity: make sure the fields we need are there
    (defined $format{"GQ"}) || die "no GQ key in FORMAT string for line:\n$line\n";
    (defined $format{"DP"}) || die "no DP key in FORMAT string for line:\n$line\n";
    (defined $format{"GT"}) || die "no GT key in FORMAT string for line:\n$line\n";
    (defined $format{"AD"}) || die "no AD key in FORMAT string for line:\n$line\n";

    # now deal with actual data fields
    while(my $data = shift(@data)) {
	# if genotype is already ./. == NOCALL, just copy the call (dropping any 
	# other FORMAT values)
	($data =~ m~^\./\.~) && ($lineToPrint .= "\t./.") && next;
	# otherwise examine content and apply filters
	my @thisData = split(/:/, $data) ;

	if ($thisData[$format{"GQ"}] < $minGQ) {
	    # GQ too low, change to NOCALL
	    $lineToPrint .= "\t./.";
	    next;
	}

	if (($thisData[$format{"DP"}] eq '.') || ($thisData[$format{"DP"}] < $minDP)) {
	    # DP too low, change to NOCALL
	    $lineToPrint .= "\t./.";
	    next;
	}
	
	# in any case we need DP > 0 (otherwise how can we make a call? and
	# this causes illegal division by zero when calculating minFracVarReads)
	if ($thisData[$format{"DP"}] == 0) {
	    # if there was a call it's bullshit, change to NOCALL
	    $lineToPrint .= "\t./.";
	    next;
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
	    my $fracVarReads = $ads[$geno2] / $thisData[$format{"DP"}] ;
	    if ($fracVarReads < $minFracVarReads) {
		# fracVarReads too low, change to NOCALL
		$lineToPrint .= "\t./.";
		next;
	    }
	}
	# else this GT is VAR1/VAR2, don't try to filter

	# other filters (eg strandDisc) would go here

	# OK data passed all filters, print and set $keepLine
	$lineToPrint .= "\t$data";
	$keepLine = 1;
    }
    # done with $line, print if at least one sample is not ./.
    ($keepLine) && (print "$lineToPrint\n");
}


