#!/usr/bin/perl

# 25/03/2018
# NTM

# Take 2 arguments: $inDir $outDir
# $inDir must contain cohort TSVs as produced by extractCohorts.pl;
# $outDir doesn't exist, it will be created and filled with one TSV
# per sample.
# For a sample, we only print lines where this sample has a non-HR 
# genotype: this genotype is printed in a new column GENOTYPE, 
# inserted right after ALT.
# Also the HV, NEGCTRL_HV, HET etc... columns are not printed.
# For every cohort, we must have a $cohort.list file in current
# dir with the list of samples.
#
# NOTE: the fact that ALT is the 3rd column, and that HV et al
# are the last 6 columns, is HARD-CODED!
# If this isn't true the script dies.

use strict;
use warnings;


(@ARGV == 2) || die "needs 2 args: an inDir and a non-existant outDir\n";
my ($inDir, $outDir) = @ARGV;
(-d $inDir) ||
    die "inDir $inDir doesn't exist or isn't a directory\n";
opendir(INDIR, $inDir) ||
    die "cannot opendir inDir $inDir\n";
(-e $outDir) && 
    die "found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "cannot mkdir outDir $outDir\n";


while (my $inFile = readdir(INDIR)) {
    ($inFile =~ (/^(\w+)\.csv/)) || next;
    my $cohort = $1;

    open(IN, "$inDir/$inFile") ||
	die "cannot open cohort datafile $inDir/$inFile\n";
    my $header = <IN>;
    chomp($header);
    ($header =~ s/^(POSITION\tREF\tALT\t)/$1GENOTYPE\t/) ||
	die "cannot add GENOTYPE to header of inFile $inFile\n";
    ($header =~ s/\tHV\tNEGCTRL_HV\tHET\tNEGCTRL_HET\tOTHER\tNEGCTRL_OTHER$//) ||
	die "cannot remove HV HET OTHER from header of inFile $inFile\n";
    # hash of filehandles open for writing, one for each sample
    # from this cohort
    # key is the sample id
    my %outFHs;
    (-f "$cohort.list") || die "cannot find file $cohort.list\n";
    open(COHORT, "$cohort.list") ||
	die "cannot open cohort file $cohort.list\n" ;
    while (my $sample = <COHORT>) {
	chomp($sample);
	my $outFile = "$outDir/$cohort.$sample.csv";
	open (my $FH, "> $outFile") || die "cannot open $outFile for writing";
	print $FH "$header\n";
	$outFHs{$sample} = $FH ;
    }
    close(COHORT);

    # now read the data
    while (my $line = <IN>) {
	chomp($line);
	my @fields = split(/\t/, $line, -1) ;
	my $toPrintStart = join("\t",@fields[0..2])."\t";
	my $toPrintEnd = join("\t",@fields[3..($#fields-6)])."\n";

	foreach my $i ($#fields-5,$#fields-3,$#fields-1) {
	    foreach my $realGenoData (split(/\|/,$fields[$i])) {
		($realGenoData =~ s/^([^~]+)~//) || 
		    die "cannot grab realGeno from $realGenoData\n";
		my $realGeno = $1;
		# replace / with \ for excel :-( 
		$realGeno =~ s~/~\\~ ;
		foreach my $sample (split(/,/,$realGenoData)) {
		    print { $outFHs{$sample} } "$toPrintStart$realGeno\t$toPrintEnd" ;
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
