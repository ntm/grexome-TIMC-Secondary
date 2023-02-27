#!/usr/bin/env perl

# 16/09/2019
# NTM

# Takes a single arg: a string defining the "transcripts of interest", currently
# 'canon' (-ical) or 'mane' (-Select)
#
# Parse the Ensembl gtf file on stdin.
# Print to stdout a tab-delimited file with columns
# "TRANSCRIPT\tGENE\tENSG\tCHROM\tSTRAND\tCDS_START\tCDS_END\tEXON_STARTS\tEXON_ENDS\n"
# containing one line per transcript of interest:
# GENE is the ENSG id if we can't find a gene name,
# STRAND is + or -,
# by convention if transcript is non-coding CDS_START==CDS_END==1,
# EXON_* columns contain comma-separated lists of coordinates, sorted
# numerically, with matching START and END coordinates satisfying
# START < END.

use strict;
use warnings;


(@ARGV == 1) ||
    die "E: needs one arg, a string defining transcripts of interest\n";

my $TOIs = shift(@ARGV);

# precise string to search for in last data column if Ensembl GTF
my $toiDef;
if ($TOIs eq 'canon') {
    $toiDef = 'tag "Ensembl_canonical";'
}
elsif ($TOIs eq 'mane') {
    $toiDef = 'tag "MANE_Select";'
}
else {
    die "E: unsuppported type of transcripts: must be 'canon' or 'mane', not $TOIs\n";
}

##############################################################

# skip header lines
foreach my $i (1..5) {
    my $line = <>;
    ($line =~ /^#!/) || die "E: skipping header line but its not a header?\n$line";
}    

# fill following data structures, each is a hash with key==$transcript and:
# value == "TRANSCRIPT\tGENE\tENSG\tCHROM\tSTRAND"
my %trans2printFirst;
# value == $chrom (for sorting the output)
my %trans2chr;
# value == arrayref \($cdsStart,$cdsEnd)
my %trans2cds;
# value == arrayref with 2 strings \($exonStarts,$exonEnds), the strings
# are comma-separated and already sorted
my %trans2exons;

while (my $line = <>) {
    chomp($line);
    my @fields = split(/\t/, $line, -1);
    (@fields == 9) || die "E: line doesn't have 9 columns:\n$line\n";

    # we only want exon, start_codon and stop_codon lines
    my $lineType = $fields[2];
    ($lineType eq "exon") || ($lineType eq "start_codon") || ($lineType eq "stop_codon") ||
	next;

    # ignore if not a TOI
    ($fields[8] =~ /$toiDef/) || next;

    # TRANSCRIPT
    ($fields[8] =~ /transcript_id "(ENST\d+)";/) ||
	die "E: cannot grab transcript_id from line:\n$line\n";
    my $transcript = $1;

    # STRAND
    my $strand = $fields[6];

    if (! defined($trans2printFirst{$transcript})) {
	# ENSG
	($fields[8] =~ /gene_id "([^"]+)";/) ||
	    die "E: cannot grab gene_id from line:\n$line\n";
	my $ensg = $1;
	# GENE
	my $gene = $ensg;
	($fields[8] =~ /gene_name "([^"]+)";/) &&
	    ($gene = $1);
	# CHROM
	my $chr = $fields[0];
	$chr =~ s/^chr//;
	# for sorting we want just the chrom number, replace X Y M|MT by 23-25
	if ($chr eq "X") { $trans2chr{$transcript} = "23" }
	elsif ($chr eq "Y") { $trans2chr{$transcript} = "24" }
	elsif (($chr eq "M") || ($chr eq "MT")) { $trans2chr{$transcript} = "25" }
	else { $trans2chr{$transcript} = $chr }
	# for printing we use chr* convention
	($chr eq "MT") && ($chr = "M");
	$chr = "chr$chr";

	$trans2printFirst{$transcript} = "$transcript\t$gene\t$ensg\t$chr\t$strand";
	$trans2cds{$transcript} = [];
	$trans2exons{$transcript} = ["",""];
    }

    if ($lineType eq "exon") {
	# we always have start < end but exons of transcritps on - strand
	# appear in reverse genomic order in the file
	if ($strand eq "+") {
	    $trans2exons{$transcript}->[0] .= "$fields[3],";
	    $trans2exons{$transcript}->[1] .= "$fields[4],";
	}
	else {
	    $trans2exons{$transcript}->[0] = "$fields[3],$trans2exons{$transcript}->[0]";
	    $trans2exons{$transcript}->[1] = "$fields[4],$trans2exons{$transcript}->[1]";
	}
    }
    elsif ($lineType eq "start_codon") {
	if ($strand eq "+") {
	    $trans2cds{$transcript}->[0] = $fields[3];
	}
	else {
	    $trans2cds{$transcript}->[1] = $fields[4];
	}
    }
    elsif ($lineType eq "stop_codon") {
	if ($strand eq "+") {
	    $trans2cds{$transcript}->[1] = $fields[4];
	}
	else {
	    $trans2cds{$transcript}->[0] = $fields[3];
	}
    }
}

# sub for sorting by chrom, then by coord of first exon, then by the list of
# all starts-stops, then by transcriptID (just to get something deterministic)
sub byChrByCoord {
    if ($trans2chr{$a} <=> $trans2chr{$b}) {
	return($trans2chr{$a} <=> $trans2chr{$b});
    }
    # else... start of first exon
    ($trans2exons{$a}->[0] =~ /^(\d+),/) ||
	die "E: in byChrByCoord cannot extract first exon start from $trans2exons{$a}->[0]\n";
    my $startA = $1;
    ($trans2exons{$b}->[0] =~ /^(\d+),/) ||
	die "E: in byChrByCoord cannot extract first exon start from $trans2exons{$b}->[0]\n";
    my $startB = $1;
    if ($startA <=> $startB) {
	return($startA <=> $startB);
    }
    # else... end of first exon
    ($trans2exons{$a}->[1] =~ /^(\d+),/) ||
	die "E: in byChrByCoord cannot extract first exon end from $trans2exons{$a}->[1]\n";
    my $endA = $1;
    ($trans2exons{$b}->[1] =~ /^(\d+),/) ||
	die "E: in byChrByCoord cannot extract first exon end from $trans2exons{$b}->[1]\n";
    my $endB = $1;
    if ($endA <=> $endB) {
	return($endA <=> $endB);
    }
    # else... stringwise compare the full lists of exon coords
    if (my $coordsCmp = (join("__",@{$trans2exons{$a}}) cmp join("__",@{$trans2exons{$b}}))) {
	return($coordsCmp);
    }
    # else stringwise compare transcriptIDs, these are different for sure
    return($a cmp $b);
}


# print header
print "TRANSCRIPT\tGENE\tENSG\tCHROM\tSTRAND\tCDS_START\tCDS_END\tEXON_STARTS\tEXON_ENDS\n";

# print data in chrom-coord order
foreach my $transcript (sort byChrByCoord keys(%trans2printFirst)) {
    print $trans2printFirst{$transcript};
    # CDS*
    if ((defined $trans2cds{$transcript}->[0]) && (defined $trans2cds{$transcript}->[1])) {
	print "\t".$trans2cds{$transcript}->[0];
	print "\t".$trans2cds{$transcript}->[1];
    }
    elsif (defined $trans2cds{$transcript}->[0]) {
	# transcript has 5' or 3' incomplete CDS, extend to last exon end
	print "\t".$trans2cds{$transcript}->[0];
	($trans2exons{$transcript}->[1] =~ /^(\d+),$/) ||
	    ($trans2exons{$transcript}->[1] =~ /,(\d+),$/) ||
	    die "E: cannot extract last exon end for $transcript: $trans2exons{$transcript}->[1]\n";
	print "\t$1";
    }
    elsif (defined $trans2cds{$transcript}->[1]) {
	# transcript has 5' or 3' incomplete CDS, extend from first exon start
	($trans2exons{$transcript}->[0] =~ /^(\d+),/) ||
	    die "E: cannot extract first exon start for $transcript: $trans2exons{$transcript}->[0]\n";
	print "\t$1";
	print "\t".$trans2cds{$transcript}->[1];
    }
    else {
	# no start nor stop codons, assume non-coding transcript, CDS* == 1 1 by convention
	# NOTE: this could be a transcript that's CDS-incomplete on both ends,
	# nevermind ignore it
	print "\t1\t1";
    }
    # EXONS*: chop final comma
    (chop($trans2exons{$transcript}->[0]) eq ',') ||
	die "E: chopped final comma from exon_starts but it's not a comma! $transcript $trans2exons{$transcript}->[0]\n";
    (chop($trans2exons{$transcript}->[1]) eq ',') ||
	die "E: chopped final comma from exon_ends but it's not a comma! $transcript $trans2exons{$transcript}->[1]\n";
    print "\t".$trans2exons{$transcript}->[0];
    print "\t".$trans2exons{$transcript}->[1];
    print "\n";
}
