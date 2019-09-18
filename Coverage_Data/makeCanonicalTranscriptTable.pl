#!/usr/bin/env perl

# 16/09/2019
# NTM


# Parse the Ensembl gtf file on stdin, and the list of canonical 
# transcripts from Ensembl contained in $canoTransFile (produced
# by getCanonicalTranscripts.pl)
# Print to stdout a tab-delimited file with columns:
# "TRANSCRIPT\tGENE\tCHROM\tSTRAND\tCDS_START\tCDS_END\tEXON_STARTS\tEXON_ENDS\n"
# containing one line per canonical transcript, 
# STRAND is + or -,
# by convention if transcript is non-coding CDS_START==CDS_END==1,
# EXON_* columns contain comma-separated lists of coordinates, sorted
# numerically, with matching START and END coordinates satisfying
# START < END.

use strict;
use warnings;


# take one arg: filename with ensembl canonical transcripts
(@ARGV == 1) || 
    die "E: needs one arg, a file listing the canonical transcripts\n";

my $ensembl_canonical = shift(@ARGV);

(-f $ensembl_canonical) ||
    die "E: argument $ensembl_canonical is not a file\n";


##############################################################
# store canonical transcripts: key==transcript_id, value==1
my %canonical = () ;

open(CANON, "gunzip -c $ensembl_canonical |") || 
    die "cannot gunzip-open $ensembl_canonical for reading\n" ;

while(my $line = <CANON>) {
    chomp $line;
    $canonical{$line} = 1 ;
}
close(CANON);


##############################################################

# skip header lines
foreach my $i (1..5) {
    my $line = <>;
    ($line =~ /^#!/) || die "E: skipping header line but its not a header?\n$line";
}    

# fill following data structures, each is a hash with key==$transcript and:
# value == "TRANSCRIPT\tGENE\tCHROM\tSTRAND"
my %trans2printFirst;
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

    # TRANSCRIPT
    ($fields[8] =~ /transcript_id "(ENST\d+)";/) ||
	die "E: cannot grab transcript_id from line:\n$line\n";
    my $transcript = $1;

    # skip if it's not a canonical transcript
    ($canonical{$transcript}) || next;

    # STRAND
    my $strand = $fields[6];

    if (! defined($trans2printFirst{$transcript})) {
	# GENE
	($fields[8] =~ /gene_name "([^"]+)";/) ||
	    die "E: cannot grab gene_name from line:\n$line\n";
	my $gene = $1;
	# CHROM: we use chr* convention
	my $chr = $fields[0];
	($chr eq "MT") && ($chr = "M");
	$chr = "chr$chr";

	$trans2printFirst{$transcript} = "$transcript\t$gene\t$chr\t$strand";
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

# print header
print "TRANSCRIPT\tGENE\tCHROM\tSTRAND\tCDS_START\tCDS_END\tEXON_STARTS\tEXON_ENDS\n";

# print data in random order
foreach my $transcript (keys(%trans2printFirst)) {
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
    # for sanity checks later:
    delete($canonical{$transcript});
}

# sanity?
# Commenting out, we exepct thousands of canonical transcripts that
# don't belong to chr1-22 or X Y M (eg they map to ALT chromosomes),
# these aren't in the Ensembl GTF.
# foreach my $transcript (sort keys(%canonical)) {
#     warn "W: $transcript is canonical but was never seen\n";
# }
