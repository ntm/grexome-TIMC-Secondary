#!/usr/bin/perl

# 07/11/2019
# NTM

# Read on stdin a table of transcripts, as produced by 
# makeCanonicalTranscriptTable.pl .
# Print to stdout the same info in BED format, one line per exon:
# chr start end name
# For name I use: $transcript_$exonNum , $exonNum starts at 1 
# in transcript order (ie taking into account the strand)

# skip header
<STDIN>;

while(my $line = <STDIN>) {
    chomp $line;
    my ($transcript,$gene,$ensg,$chr,$strand,$cdsS,$cdsE,$exonS,$exonE) = split(/\t/,$line);

    my @exonsS = split(/,/,$exonS);
    my @exonsE = split(/,/,$exonE);
    (@exonsS == @exonsE) || die "E: different numbers of exon starts and ends in:\n$line\n";

    foreach my $i (0..$#exonsS) {
        # if on '+' strand increment so we start at 1
        my $exonNum = $i + 1;
        # if transcript is on - strand we count down
        ($strand eq '-') && ($exonNum = scalar(@exonsS) - $i);
        print "$chr\t$exonsS[$i]\t$exonsE[$i]\t$transcript"."_$exonNum\n";
    }
}
