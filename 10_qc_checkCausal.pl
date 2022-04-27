#!/usr/bin/perl

# 25/06/2021
# NTM

# QC script: look for severe biallelic variants affecting  the canonical
# transcript of each "causal gene" from the samples metadata file.


use strict;
use warnings;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use Getopt::Long;
use POSIX qw(strftime);

use lib "$RealBin";
use grexome_metaParse qw(parseSamples parseCandidateGenes);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## options / params from the command-line

# samples, pathologies and candidateGenes XLSX files, empty defaults
my ($samplesFile, $pathologies, $candidatesFiles) = ("","","");

# dir containing samples CSV files produced by grexome-TIMC-Secondary,
# default to current dir
my $inDir = ".";

# help: if true just print $USAGE and exit
my $help = '';


my $USAGE = "\nFor each sample that has a causal gene in the metadata XLSX, examine the SAMPLES 
results CSV files in indir and report if and how the causal gene's canonical transcript is hit.
Print our findings to stdout.
Search for samples not referenced in the metadata XLSX that might carry causal genes.
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--samplesFile : samples metadata xlsx file, with path
--pathologies string [optional] : pathologies metadata xlsx file, with path
--candidateGenes string [optional] : comma-separated list of xlsx files holding known candidate genes, with paths
--indir [$inDir] : dir containing samples CSV files produced by grexome-TIMC-Secondary
--help : print this USAGE";

GetOptions ("samplesFile=s" => \$samplesFile,
        "pathologies=s" => \$pathologies,
	    "candidateGenes=s" => \$candidatesFiles,
	    "indir=s" => \$inDir,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($samplesFile) || die "E $0: you must provide a samplesFile file\n$USAGE\n";
(-f $samplesFile) || die "E $0: the supplied samplesFile file doesn't exist\n";

(-d $inDir) ||
    die "E $0: inDir $inDir doesn't exist or isn't a directory\n$USAGE\n";
opendir(INDIR, $inDir) ||
    die "E $0: cannot opendir inDir $inDir\n";

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


#########################################################
# parse useful info from samples metadata file:
# hashref, key==sample id, value is the $cohort this sample belongs to (ie pathologyID)
my $sample2cohortR;
# hashref, key==sample id, value == HGNC gene name of the known causal gene
my $sample2causalR;

{
    my @parsed = &parseSamples($samplesFile);
    $sample2cohortR = $parsed[0];
    $sample2causalR = $parsed[3];
}

# $parseCandidateGenes actually returns a hashref, key==$cohort, value is
# a hashref whose keys are gene names and values are the "Confidence scores"
my $knownCandsR;
$knownCandsR = &parseCandidateGenes($candidatesFiles, $samplesFile, $pathologies);


#########################################################
# examine each sample CSV file

# total number of samples
my $nbSamples = scalar(keys(%$sample2cohortR));
# number of samples with a causal gene in metadata
my $nbCausal = scalar(keys(%$sample2causalR));
# arrays of strings to print for samples whose "causal gene" is hit by
# biallelic HIGH, biallelic MODHIGH, biallelic MODERATE/LOW, monoallelic,
# or no variant at all
my @causalHigh = ();
my @causalModHigh = ();
my @causalLow = ();
my @mono = ();
my @noVar = ();

# arrays of strings to print for samples affected by a causal gene not 
# identified in metadata hits by HV HIGH , HV MODHIGH or biallelic MODHIGH
my @candidateHVHigh = ();
my @candidateHVModHigh = ();
my @candidateBiallelicModHigh= ();

foreach my $inFile (sort(readdir(INDIR))) {
    ($inFile =~ /^\./) && next;
    my $sample;
    foreach my $s (keys(%$sample2cohortR)) {
	if ($inFile =~ /^$sample2cohortR->{$s}\.$s\./) {
	    $sample = $s;
	    last;
	}
    }
    ($sample) || 
	((print "W: inFile $inFile doesn't seem to correspond to any sample, skipping it\n") && next);

    my $patho = $sample2cohortR->{$sample};
    # empty $sample2cohortR as we go, for sanity testing and speed
    delete($sample2cohortR->{$sample});

    # process effective only for samples with causal genes in the metadata file
    if(defined($sample2causalR->{$sample})){
        my $causal = $sample2causalR->{$sample};

        # find lines affecting a CANONICAL transcript of gene $causal in $inFile
        open(INFILE, "$inDir/$inFile") ||
        die "E $0: cannot open inFile $inDir/$inFile\n";

        # indexes of columns of interest
        my ($symbolCol, $biallelCol, $canonCol) = (-1,-1,-1);
        # header: grab column indexes of interest
        my $header = <INFILE>;
        chomp($header);
        my @header = split(/\t/,$header);
        foreach my $i (0..$#header) {
        ($header[$i] eq 'SYMBOL') && ($symbolCol = $i);
        ($header[$i] eq 'BIALLELIC') && ($biallelCol = $i);
        ($header[$i] eq 'CANONICAL') && ($canonCol = $i);
        ($biallelCol > -1) && ($symbolCol > -1) && ($canonCol > -1) && last;
        }
        (($biallelCol > -1) && ($symbolCol > -1) && ($canonCol > -1)) || 
        die "E $0: cannot find required column headers in $inFile:\n$header\n";

        # find lines affecting a CANONICAL transcript of $causal in $inFile,
        # if no lines $foundCausal stays false
        my $foundCausal = 0;
        while (my $line = <INFILE>) {
        chomp($line);
        my @line = split(/\t/,$line);
        ($line[$symbolCol] eq "' $causal") || next;
        ($line[$canonCol] eq 'YES') || next;

        if ($line[$biallelCol] eq 'HIGH') {
            push(@causalHigh, "$patho\t$sample\t$causal\tHIGH");
        }
        elsif ($line[$biallelCol] eq 'MODHIGH') {
            push(@causalModHigh, "$patho\t$sample\t$causal\tMODHIGH");
        }
        elsif ($line[$biallelCol] eq 'NO') {
            push(@mono, "$patho\t$sample\t$causal\tMONO-ALLELIC");
        }
        else {
            push(@causalLow, "$patho\t$sample\t$causal\tLOW-MODER");
        }
        $foundCausal = 1;
        last;
        }
        close(INFILE);
        ($foundCausal) ||
        push(@noVar, "$patho\t$sample\t$causal\tNOVARIANT");
    }
    # searching causal genes in samples not listed in the metadata
    else{
        open(INFILE, "$inDir/$inFile") ||
        die "E $0: cannot open inFile $inDir/$inFile\n";

        # indexes of columns of interest 
        my ($symbolCol,$genotypeCol, $biallelCol, $canonCol) = (-1,-1,-1,-1);
        # header: grab column indexes of interest
        my $header = <INFILE>;
        chomp($header);
        my @header = split(/\t/,$header);
        foreach my $i (0..$#header) {
        ($header[$i] eq 'SYMBOL') && ($symbolCol = $i);
        ($header[$i] eq 'GENOTYPE') && ($genotypeCol = $i);
        ($header[$i] eq 'BIALLELIC') && ($biallelCol = $i);
        ($header[$i] eq 'CANONICAL') && ($canonCol = $i);
        ($biallelCol > -1) && ($genotypeCol > -1) && ($symbolCol > -1) && ($canonCol > -1) && last;
        }
        ($biallelCol > -1) && ($genotypeCol > -1) && ($symbolCol > -1) && ($canonCol > -1) || 
        die "E $0: cannot find required column headers in $inFile:\n$header\n";


        # Store the name of the gene with a variation already stored in the 
        # tables @candidateHVHigh or @candidateHVModHigh or @candidateBiallelicModHigh
        my $GeneAlreadySeen="";

        # variants counter for the same gene for the same sample.
        # initialized to 1 => used when a variation is retained 
        my $varNB=1;

        # variation genotype 
        # filtering if HV_HIGH present we do not want HV_MODHIGH or Biallelic_MODHIGH for the same patient.
        # status: 0 -> initialization, 1-> HV_HIGH observed, 2-> HV_MODHIGH,-1-> BIALLELIC_MODHIGH 
        # TODO make a better implementation (eg arrray of variations?)
        my $genoSeen=0;
        
        while (my $line = <INFILE>) {
        chomp($line);
        my @line = split(/\t/,$line);

        #selection variations affecting only canonical genes
        ($line[$canonCol] eq 'YES') || next;

        # removal of "' GeneName" formatting for "GeneName"
        # necessary for comparisons with the GeneName contained in knownCandsR
        my $gene=substr $line[$symbolCol], 2;

        # pathology and gene exist in the hash ref knownCandsR 
        # otherwise next variant
        (defined($knownCandsR->{$patho}->{$gene})) || next;
        my $score=$knownCandsR->{$patho}->{$gene};

        if ($line[$genotypeCol] eq 'HV') {
            my $varStatus="$line[$genotypeCol]_$line[$biallelCol]";
            if ($line[$biallelCol] eq 'HIGH'){
                if ($gene ne $GeneAlreadySeen) {
                    push(@candidateHVHigh,"$patho:$score\t$sample\t$gene\t$varStatus\t$varNB");
                    $GeneAlreadySeen=$gene;
                    $genoSeen=1;
                }
                else{
                    ($genoSeen==1) || die "E : previous variants for $sample on $gene have a different genotype than HV_HIGH.$genoSeen Please check.\n";
                    $varNB+=1;
                    $candidateHVHigh[-1]="$patho:$score\t$sample\t$gene\t$varStatus\t$varNB";
                }
            }
            elsif (($line[$biallelCol] eq 'MODHIGH') && ($genoSeen != 1)) {
                if ($gene ne $GeneAlreadySeen) {
                    push(@candidateHVModHigh,"$patho:$score\t$sample\t$gene\t$varStatus\t$varNB");
                    $GeneAlreadySeen=$gene;
                    $genoSeen=2;
                }
                else{
                    ($genoSeen==2) || die "E : previous variants for $sample on $gene have a different genotype than HV_MODHIGH.$genoSeen Please check.\n";
                    $varNB+=1;
                    $candidateHVModHigh[-1]="$patho:$score\t$sample\t$gene\t$varStatus\t$varNB";
                }
            }
            else{
                next;
            } 
        }
        elsif (($line[$genotypeCol] eq 'HET') && ($line[$biallelCol] eq 'MODHIGH') && ($genoSeen <= 0)){
            my $varStatus="BIALLELIC_$line[$biallelCol]";
            if ($gene ne $GeneAlreadySeen) {
                push(@candidateBiallelicModHigh,"$patho:$score\t$sample\t$gene\t$varStatus\t$varNB");
                $GeneAlreadySeen=$gene;
                $genoSeen=-1;
            }
            else{
                ($genoSeen==-1) ||  die "E : previous variants for $sample on $gene have a different genotype than HET_MODHIGH.$genoSeen Please check.\n";
                $varNB+=1;
                $candidateBiallelicModHigh[-1]="$patho:$score\t$sample\t$gene\t$varStatus\t$varNB";
            }
        }
        else{
            next;
        } 
        }
        close(INFILE);
    }
}

closedir(INDIR);

# sanity: every sample from metadata should have been seen, assuming
# we didn't analyze a subcohort/incomplete infile
(keys(%$sample2cohortR)) &&
    (print "W: ".scalar(keys(%$sample2cohortR))." samples from metadata have no samples.csv files\n");

print "Total number of samples in metadata: $nbSamples\n";
print "Samples with a causal gene in metadata: $nbCausal\n\n";
print "Examining the samples CSV files in $inDir (canonical transcripts only), we found:\n";
print "samples whose causal gene is BIALLELIC HIGH: ".scalar(@causalHigh)."\n";
print "samples whose causal gene is BIALLELIC MODHIGH: ".scalar(@causalModHigh)."\n";
print "samples whose causal gene is BIALLELIC MODERATE or LOW: ".scalar(@causalLow)."\n";
print "samples whose causal gene is only hit by MONOALLELIC: ".scalar(@mono)."\n";
print "samples whose causal gene has NO VARIANT: ".scalar(@noVar)."\n";
print "samples carrying a causal gene but not previously observed: ".(scalar(@candidateHVHigh)+scalar(@candidateHVModHigh)+scalar(@candidateBiallelicModHigh))."\n";
print"\n";

if (@causalHigh) {
    print "List of samples whose causal gene is BIALLELIC HIGH:\n";
    print join("\n", @causalHigh);
    print "\n\n";
}
if (@causalModHigh) {
    print "List of samples whose causal gene is BIALLELIC MODHIGH:\n";
    print join("\n", @causalModHigh);
    print "\n\n";
}
if (@causalLow) {
    print "List of samples whose causal gene is BIALLELIC LOW or MODERATE:\n";
    print join("\n", @causalLow);
    print "\n\n";
}
if (@mono) {
    print "List of samples whose causal gene is only hit by ONE MONOALLELIC variant:\n";
    print join("\n", @mono);
    print "\n\n";
}
if (@noVar) {
    print "List of samples whose causal gene has NO VARIANT:\n";
    print join("\n", @noVar);
    print "\n\n";
}

if (@candidateHVHigh){
    print "List of samples samples affected by a causal gene is HV HIGH not identified in metadata:\n";
    print join("\n", @candidateHVHigh);
    print "\n\n";
}

if (@candidateHVModHigh){
    print "List of samples samples affected by a causal gene is HV MODHIGH not identified in metadata:\n";
    print join("\n", @candidateHVModHigh);
    print "\n\n";
}

if (@candidateBiallelicModHigh){
    print "List of samples samples affected by a causal gene is BIALLELIC MODHIGH not identified in metadata:\n";
    print join("\n", @candidateBiallelicModHigh);
    print "\n\n";
}

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";