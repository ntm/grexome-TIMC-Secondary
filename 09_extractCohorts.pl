#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2025
#
# This file is part of grexome-TIMC-Secondary, written by Nicolas Thierry-Mieg
# (CNRS, France) Nicolas.Thierry-Mieg@univ-grenoble-alpes.fr
#
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.
############################################################################################


# 25/03/2018
# NTM

# Takes as arguments a $samplesFile xlsx file, a $pathologies xlsx file, a
# comma-separated list of xlsx candidateGenes files, an $outDir and
# a $tmpDir that don't exist; 
# reads on stdin a fully annotated TSV file;
# makes $outDir and creates in it one gzipped TSV file per cohort.
# The cohorts are defined in $pathologies and/or $samplesFile.
# For each sample, any identified causal (mutation in a) gene 
# is grabbed from $samplesFile.
# "Compatible" cohorts (ie belonging to a common compatibility group in the
# pathologies XLSX file) appear in COMPAT columns and are NOT used
# as negative controls for each other.
# For optimal performance $tmpDir should be on a RAMDISK (eg tmpfs).
#
# A new KNOWN_CANDIDATE_GENE column is inserted right after SYMBOL:
# it holds a comma-separated list of "$pathologyID:$score" strings, listing
# each pathologyID for which SYMBOL is a known candidate or causal gene with
# confidence score $score, as produced by &parsePathologies().
#
# The HV/HET/OTHER/HR columns are removed and used to produce the following
# columns in each $cohort outfile, in this order and starting where HV was:
# - $cohort_HV, $cohort_HET -> samples from $cohort (except those having
#   a causal variant in another gene)
# - $cohort_OTHERCAUSE_HV, $cohort_OTHERCAUSE_HET -> samples from $cohort
#   that have a causal variant in another gene
# - COMPAT_HV, COMPAT_HET -> samples from compatible cohorts
# - NEGCTRL_HV, NEGCTRL_HET -> samples from all other cohorts
# - OTHERGENO -> samples with OTHER genotypes (from any cohort)
#
# In addition, new COUNT_* columns are created for each of the above columns,
# in the same order, followed by a single COUNT_HR column (counting HR samples
# from any cohort without caring about $compatible or $causalGene).
# The COUNTs are inserted right after the new KNOWN_CANDIDATE_GENE column.
#
# Lines where no samples from the cohort are HV|HET (for this alt allele)
# are skipped. We rely on the fact that vcf2tsv.pl moved HV/HET genotypes
# concerning other alleles to OTHER (but we check it).


use strict;
use warnings;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use Getopt::Long;
use POSIX qw(strftime);
use Parallel::ForkManager;

use lib "$RealBin";
use grexome_metaParse qw(parsePathologies parseSamples parseCandidateGenes);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## hard-coded stuff that shouldn't change much

# adaptive strategy parameters for $batchSize: without --batchSize, $batchSize will
# be adjusted in order to process ($jobs-1) batches in [$batchTimeLow,$batchTimeHigh] seconds
my $batchSizeAdaptive = 1;
# aim for between 2m and 10m for ($jobs-1) batches
my ($batchTimeLow, $batchTimeHigh) = (120,600);


#############################################
## options / params from the command-line

# samples, pathologies and candidateGenes XLSX files, empty defaults
my ($samplesFile, $pathologies, $candidatesFiles) = ("","","");

# outDir and tmpDir, also no defaults
my ($outDir, $tmpDir);

# number of jobs
my $jobs = 16;

# max number of lines to read in a single batch. Each batch is then
# processed by a worker thread. This is a performance tuning param,
# higher values shold be faster but increase RAM requirements.
# We have used values between 10k and 100k.
# Providing it with --batchSize hard-sets $batchSize, could be useful if you lack RAM.
# Without --batchSize we start with a low value and use an adaptive strategy
my $batchSize;

# help: if true just print $USAGE and exit
my $help = '';


my $USAGE = "\nParse on STDIN a fully annotated TSV file as produced by steps 1-8 of this secondaryAnalysis pipeline; create in outDir one gzipped TSV file per cohort.\n
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--samples string : samples metadata xlsx file, with path
--pathologies string [optional] : pathologies metadata xlsx file, with path
--candidateGenes string [optional] : comma-separated list of xlsx files holding known candidate genes, with paths
--outdir string : subdir where resulting cohort files will be created, must not pre-exist
--tmpdir string : subdir where tmp files will be created (on a RAMDISK if possible), must not pre-exist and will be removed after execution
--jobs [$jobs] : number of parallel jobs=threads to run
--batchSize [adaptive] : size of each batch, lower decreases RAM requirements and performance,
             defaults to an optimized adaptive strategy, you shouldn\'t need to specify this except
             if you are running out of RAM (suggested reasonable values: between 10000 and 100000)
--help : print this USAGE";

GetOptions ("samples=s" => \$samplesFile,
            "pathologies=s" => \$pathologies,
            "candidateGenes=s" => \$candidatesFiles,
            "outdir=s" => \$outDir,
            "tmpdir=s" => \$tmpDir,
            "jobs=i" => \$jobs,
            "batchSize=i" => \$batchSize,
            "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($samplesFile) || die "E $0: you must provide a samples file\n";
(-f $samplesFile) || die "E $0: the supplied samples file doesn't exist\n";

($outDir) || die "E $0: you must provide an outDir\n";
(-e $outDir) && 
    die "E $0: found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "E $0: cannot mkdir outDir $outDir\n";

($tmpDir) || die "E $0: you must provide a tmpDir\n";
(-e $tmpDir) && 
    die "E $0: found argument $tmpDir but it already exists, remove it or choose another name.\n";
mkdir($tmpDir) || die "E $0: cannot mkdir tmpDir $tmpDir\n";

if ($jobs <= 2) {
    #  need one thread for eatTmpFiles and at least one worker thread
    warn "W $0: you set jobs=$jobs but we need at least 2 jobs, setting jobs=2 and proceeding\n";
    $jobs = 2;
}

if ($batchSize) {
    # disable adaptive strategy
    $batchSizeAdaptive = 0;
}
else {
    # default inital value
    $batchSize = 20000;
}

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


#########################################################
# parse all provided metadata files

# %knownCandidateGenes: key==$gene, value is a comma-separated sorted list
# of "$patho:$score" pairs
my %knownCandidateGenes;

{
    # $parseCandidateGenes actually returns a hashref, key==$cohort, value is
    # a hashref whose keys are gene names and values are the "Confidence scores"
    my $knownCandsR;
    if ($pathologies) {
        $knownCandsR = &parseCandidateGenes($candidatesFiles, $samplesFile, 0, $pathologies);
    }
    else {
        $knownCandsR = &parseCandidateGenes($candidatesFiles, $samplesFile, 0);
    }

    foreach my $c (sort keys(%$knownCandsR)) {
        foreach my $gene (keys(%{$knownCandsR->{$c}})) {
            my $score = $knownCandsR->{$c}->{$gene};
            if ($knownCandidateGenes{$gene}) {
                $knownCandidateGenes{$gene} .= ",$c:$score";
            }
            else {
                $knownCandidateGenes{$gene} = "$c:$score";
            }
        }
    }
}


# If $pathologies was provided, populate $compatibleR
# $compatibleR: hashref, key is a cohort name, value is a hashref
# with keys == cohorts that are compatible with this cohort, value==1
my $compatibleR;
if ($pathologies) {
    (-f $pathologies) || die "E $0: the supplied pathologies file $pathologies doesn't exist\n";
    $compatibleR = &parsePathologies($pathologies);
}
# else $compatibleR stays undef

# parse useful info from samples metadata file:
# hashref, key==sample id, value is the $cohort this sample belongs to (ie pathologyID)
my $sample2cohortR;
# hashref, causal gene, key==sample id, value == HGNC gene name
my $sample2causalR;

{
    my @parsed;
    if ($pathologies) {
        @parsed = &parseSamples($samplesFile, $pathologies);
    }
    else {
        @parsed = &parseSamples($samplesFile);
    }
    $sample2cohortR = $parsed[0];
    $sample2causalR = $parsed[3];
}

# cohort names, sorted and non-redundant
my @cohorts;
{
    my %seen;
    foreach my $c (values (%$sample2cohortR)) {
        $seen{$c} = 1;
    }
    @cohorts = sort(keys(%seen));
}

#########################################################

# array of filehandles open for writing, one for each cohort, same indexes as @cohorts
my @outFHs;
foreach my $cohorti (0..$#cohorts) {
    # sanity check (eg spaces in a cohort name break things)
    ($cohorts[$cohorti] =~ /^\w+$/) || die "E $0: cohort names must be alphanum only\n";
    my $outFile = "$outDir/$cohorts[$cohorti].csv.gz";
    open (my $FH, "| gzip -c > $outFile") || die "E $0: cannot gzip-open $outFile for writing";
    $outFHs[$cohorti] = $FH ;
}

#########################################################
# headers

my $header = <STDIN>;
chomp($header);
my @headers = split(/\t/, $header);

# useful columns from infile: SYMBOL and genotypes
my $symbolCol;
# key is one of "HV","HET","OTHER","HR", value is the column index in infile
my %genoCols;
foreach my $i (0..$#headers) {
    ($headers[$i] eq "SYMBOL") && ($symbolCol = $i);
    foreach my $geno ("HV","HET","OTHER","HR") {
        ($headers[$i] eq $geno) && ($genoCols{$geno} = $i);
    }
}
($symbolCol) || die "E $0: could not find SYMBOL in headers\n";
foreach my $geno ("HV","HET","OTHER","HR") {
    ($genoCols{$geno}) || die "E $0: cound not find $geno in headers\n";
}

# print new headers
foreach my $cohorti (0..$#cohorts) {
    # we always want to keep the first column (chr:coord) and 
    # this simplifies things (no \t)
    my $toPrint = "$headers[0]";
    foreach my $i (1..$#headers) {
        if ($i == $symbolCol) {
            $toPrint .= "\t$headers[$i]";
            # KNOWN_CANDIDATE_GENE and COUNTs go right after SYMBOL
            $toPrint .= "\tKNOWN_CANDIDATE_GENE";
            foreach my $group ($cohorts[$cohorti], $cohorts[$cohorti]."_OTHERCAUSE",
                               "COMPAT", "NEGCTRL") {
                foreach my $geno ("HV", "HET") {
                    $toPrint .= "\tCOUNT_$group"."_$geno";
                }
            }
            # OTHERGENO and HR get a single count for all cohorts
            $toPrint .= "\tCOUNT_OTHERGENO";
            $toPrint .= "\tCOUNT_HR";
        }
        elsif ($i == $genoCols{"HV"}) {
            # HV gets replaced by the new cohort-genotype categories
            foreach my $group ($cohorts[$cohorti], $cohorts[$cohorti]."_OTHERCAUSE",
                               "COMPAT", "NEGCTRL") {
                foreach my $geno ("HV", "HET") {
                    $toPrint .= "\t$group"."_$geno";
                }
            }
            # OTHERGENOs go in a single column for all cohorts, HRs are not listed
            $toPrint .= "\tOTHERGENO";
        }
        elsif (($i == $genoCols{"HET"}) || ($i == $genoCols{"OTHER"}) || ($i == $genoCols{"HR"})) {
            # NOOP, all sample lists taken care of above
        }
        else {
            # all other columns are kept as-is
            $toPrint .= "\t$headers[$i]";
        }
    }
    print { $outFHs[$cohorti] } "$toPrint\n";
    # flush filehandle before starting our eatTmpFiles job
    $outFHs[$cohorti]->flush();
}

#########################################################
# Parse data lines in batches of $batchSize lines, in parallel

# create fork manager
my $pm = new Parallel::ForkManager($jobs);

# $childFailed will become non-zero if at least one child died
my $childFailed = 0;
# Set up a callback so the parent knows if a child dies
$pm->run_on_finish( sub { ($_[1]) && ($childFailed=1) });

# spawn a child process that waits for workers to finish producing batches,
# and prints the tmpfiles to @outFHs in correct order, cleaning up behind 
# itself. 
if (! $pm->start) {
    &eatTmpFiles($tmpDir,\@cohorts,\@outFHs);
    $pm->finish;
}


# boolean flag, true iff current batch is the last
my $lastBatch = 0;
# number of the current batch
my $batchNum = 0;
# timestamp for batchSizeAdaptive strategy
my $timestampAdaptive = time();

while (!$lastBatch) {
    if ($childFailed) {
        $now = strftime("%F %T", localtime);
        die "E $now: $0 FAILED - some child died, no point going on\n";
    }
    $batchNum++;
    if (($batchSizeAdaptive) && ($batchNum % ($jobs-1) == 0) && ($batchNum >= $jobs)) {
        # jobs-1 because we want #workerThreads, excluding the eatTmpFiles job
        # $batchNum >= $jobs so we don't adjust on first batch of workers
        my $newTime = time();
        my $elapsed = $newTime - $timestampAdaptive;
        if ($elapsed < $batchTimeLow) {
            # increase batchSize by factor (1.2 * $btLow / $elapsed)
            # avoid divByZero
            $elapsed += 1;
            $batchSize = int(1.2 * $batchSize * $batchTimeLow / $elapsed);
            $now = strftime("%F %T", localtime);
            warn "I $now: $0 - batchNum=$batchNum, adjusting batchSize up to $batchSize\n";
        }
        elsif ($elapsed > $batchTimeHigh) {
            # decrease batchSize by factor $btHigh / (1.2 * $elapsed) 
            $batchSize = int($batchSize * $batchTimeHigh / $elapsed / 1.2);
            $now = strftime("%F %T", localtime);
            warn "I $now: $0 - batchNum=$batchNum, adjusting batchSize down to $batchSize\n";
        }
        $timestampAdaptive = $newTime;
    }

    my @lines = ();
    foreach my $i (1..$batchSize) {
        if (my $line = <STDIN>) {
            chomp($line);
            push(@lines,$line);
        }
        else {
            # no more lines
            $lastBatch = 1;
            last;
        }
    }

    # let worker threads take it from there
    $pm->start && next;
    # NOTE: IF YOU CHANGE the tmp filenames below ($tmpOut, $tmpOutFlag),
    # you MUST EDIT &eatTmpFiles()

    # create tmp output filehandles for this batch
    my @tmpOutFHs = ();
    foreach my $i (0..$#cohorts) {
        my $tmpOut = "$tmpDir/$batchNum.$cohorts[$i].tsv";
        open(my $outFH, "> $tmpOut") || die "E $0: cannot open $tmpOut for writing\n";
        push(@tmpOutFHs, $outFH);
    }
    # process this batch
    &processBatch(\@lines,\%knownCandidateGenes,$sample2cohortR,\@cohorts,
                  $sample2causalR,$compatibleR,$symbolCol,\%genoCols,\@tmpOutFHs);

    # done, close tmp FHs and create flag-file
    foreach my $outFH (@tmpOutFHs) {
        close($outFH) || die "E $0: cannot close tmp outFH $outFH\n";
    }
    my $tmpOutFlag = "$tmpDir/$batchNum.done";
    open(OUTFLAG, "> $tmpOutFlag") || die "E $0: cannot open flagfile $tmpOutFlag for writing\n";
    print OUTFLAG "$batchNum\n";
    close(OUTFLAG);
    $pm->finish;
}

# some children are still processing batches, but we know the last
# batchNum that will ever exist, tell &eatTmpFiles() so it can exit
# (of course if you change $tmpOutLast you have to edit &eatTmpFiles)
my $tmpOutLast = "$tmpDir/lastBatch";
open(OUTLAST, "> $tmpOutLast") || die "E $0: cannot open tmp-last-file $tmpOutLast for writing\n";
print OUTLAST "$batchNum\n";
close OUTLAST;

$pm->wait_all_children;

if ($childFailed) {
    $now = strftime("%F %T", localtime);
    die "E $now: $0 FAILED\n";
}

foreach my $fh (@outFHs) {
    close($fh);
}

$now = strftime("%F %T", localtime);
rmdir($tmpDir) || 
    die "E $now: $0 - all done but cannot rmdir tmpDir $tmpDir, why? $!\n";
warn "I $now: $0 - ALL DONE, completed successfully!\n";



#############################################
## subs

# process a batch of lines
# args:
# - ref to array of chomped lines
# - ref to %knownCandidateGenes
# - refs to %sample2cohort, to @cohorts, and to %sample2causal
# - ref to %compatible
# - $symbolCol, the column index of the SYMBOL column (in data lines)
# - ref to %genoCols whose keys are "HV","HET,"OTHER","HR" and values are the
#   corresponding columns (in data lines)
# - $tmpOutFilesR, ref to array of filehandles open for writing, one for each cohort,
#   same indexes as @cohorts
sub processBatch {
    (@_ == 9) || die "E $0: processBatch needs 9 args\n";
    my ($linesR,$knownCandidateGenesR,$sample2cohortR,$cohortsR,$sample2causalR,
        $compatibleR,$symbolCol,$genoColsR,$tmpOutFilesR) = @_;

    foreach my $line (@$linesR) {
        my @fields = split(/\t/, $line, -1) ;

        # $symbol doesn't depend on cohorts
        my $symbol = $fields[$symbolCol];

        # array of references (to sampleList arrays), one per cohort, same order
        # as in $cohortsR:
        # each element of @cohort2samplelists is a ref to a sampleLists array, ie
        # for each $cohort we build an array of 8 lists of samples (from the HV and
        # HET columns), in the following order and as defined at the top of this file:
        # $cohort_HV, $cohort_HET, $cohort_OTHERCAUSE_HV, $cohort_OTHERCAUSE_HET,
        # COMPAT_HV, COMPAT_HET, NEGCTRL_HV, NEGCTRL_HET
        # Each "list of samples" is actually a "genoData", ie it could look like:
        # 0/1~sample0112[39:0.95],sample0129[43:1.00]
        # If there is no sample in the category the array element remains empty.
        my @cohort2sampleLists;

        # similary for each cohort store a ref to an array of 8 COUNTs corresponding
        # to the sampleLists
        my @cohort2SLcounts;

        # initialize both arrays with refs to arrays holding empty strings / zeroes
        foreach my $cohorti (0..$#$cohortsR) {
            $cohort2sampleLists[$cohorti] = ["","","","","","","",""];
            $cohort2SLcounts[$cohorti] = [0,0,0,0,0,0,0,0];
        }
        
        # parse the HV and HET data fields and fill @cohort2samplelists,
        # @genoNames MUST BE in same order as in @cohort2sampleLists
        my @genoNames = ("HV","HET");
        foreach my $gni (0..$#genoNames) {
            my $genoData = $fields[$genoColsR->{$genoNames[$gni]}];
            # skip if no sample has this genotype
            ($genoData) || next;
            # sanity: at most one genotype (except for OTHER column)
            ($genoData =~ /\|/) &&
                die "E $0: more than one genotype for geno $genoNames[$gni], impossible. Line:\n$line\n";
            ($genoData =~ /^(\d+\/\d+)~([^~\|]+)$/) ||
                die "E $0: cannot parse GENOS data $genoData in line:\n$line\n";
            # $geno is the genotype (eg 1/1 or 0/2)
            my $geno = $1;
            my @samples = split(/,/,$2);
            foreach my $sample (@samples) {
                my $sampleID = $sample;
                # remove trailing [DP:AF] / [GQ:FR:BP] if it's there
                $sampleID =~ s/\[[^\]]+\]$//;
                # sanity check
                (defined $sample2cohortR->{$sampleID}) ||
                    die "E $0: sampleID $sampleID doesn't exist, sample was $sample\n";

                foreach my $cohorti (0..$#$cohortsR) {
                    my $cohort = $cohortsR->[$cohorti];
                    if ($sample2cohortR->{$sampleID} eq $cohort) {
                        # $sample belongs to $cohort
                        if ((! defined $sample2causalR->{$sampleID}) || ($sample2causalR->{$sampleID} eq $symbol)) {
                            # sample has no causal gene or it's the current gene: add to $cohort_HV or HET
                            if ($cohort2sampleLists[$cohorti]->[$gni]) {
                                $cohort2sampleLists[$cohorti]->[$gni] .= ",$sample";
                            }
                            else {
                                $cohort2sampleLists[$cohorti]->[$gni] = "$geno~$sample";
                            }
                            $cohort2SLcounts[$cohorti]->[$gni]++;
                        }
                        else {
                            # sample has a causal variant in another gene: OTHERCAUSE is at index 2+$gni
                            if ($cohort2sampleLists[$cohorti]->[2+$gni]) {
                                $cohort2sampleLists[$cohorti]->[2+$gni] .= ",$sample";
                            }
                            else {
                                $cohort2sampleLists[$cohorti]->[2+$gni] = "$geno~$sample";
                            }
                            $cohort2SLcounts[$cohorti]->[2+$gni]++;
                        }
                    }
                    elsif (defined $compatibleR->{$sample2cohortR->{$sampleID}}->{$cohort}) {
                        # $sample belongs to a cohort compatible with $cohort, store at index 4+$gni
                        if ($cohort2sampleLists[$cohorti]->[4+$gni]) {
                            $cohort2sampleLists[$cohorti]->[4+$gni] .= ",$sample";
                        }
                        else {
                            $cohort2sampleLists[$cohorti]->[4+$gni] = "$geno~$sample";
                        }
                        $cohort2SLcounts[$cohorti]->[4+$gni]++;
                    }
                    else {
                        # $sample belongs to another cohort, use as NEGCTRL ie index 6+$gni
                        if ($cohort2sampleLists[$cohorti]->[6+$gni]) {
                            $cohort2sampleLists[$cohorti]->[6+$gni] .= ",$sample";
                        }
                        else {
                            $cohort2sampleLists[$cohorti]->[6+$gni] = "$geno~$sample";
                        }
                        $cohort2SLcounts[$cohorti]->[6+$gni]++;
                    }
                }
            }
        }

        # OTHERGENO is the same for every $cohort, it's simply copied:
        my $otherGeno = $fields[$genoColsR->{"OTHER"}];
        # COUNT the otherGeno samples: set to zero if empty, otherwise 
        # we have |-separated lists of different genotypes and each
        # of these is a  ,-separated lists of samples... therefore every
        # sample is followed by , or | except the last sample
        my $otherGenoCount = 0;
        if ($otherGeno) {
            my @otherGenoCount = ($otherGeno =~ /[,|]/g);
            $otherGenoCount = 1 + scalar(@otherGenoCount);
        }

        # COUNT_HR is also the same for every cohort, use the same method
        my $hrCount = 0;
        if ($fields[$genoColsR->{"HR"}]) {
            my @hrCount = ($fields[$genoColsR->{"HR"}] =~ /[,|]/g);
            $hrCount = 1 + scalar(@hrCount);
        }

        # OK, now print stuff for every cohort that has at least one sample in
        # $cohort_HV or $cohort_HET or $cohort_OTHERCAUSE_*
        foreach my $cohorti (0..$#$cohortsR) {
            my $cohort = $cohortsR->[$cohorti];

            # skip if no HV or HET sample (possibly with other known causal variant)
            # in this cohort
            ($cohort2SLcounts[$cohorti]->[0] > 0) || ($cohort2SLcounts[$cohorti]->[1] > 0) ||
                ($cohort2SLcounts[$cohorti]->[2] > 0) || ($cohort2SLcounts[$cohorti]->[3] > 0) || next;

            # OK we have some data to print for $cohort
            # always keep first field 
            my $toPrint = "$fields[0]";
            foreach my $i (1..$#fields) {
                if ($i == $symbolCol) {
                    # prepend apostrophe-space to gene names to avoid excel corrupting everything
                    $toPrint .= "\t\' $fields[$i]\t";

                    # KNOWN_CAND column:
                    if ($knownCandidateGenesR->{$fields[$i]}) {
                        $toPrint .= $knownCandidateGenesR->{$fields[$i]};
                    }
                    # else: not a known candidate for any patho, leave empty
                    
                    # print all cohort-specific COUNTs
                    $toPrint .= "\t".join("\t",@{$cohort2SLcounts[$cohorti]});
                    # also print OTHERGENO and HR counts
                    $toPrint .= "\t$otherGenoCount\t$hrCount";
                }
                elsif ($i == $genoColsR->{"HV"}) {
                    # HV -> print all 8 sampleLists followed by OTHERGENO column
                    $toPrint .= "\t".join("\t",@{$cohort2sampleLists[$cohorti]});
                    $toPrint .= "\t$otherGeno";
                }
                elsif (($i == $genoColsR->{"HET"}) || ($i == $genoColsR->{"OTHER"}) ||
                       ($i == $genoColsR->{"HR"})) {
                    # NOOP: skip these columns
                }
                else {
                    # print other columns as-is
                    $toPrint .= "\t$fields[$i]";
                }
            }
            print { $tmpOutFilesR->[$cohorti] } "$toPrint\n";
        }
    }
}


###############
# this function waits for flagfiles to be created and "eats" the 
# corresponding tmpFiles in order, starting at 1.
# "eating" means print lines to the relevant $outFH and remove the tmpfiles.
# We also watch for $tmpOutLast, a file that will tell us the last batch number
# to wait for.
sub eatTmpFiles {
    (@_ == 3) || die "E $0: eatTmpFiles needs 3 args.\n";
    my ($tmpDir,$cohortsR,$outFHsR) = @_;

    # NOTE: all tmp filenames (eg $tmpOutLast) are hard-coded here and 
    # MUST MATCH those created by the worker threads.

    # when created, $tmpOutLast will contain the number of the last batch
    my $tmpOutLast = "$tmpDir/lastBatch";
    # $lastBatch remains undefined until we read it from $tmpOutLast
    my $lastBatch;

    # next batch to eat
    my $nextBatch = 1;

    while(1) {
        my $tmpOutFlag = "$tmpDir/$nextBatch.done";

        if (-e $tmpOutFlag) {
            # next batch of tmp files are ready
            foreach my $i (0..$#$cohortsR) {
                my $tmpFile = "$tmpDir/$nextBatch.$cohortsR->[$i].tsv";
                open (IN, $tmpFile) || 
                    die "E $0: in eatTmpFiles, flagfile $tmpOutFlag exists but cant read tmpFile $tmpFile: $!\n";
                while(<IN>) {
                    print {$outFHsR->[$i]} $_;
                }
                close(IN);
                (unlink($tmpFile) == 1) ||
                    die "E $0: in eatTmpFiles, done with tmpFile $tmpFile but cannot unlink it: $!\n";
            }
            (unlink($tmpOutFlag) == 1) ||
                die "E $0: in eatTmpFiles, done with files for batch $nextBatch but cannot unlink tmpOutFlag: $!\n";

            my $now = strftime("%F %T", localtime);
            # progress log: one INFO message every 10 batches
            ($nextBatch % 10) || warn("I $now: $0 - done processing batch $nextBatch\n");
            $nextBatch++;
            next;
        }

        elsif (-e $tmpOutLast) {
            open (IN, $tmpOutLast) || 
                die "E $0: in eatTmpFiles, cannot open tmpOutLast $tmpOutLast although it exists: $!\n";
            $lastBatch = <IN>;
            chomp($lastBatch);
            close(IN);
            unlink($tmpOutLast) || 
                die "E $0: in eatTmpFiles, cannot unlink tmpOutLast $tmpOutLast: $!\n";
            next;
        }

        elsif ((defined $lastBatch) && ($lastBatch < $nextBatch)) {
            # all done, return so this process can finish
            return();
        }

        else {
            # wait a few seconds before looking again
            sleep(10);
            next;
        }
    }
}

