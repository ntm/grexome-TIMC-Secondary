#!/usr/bin/perl

# NTM
# 09/07/2019 (as _strelka version, but starting from the older GATK version)


# Parses on stdin a Strelka or GATK GVCF file with one or more sample data columns;
# Args: see $USAGE.
# Prints to stdout a VCF (default) or GVCF (with --keepHR) file where:
# - samples that don't appear in $samplesFile "sampleID" column are removed
#   (allows to discard data for samples that are in the GVCF but were
#   obsoleted as dupes);
# - ignore all samples except the $samplesOfInterest, if specified;
# - QUAL is cleared to '.';
# - INFO is cleared to '.' (except for non-variant blocks with --keepHR, where 
#      INFO is kept as-is since END= is required)
# - phased genotypes x|y are replaced by unphased x/y;
# - hemizygous calls x (strelka) or x/* or */x (gatk) are replaced by HV x/x;
# - the new and useless GATK "weAreInAHomoDel" calls */* are replaced by ./.;
# - the variant calls in data columns are replaced by ./. if a call-condition
#   is not met (see "heuristics"), or if previous call was '.' (the Strelka NOCALL);
# - without --keepHR, lines where every sample is now ./. or 0/0 are skipped;
# - with --keepHR, lines where every sample is now ./. are skipped;
# - AF is moved (if it pre-exists) or added (otherwise) to FORMAT right after GT, 
#   and every 0/x or x/x call gets for AF the fraction of variant reads (rounded
#   to 2 decimals), HR and x/y calls get '.';
# - fix blatantly wrong genotype calls, see "heuristics" below.
# - ALTs that are not called in any sample (after fixing errors) are removed, and
#   DATA values are adjusted accordingly: GT is fixed (new ALT indexes), AD/ADF/ADR and
#   PL values for removed ALTs are discarded
# - work around strelka "feature": indel positions can be preceded by an HR call
#   that includes the first base of the indel (HR call at $pos or non-variant block
#   with END=$pos followed by indel call at the same $pos).


use strict;
use warnings;
use File::Basename qw(basename);
use FindBin qw($RealBin);
use Getopt::Long;
use POSIX qw(strftime);
use Parallel::ForkManager;

use lib "$RealBin";
use grexome_metaParse qw(parseSamples);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## hard-coded stuff that shouldn't change much

# max number of lines to read in a single batch. Each batch is then
# processed by a worker thread.
# Reduce if you are filling up $tmpDir (which should be on a RAMDISK),
# increase if jobs are almost instantaneous (because you are then 
# wasting time in parallelization overhead)
my $batchSize = 500000;


# heuristics for fixing low-quality or blatantly wrong genotype calls 
# [$dp below represents max(DP,sumOfADs), except for non-variant blocks where it is MIN_DP]:
# if $dp < $minDP , any call becomes NOCALL
# if max(GQ,GQX) < $minGQ , any call becomes NOCALL
# if AF < $minAF and call was REF/VAR or VAR/VAR, call becomes NOCALL
# if $dp >= $minDP_HV and AF >= $minAF_HV , call becomes HV
# if $dp >= $minDP_HET and AF >= $minAF_HET and AF <= $maxAF_HET, call becomes HET
my %filterParams = (
    "minDP" => 10,
    "minGQ" => 20,
    "minAF" => 0.15,
    "minDP_HV" => 20,
    "minAF_HV" => 0.85,
    "minDP_HET" => 20,
    "minAF_HET" => 0.25,
    "maxAF_HET" => 0.75);


#############################################
## options / params from the command-line

# samples metadata XLSX, no default
my $samplesFile;

# comma-separated list of samples of interest, if not empty all
# other samples are skipped. If empty every sample is kept.
my $samplesOfInterest;

# keepHR: if true keep lines even if they have only homoref calls 
# (ie produce GVCF), otherwise HR-only lines are skipped (produce VCF)
my $keepHR = '';

# for multi-threading, need to create a tmpDir. It will be removed
# when we are done and must not pre-exist (but it's parent must exist).
# To improve performance it should be on a ramdisk.
my $tmpDir = "tmpdir_filterBadCalls/";

# number of parallel jobs to run
my $numJobs = 16;

# if $verbose > 0 print more info to stderr (currently about 
# fixed genotype calls), larger values increase verbosity
my $verbose = 0;

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "Parse a Strelka or GATK4 GVCF on stdin, print to stdout a similar GVCF or VCF where:
- calls that fail basic quality filters are changed to NOCALL,
- calls that are blatantly wrong are fixed,
- ALTs that are never called are removed,
- lines are only printed if at least one sample still has some genotype call (including HomoRefs with --keepHR, excluding HomoRefs without).
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--samplesFile string [no default] : samples metadata xlsx file, with path
--samplesOfInterest string [default = all samples in xlsx] : comma-separated list of sampleIDs of interest (other samples are ignored)
--keepHR : keep lines even if the only genotype call is homoref
--tmpdir string [default = $tmpDir] : subdir where tmp files will be created (on a RAMDISK if possible), it must not pre-exist (but it's 
	 parent must exist) and it will be removed after execution
--jobs N [default = $numJobs] : number of parallel jobs=threads to run
--verbose N [default 0] : if > 0 increase verbosity on stderr
--help : print this USAGE";

# construct string with full command-line for adding to headers, must be
# done before GetOptions
my $addToHeader = "$0 ".join(" ",@ARGV);
chomp($addToHeader);
$addToHeader .= " > ".`readlink -f /proc/$$/fd/1` ;
chomp($addToHeader);
$addToHeader .= " 2> ".`readlink -f /proc/$$/fd/2` ;
chomp($addToHeader);
$addToHeader = "##filterBadCalls=<commandLine=\"$addToHeader\">\n";

GetOptions ("samplesFile=s" => \$samplesFile,
	    "samplesOfInterest=s" => \$samplesOfInterest,
	    "keepHR" => \$keepHR,
	    "jobs=i" => \$numJobs,
	    "verbose=i" => \$verbose,
	    "tmpdir=s" => \$tmpDir,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($samplesFile) || die "E $0: you must provide a samplesFile file\n";
(-f $samplesFile) || die "E $0: the supplied samplesFile file doesn't exist\n";

(-e $tmpDir) && 
    die "E $0: tmpDir $tmpDir exists, please remove or rename it, or provide a different one with --tmpdir\n";
mkdir($tmpDir) || 
    die "E $0: cannot mkdir tmpDir $tmpDir\n";

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


#########################################################
# parse samples metadata file to grab sampleIDs, limit to samples of interest

# key==existing sample of interest, value==1
my %samples = ();

# just use the first hashref from parseSamples, ignoring pathologyIDs
my ($s2pathoR) = &parseSamples($samplesFile);
foreach my $s (keys %$s2pathoR) {
    $samples{$s} = 1;
}

if ($samplesOfInterest) {
    # make sure every listed sample is in %samples and promote it's value to 2
    foreach my $soi (split(/,/, $samplesOfInterest)) {
	($samples{$soi}) ||
	    die "E $0: processing samplesOfInterest: a specified sample $soi does not exist in the samples metadata file\n";
	($samples{$soi} == 1) ||
	    warn "W $0: processing samplesOfInterest: sample $soi was specified twice, is that a typo?\n";
	$samples{$soi} = 2;
    }
    # now demote the SOIs to 1 and ignore all other samples
    foreach my $s (keys %samples) {
	if ($samples{$s} == 2) {
	    $samples{$s} = 1;
	}
	else {
	    delete($samples{$s});
	}
    }
}
    
#############################################
# deal with headers

# array, same number of elements as there are columns in the #CHROM line
# (and hence in each data line), value is true iff column must be skipped
# (corresponding to samples that no longer exist or are obsoleted in samplesFile, 
# eg they were dupes of other samples with better sequencing)
my @skippedCols = ();

# accumulator for header, so we can avoid printing anything if stdin doesn't 
# contain any samples (of interest)
my $headerToPrint = "";

# parse header, just copy it except #CHROM
while(my $line = <STDIN>) {
    if ($line =~ /^##/) {
	$headerToPrint .= $line;
    }
    elsif ($line =~ /^#CHROM/) {
	# add ##comment with full command line run
	$headerToPrint .= $addToHeader;
	# remove samples that don't exist in $samplesFile (anymore) or are not of interest
	# also remember if at least one sample remains
	my $goodSamples = 0;
	chomp($line);
	my @fields = split(/\t/,$line);
	foreach my $i (reverse(9..$#fields)) {
	    # reverse so we can splice bad columns out
	    if (! $samples{$fields[$i]}) {
		splice(@fields,$i,1);
		$skippedCols[$i] = 1;
	    }
	    else {
		$goodSamples = 1;
	    }
	}
	if (! $goodSamples) {
	    # STDIN doesn't contain any samples of interest
	    warn "W $0: no valid samples (of interest) here, nothing to do\n";
	    rmdir($tmpDir) || 
		die "E $0: no samples of interest but cannot rmdir tmpDir $tmpDir, why? $!\n";
	    exit(0);
	}
	else {
	    $headerToPrint .= join("\t",@fields)."\n";
	    print $headerToPrint;
	    last;
	}
    }
    else {
	die "E $0: parsing header, found bad line:\n$line";
    }
}


# flush stdout before starting our eatTmpFiles job
STDOUT->flush();


#############################################
# read data lines

# create fork manager
my $pm = new Parallel::ForkManager($numJobs);

# spawn a child process that waits for workers to finish producing batches,
# and prints the tmpfiles to stdout in correct order, cleaning up behind itself
if (! $pm->start) {
    &eatTmpFiles($tmpDir);
    $pm->finish;
}

# number of the current batch
my $batchNum = 0;

# for making sure a batch doesn't start with an indel, stores the chomped
# first line to place in the next batch (or '' if there are no more lines)
my $lineForNextBatch = <STDIN>;
chomp($lineForNextBatch);

while ($lineForNextBatch) {
    $batchNum++;
    my @lines = ($lineForNextBatch);
    $lineForNextBatch = '';
    my $eaten = 0;
    while (my $line = <STDIN>) {
	chomp($line);
	if (($eaten >= $batchSize) && ($line =~ /^[^\t]+\t[^\t]+\t[^\t]+\t.\t([^\t]+)\t/)) {
	    # batch already has enough lines and has a single-char REF...
	    my $alts = $1;
	    if (($alts eq '<NON_REF>') || ($alts =~ /^.,<NON_REF>$/) || ($alts =~ /^.$/)) {
		# non-variant position or single-char ALT, can't be an indel
		$lineForNextBatch = $line;
		last;
	    }
	}
	# if we didn't "last" batch isn't full yet, or $line might be an indel -> in both cases eat it
	push(@lines,$line);
	$eaten++;
    }

    # let worker threads take it from there
    $pm->start && next;
    # NOTE: IF YOU CHANGE the tmp filenames below ($tmpOut, $tmpOutFlag),
    # you MUST EDIT &eatTmpFiles()

    # create tmp output filehandle for this batch
    my $tmpOut = "$tmpDir/$batchNum.vcf";
    open(my $tmpOutFH, "> $tmpOut") || die "E $0: cannot open $tmpOut for writing\n";

    # process this batch
    &processBatch(\@lines,$tmpOutFH,\%filterParams,\@skippedCols,$keepHR,$verbose);

    # done, close tmp FH and create flag-file
    close($tmpOutFH) || die "E $0: cannot close tmp outFH $tmpOutFH\n";
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

$now = strftime("%F %T", localtime);
rmdir($tmpDir) || 
    die "E $now: $0 - all done but cannot rmdir tmpDir $tmpDir, why? $!\n";
warn "I $now: $0 - ALL DONE, completed successfully!\n";


#############################################
## subs


# process a batch of lines
# args:
# - ref to array of chomped lines
# - outFH open filehandle to print to
# - hashref with filter params
# - ref to array saying which columns to skip
# - $keepHR: HR-only lines are skipped if false, kept if true
# - $verbose, 0 is quiet, increase value for more verbosity
# WARNING: if a batch ends with a non-variant call or block and the
# next batch starts with an indel, the prevToPrint mechanism to work
# around a strelka bug won't work at the batch boundary. To avoid this
# make sure batches don't start with an indel position.
sub processBatch {
    (@_ == 6) || die "E $0: processBatch needs 6 args\n";
    my ($linesR,$outFH,$filterParamsR,$skippedColsR,$keepHR,$verbose) = @_;

    # counters for number of blatant errors fixed to HV or HET, and for fixed DPs
    my $fixedToHV = 0;
    my $fixedToHET = 0;
    my $fixedDP = 0;

    # delay printing lines so we can decrement END= if needed
    # (to work-around strelka bug: indels can be preceded by HR calls
    # at the same POS or by non-variant blocks whose END= goes one too far)
    # => our solution: delete the HR call at the same POS as a subsequent indel
    my @prevToPrint = ();
    
    foreach my $line (@$linesR) {
	# $keepLine: boolean, true if at least one non-'*' ALT is called for at 
	# least one sample after cleaning and filtering
	my $keepLine = 0;
	my @data = split(/\t/, $line);
	(@data >= 10) || die "E $0: no sample data in line?\n$line\n";

	# BEFORE ANYTHING ELSE: deal with @prevToPrint
	if (@prevToPrint) {
	    if (($prevToPrint[0] eq $data[0]) && ($prevToPrint[1] == $data[1]) && 
		(($prevToPrint[4] eq '.') || ($prevToPrint[4] eq '<NON_REF>'))) {
		# prevToPrint was HR line at same POS as $line, don't print prev -> NOOP
	    }
	    else {
		# if prev was a non-var block on same chrom ending at current POS, decrement END=
		my $thisPos = $data[1];
		my $prevEnd = $thisPos - 1;
		$prevToPrint[7] =~ s/END=$thisPos;/END=$prevEnd;/;
		# whether we substituted or not, print prev
		print $outFH join("\t", @prevToPrint)."\n";
	    }
	    # in all cases clear prev
	    @prevToPrint = ();
	}
	
	# if not --keepHR and there is no ALT in line, skip immediately
	(! $keepHR) && (($data[4] eq '.') || ($data[4] eq '<NON_REF>')) && next;
	# GATK4 produces useless lines where there are NO sequencing reads
	# (where FORMAT is eg GT or GT:GQ:PL), skip them immediately:
	# any line with supporting reads must have a DP or AD field
	($data[8] =~ /:DP:/) || ($data[8] =~ /:AD:/) ||
	    ($data[8] =~ /:DP$/) || ($data[8] =~ /:AD$/) || next;

	# after parsing the line, $altsCalled[$i] will be true iff $alts[$i] is part of
	# a called geno for at least one sample
	my @alts = split(/,/,$data[4]);
	my @altsCalled;
	# grab alleleNum of ALT '*' if it's present
	my $starNum = -1;
	foreach my $alti (0..$#alts) {
	    ($alts[$alti] eq '*') && ($starNum = $alti + 1) && last;
	}
	# first 9 fields are copied as-is except: QUAL is cleared, INFO is cleared except 
	# if it contains END=, and AF is moved or added to FORMAT after GT
	# NOTE: ALT may be modified later (remove uncalled ALTs and normalize remaining ALTs)
	my @lineToPrint = @data[0..4];
	# clear QUAL, copy FILTER
	push(@lineToPrint, '.', $data[6]);
	# copy INFO if it contains END=, clear otherwise
	if ($data[7] =~ /^END=/) {
	    push(@lineToPrint, $data[7]);
	}
	else {
	    push(@lineToPrint, '.');
	}
	my $format = $data[8];
	my $newFormat = $format;
	# if AF isn't already right after GT we move it / create it there
	if ($newFormat !~ /^GT:AF:/) {
	    # remove it if it existed
	    $newFormat =~ s/:AF:/:/;
	    # in any case add it after GT
	    ($newFormat =~ s/^GT:/GT:AF:/)  || 
		die "E $0: cannot add AF after GT in format: $format\n";
	}
	push(@lineToPrint, $newFormat);
	# %format: key is a FORMAT key (eg DP), value is the index of that key in $format
	my %format;
	{
	    my @format = split(/:/, $format);
	    foreach my $i (0..$#format) {
		$format{$format[$i]} = $i ;
	    }
	}
	# sanity: make sure the fields we need are there
	(defined $format{"GT"}) || die "E $0: no GT key in FORMAT string for line:\n$line\n";
	(defined $format{"GQ"}) || (defined $format{"GQX"}) ||die "E $0: no GQ or GQX key in FORMAT string for line:\n$line\n";
	(defined $format{"AD"}) || (defined $format{"DP"}) || die "E $0: no AD or DP key in FORMAT string for line:\n$line\n";

	# now deal with actual data fields
	foreach my $i (9..$#data) {
	    ($skippedColsR->[$i]) && next;
	    my $thisData = $data[$i];
	    # if genotype is already '.' or './.' == NOCALL, just use ./.
	    if (($thisData =~ m~^\.$~) || ($thisData =~ m~^\.:~) || ($thisData =~ m~^\./\.~)) {
		push(@lineToPrint, './.') ;
		next;
	    }
	    # also if call is */* or *|* , just replace with ./.
	    if ($thisData =~ m~^$starNum[/|]$starNum:~) {
		push(@lineToPrint, './.') ;
		next;
	    }

	    # otherwise examine content and apply filters
	    my @thisData = split(/:/, $thisData) ;

	    # calculate $gq = max(GQ,GQX) making sure things are defined.
	    # $gq stays at -1 if GQ and GQX are both undef or '.'
	    my $gq = -1;
	    if ((defined $format{"GQ"}) && (defined $thisData[$format{"GQ"}]) && ($thisData[$format{"GQ"}] ne '.')) {
		$gq = $thisData[$format{"GQ"}];
	    }
	    if ((defined $format{"GQX"}) && (defined $thisData[$format{"GQX"}]) &&
		($thisData[$format{"GQX"}] ne '.') && ($thisData[$format{"GQX"}] > $gq)) {
		$gq = $thisData[$format{"GQX"}];
	    }
	    if ($gq < $filterParamsR->{"minGQ"}) {
		# GQ and GQX (if it exists) are both undef or too low, change to NOCALL
		push(@lineToPrint, './.') ;
		next;
	    }

	    # clean up GTs -> obtain unphased, biallelic, sorted genotype
	    $thisData[$format{"GT"}] = &cleanGT($thisData[$format{"GT"}], $starNum);

	    # if '*' is in ALTs we want to ignore any reads attributed to it in DP and AD,
	    # so our DP and AF are correct (and to avoid incorrectly "fixing" calls)
	    # NOTE: we don't touch anything else (eg GQ, PL, SB, and GATK can't output ADF/ADR currently)
	    if (($starNum != -1) && (defined $format{"AD"}) && (defined $thisData[$format{"AD"}]) &&
		($thisData[$format{"AD"}]  =~ /^[\d,]+$/)) {
		# '*' is in ALTs and AD is defined and has data
		my @ADs = split(/,/,$thisData[$format{"AD"}]);
		# sanity
		(@ADs == 1 + @alts) ||
		    die "E $0: WTF! AD has wrong number of values in $thisData from line:\n$line\n";
		# sanity: '*' means this is GATK and if we have AD in GATK we must have DP
		((defined $format{"DP"}) && (defined $thisData[$format{"DP"}]) && ($thisData[$format{"DP"}] =~ /^\d+$/)) ||
		    die "E $0: WTF! we have AD but not DP in $thisData from line:\n$line\n";
		# decrement DP as needed and set AD for '*' to 0
		$thisData[$format{"DP"}] -= $ADs[$starNum];
		$ADs[$starNum] = 0;
		$thisData[$format{"AD"}] = join(',', @ADs);
	    }

	    # grab the depth (DP or sumOfADs, whichever is defined and higher, except use MIN_DP if in a non-var block)
	    my $thisDP = -1;
	    if ((defined $format{"DP"}) && (defined $thisData[$format{"DP"}]) && ($thisData[$format{"DP"}] ne '.')) {
		$thisDP = $thisData[$format{"DP"}];
	    }
	    if ((defined $format{"AD"}) && (defined $thisData[$format{"AD"}]) && ($thisData[$format{"AD"}] =~ /^[\d,]+$/)) {
		my $sumOfADs = 0;
		foreach my $ad (split(/,/,$thisData[$format{"AD"}])) {
		    $sumOfADs += $ad;
		}
		if ($thisDP < $sumOfADs) {
		    # if DP exists there's no sensible reason for it to be smaller than sumOfADs
		    # (but we have seen it happen with GATK), fix it
		    if ($thisDP > -1) {
			$thisData[$format{"DP"}] = $sumOfADs;
			$fixedDP++;
			if ($verbose >= 2) {
			    warn "I $0: fix DP in: $data[0]:$data[1] $data[3] > $data[4] sample ".($i-9)." DP=$thisDP sumOfADs=$sumOfADs\n";
			}
		    }
		    $thisDP = $sumOfADs;
		}
	    }
	    if ((defined $format{"MIN_DP"}) && (defined $thisData[$format{"MIN_DP"}]) && ($thisData[$format{"MIN_DP"}] ne '.')) {
		$thisDP = $thisData[$format{"MIN_DP"}];
	    }
	    # if depth too low or undefined for this sample, change to NOCALL
	    if ($thisDP < $filterParamsR->{"minDP"}) {
		push(@lineToPrint, './.') ;
		next;
	    }

	    # calculate AF, for fracVarReads filter
	    # NOTE: with GATK I had some issues with DP=0 calls, causing illegal divisions
	    # by zero when calculating AF, but that is now skipped above
	    my $af;
	    my ($geno1,$geno2) = split(/\//, $thisData[$format{"GT"}]);
	    if ((defined $format{"AF"}) && (defined $thisData[$format{"AF"}])) {
		# AF was already there, just reuse
		$af = $thisData[$format{"AF"}];
	    }
	    elsif (($geno2 != 0) && (($geno1 == 0) || ($geno1 == $geno2))) {
		# 0/x HET or x/x HV, AD should always be there
		if ((!defined $thisData[$format{"AD"}]) || ($thisData[$format{"AD"}] !~ /^[\d,]+$/)) {
		    die "E $0: GT is HET or HV but we don't have AD or AD data is blank in:\n$line\n";
		}
		my @ads = split(/,/, $thisData[$format{"AD"}]);
		# $geno2 is always the index of the VAR (thanks to sorting in cleanGT)
		my $fracVarReads = $ads[$geno2] / $thisDP ;
		# round AF to nearest float with 2 decimals
		$af = sprintf("%.2f",$fracVarReads);
	    }
	    else {
		# AF doesn't pre-exist and this is HR or x/y, set AF='.'
		$af = '.';
	    }

	    if (($af ne '.') && ($af < $filterParamsR->{"minAF"})) {
		# AF too low, change to NOCALL
		push(@lineToPrint, './.') ;
		next;
	    }

	    # we have $thisDP and $af , fix blatantly wrong calls
	    if (($thisDP >= $filterParamsR->{"minDP_HV"}) && ($geno1 == 0) &&
		($af ne '.') && ($af >= $filterParamsR->{"minAF_HV"})) {
		# change to HV
		$thisData[$format{"GT"}] = "$geno2/$geno2";
		$fixedToHV++;
		if ($verbose >= 2) {
		    # warn with chrom pos ref > alts sample dp af
		    warn "I $0: fix to HV, $data[0]:$data[1] $data[3] > $data[4] sample ".($i-9)." DP=$thisDP AF=$af\n";
		}
	    }
	    if (($thisDP >= $filterParamsR->{"minDP_HET"}) && ($geno1 != 0) && ($af ne '.') && 
		($af >= $filterParamsR->{"minAF_HET"}) && ($af <= $filterParamsR->{"maxAF_HET"})) {
		# change to HET
		$thisData[$format{"GT"}] = "0/$geno2";
		$fixedToHET++;
		if ($verbose >= 2) {
		    # warn with chrom pos ref > alts sample dp af
		    warn "I $0: fix to HET, $data[0]:$data[1] $data[3] > $data[4] sample ".($i-9)." DP=$thisDP AF=$af\n";
		}
	    }

	    # other filters (eg strandDisc) could go here

	    # OK data passed all filters, AF may need to be moved or added
	    if ((defined $format{"AF"}) && (defined $thisData[$format{"AF"}])) {
		if ($format{"AF"} != 1) {
		    # remove from previous position, and add back in second position
		    splice(@thisData, $format{"AF"}, 1);
		    splice(@thisData, 1, 0, $af);
		}
		# else AF was already in second position, noop
	    }
	    else {
		# AF wasn't present, add it in second position
		splice(@thisData, 1, 0, $af);
	    }

	    push(@lineToPrint, join(':',@thisData));

	    if ($keepHR || ($thisData[$format{"GT"}] ne '0/0')) {
		$keepLine = 1;
	    }

	    # remember any called ALT allele
	    ($geno1 > 0) && ($altsCalled[$geno1 - 1] = 1);
	    ($geno2 > 0) && ($altsCalled[$geno2 - 1] = 1);
	}
	
	# done parsing $line
	if ($keepLine) {
	    # remove any uncalled ALTs and fix DATA
	    &removeUncalledALTs(\@lineToPrint, \@altsCalled);
	    # save for printing
	    @prevToPrint = @lineToPrint;
	}
	# else this line is discarded -> NOOP
    }
    # print last line if needed
    (@prevToPrint) && (print $outFH join("\t", @prevToPrint)."\n");
    # INFO with number of fixed calls in this batch, we don't care that this comes
    # out of order to stderr but don't log if we already printed each fixed call
    if ($verbose == 1) {
	($fixedToHV) && (warn "I $0: fixed $fixedToHV calls from HET to HV\n");
	($fixedToHET) && (warn "I $0: fixed $fixedToHET calls from HV to HET\n");
 	($fixedDP) && (warn "I $0: fixed $fixedDP DP values to sumOfADs (was larger)\n");
   }
}

###############
# clean up GTs:
# args: a GT string, the alleleNum of the '*' allele (or -1 if it's not present)
# returns a "cleaned up" version of the GT: unphased, biallelic, sorted
sub cleanGT {
    (@_ == 2) ||die "E $0: cleanGT needs 2 args.\n";
    my ($gt, $starNum) = @_;
    
    # Strelka and GATK make some phased calls sometimes, homogenize as unphased
    $gt =~ s~\|~/~ ;
    # Strelka makes some hemizygous calls as 'x' (eg when the position
    # is in a HET deletion), makes sense but still, homogenize as HOMO
    $gt =~ s~^(\d+)$~$1/$1~;
    # grab geno
    my ($geno1,$geno2) = split(/\//, $gt);
    ((defined $geno1) && (defined $geno2)) ||
	die "E $0: a sample's genotype cannot be split: $gt\n";
    # GATK hemizygous calls (eg under a HET DEL) appear as x/* or */x, fix to x/x
    if ($geno2 == $starNum) {
	# */* shouldn't exist (replaced by ./. and skipped before calling &cleanGT)
	($geno1 == $starNum) &&
	    die "E $0: don't call this sub with genotype */* , just skip these useless calls\n"; 
	$geno2 = $geno1;
    }
    elsif ($geno1 == $starNum) {
	$geno1 = $geno2;
    }
    # make sure alleles are in sorted order
    if ($geno2 < $geno1) {
	my $genot = $geno1;
	$geno1 = $geno2;
	$geno2 = $genot;
    }
    # OK, return clean GT
    return("$geno1/$geno2");
}


###############
# removeUncalledALTs, args:
# - arrayref, tab-splitted VCF line
# - arrayref of ALT indexes that were called in at least one sample
#
# -> discard any uncalled ALTs from the ALT column, set to '.' if there are no more ALTs
# -> adjust every GT (new ALT indexes)
# -> discard AD/ADF/ADR and PL values for removed ALTs
sub removeUncalledALTs {
    (@_ == 2) ||die "E $0: removeUncalledALTs needs 2 args.\n";
    my ($lineToPrintR, $altsCalledR) = @_;


}


###############
# this function waits for flagfiles to be created and "eats" the 
# corresponding tmpFiles in order, starting at 1.
# "eating" means print lines to stdout and remove the tmpfiles.
# we also watch for $tmpOutLast, a file that will tell us
# the last batch number to wait for
sub eatTmpFiles {
    (@_ == 1) || die "E $0: eatTmpFiles needs 1 arg.\n";
    my ($tmpDir) = @_;

    # NOTE: all tmp filenames (eg $tmpOut) are hard-coded here and 
    # MUST MATCH those created by the main thread.

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
	    my $tmpOut = "$tmpDir/$nextBatch.vcf";
		open (IN, $tmpOut) || 
		    die "E $0: in eatTmpFiles, flagfile $tmpOutFlag exists but cant read tmpFile $tmpOut: $!\n";
	    while(<IN>) {
		print $_;
	    }
	    close(IN);
	    (unlink($tmpOut,$tmpOutFlag) == 2) ||
		die "E $0: in eatTmpFiles, done with tmpFile $tmpOut and flagfile $tmpOutFlag but cannot unlink them: $!\n";
	    # progress log: one INFO message every 10 batches
	    if (! $nextBatch % 10) {
		my $now = strftime("%F %T", localtime);
		warn "I $now: $0 - done processing batch $nextBatch\n";
	    }
	    $nextBatch++;
	    next;
	}

	elsif ((-e $tmpOutLast) && (! -z $tmpOutLast)) {
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
	    sleep(2);
	    next;
	}
    }
}

