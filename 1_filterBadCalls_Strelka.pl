#!/usr/bin/perl

# NTM
# 09/07/2019 (but starting from 1_filterBadCalls.pl which is older)


# Parses on stdin a Strelka GVCF file with one or more sample data columns;
# Args: see $USAGE.
# Prints to stdout a VCF file where:
# - samples that don't appear in $metadata "sampleID" column are removed
#   (allows to discard data for samples that are in the GVCF but were
#   obsoleted as dupes);
# - ignore all samples except the $samplesOfInterest, if specified;
# - non-variant lines are removed;
# - the variant calls in data columns are replaced by ./. if a call-condition
#   is not met (see "heuristics"), or if previous call was '.' (the Strelka NOCALL);
# - lines where every sample is now ./. or 0/0 are skipped;
# - AF is added to FORMAT right after GT, and every 0/x or x/x call gets
#   for AF the fraction of variant reads (rounded to 2 decimals), HR
#   and x/y calls get '.';
# - fix blatantly wrong genotype calls, see "heuristics" below.

use strict;
use warnings;
use Getopt::Long;
use Spreadsheet::XLSX;
use POSIX qw(strftime);
use Parallel::ForkManager;


#############################################
## hard-coded stuff that shouldn't change much

# max number of lines to read in a single batch. Each batch is then
# processed by a worker thread.
# Reduce if you are filling up $tmpDir (which should be on a RAMDISK),
# increase if jobs are almost instantaneous (because you are then 
# wasting time in parallelization overhead)
my $batchSize = 500000;


# heuristics for fixing low-quality or blatantly wrong genotype calls 
# [$dp below represents max(DP,sumOfADs)]:
# if $dp < $minDP , any call becomes NOCALL
# if GQX < $minGQX , any call becomes NOCALL
# if AF < $minAF and call was REF/VAR or VAR/VAR, call becomes NOCALL
# if $dp >= $minDP_HV and AF >= $minAF_HV , call becomes HV
# if $dp >= $minDP_HET and AF >= $minAF_HET and AF <= $maxAF_HET, call becomes HET
my %filterParams = (
    "minDP" => 10,
    "minGQX" => 20,
    "minAF" => 0.15,
    "minDP_HV" => 20,
    "minAF_HV" => 0.85,
    "minDP_HET" => 20,
    "minAF_HET" => 0.25,
    "maxAF_HET" => 0.75);


#############################################
## options / params from the command-line

# metadata XLSX, no default
my $metadata;

# comma-separated list of samples of interest, if not empty all
# other samples are skipped. If empty every sample is kept.
my $samplesOfInterest;

# for multi-threading, need to create a tmpDir. It will
# be removed when we are done and must not pre-exist.
# To improve performance it should be on a ramdisk.
my $tmpDir = "tmpdir_filterBadCalls/";

# number of parallel jobs to run
my $numJobs = 16;

# if $verbose > 0 print more info to stderr (currently about 
# fixed genotype calls), larger values increase verbosity
my $verbose = 0;

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "Parse a Strelka GVCF on stdin, print to stdout a similar GVCF where:
- calls that fail basic quality filters are changed to NOCALL,
- calls that are blatantly wrong are fixed,
- lines are only printed if at least one sample has a non-HR genotype call.\n
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--metadata string [no default] : patient metadata xlsx file, with path
--samplesOfInterest string [default = all samples in metadata xlsx] : comma-separated list of sampleIDs of interest (other samples are ignored)
--tmpdir string [default = $tmpDir] : subdir where tmp files will be created (on a RAMDISK if possible), must not pre-exist and will be removed after execution
--jobs N [default = $numJobs] : number of parallel jobs=threads to run
--verbose N [default 0] : if > 0 increase verbosity on stderr
--help : print this USAGE";

GetOptions ("metadata=s" => \$metadata,
	    "samplesOfInterest=s" => \$samplesOfInterest,
	    "jobs=i" => \$numJobs,
	    "verbose=i" => \$verbose,
	    "tmpdir=s" => \$tmpDir,
	    "help" => \$help)
    or die("Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) &&
    die "$USAGE\n\n";

($metadata) || die "E: you must provide a metadata file\n";
(-f $metadata) || die "E: the supplied metadata file doesn't exist\n";

(-e $tmpDir) && 
    die "E: tmpDir $tmpDir exists, please remove or rename it, or provide a different one with --tmpdir\n";
mkdir($tmpDir) || 
    die "E: cannot mkdir tmpDir $tmpDir\n";

    
#########################################################
# parse patient metadata file to grab sampleIDs, limit to samples of interest

# key==existing sample of interest, value==1
my %samples = ();

{
    my $workbook = Spreadsheet::XLSX->new("$metadata");
    (defined $workbook) ||
	die "E when parsing xlsx\n";
    ($workbook->worksheet_count() == 1) ||
	die "E parsing xlsx: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($sampleCol) = (-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	# if column has no header just ignore it
	(defined $cell) || next;
	($cell->value() eq "sampleID") &&
	    ($sampleCol = $col);
    }
    ($sampleCol >= 0) ||
	die "E parsing xlsx: no column title is sampleID\n";

    foreach my $row ($rowMin+1..$rowMax) {
	my $sample = $worksheet->get_cell($row, $sampleCol)->value;
	# skip "none" lines
	($sample eq "none") && next;
	(defined $samples{$sample}) && 
	    die "E parsing xlsx: have 2 lines with sample $sample\n";
	$samples{$sample} = 1;
    }
}

if ($samplesOfInterest) {
    # make sure every listed sample is in %samples and promote it's value to 2
    foreach my $soi (split(/,/, $samplesOfInterest)) {
	($samples{$soi}) ||
	    die "E processing samplesOfInterest: a specified sample $soi does not exist in the metadata file\n";
	($samples{$soi} == 1) ||
	    warn "W processing samplesOfInterest: sample $soi was specified twice, is that a typo?\n";
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

my $now = strftime("%F %T", localtime);
warn "I: $now - starting to run: ".join(" ", $0, @ARGV)."\n";

# array, same number of elements as there are columns in the #CHROM line
# (and hence in each data line), value is true iff column must be skipped
# (corresponding to samples that no longer exist in the metadata file, 
# eg they were dupes of other samples with better sequencing)
my @skippedCols = ();

# parse header, just copy it except we remove samples that don't exist anymore
while(my $line = <STDIN>) {
    if ($line =~ /^##/) {
	print $line;
    }
    elsif ($line =~ /^#CHROM/) {
	# add ##comment with full command line run
	my $com = qx/ps -o args= $$/;
	chomp($com);
	$com .= " < ".`readlink -f /proc/$$/fd/0` ;
	chomp($com);
	$com .= " > ".`readlink -f /proc/$$/fd/1` ;
	chomp($com);
	$com .= " 2> ".`readlink -f /proc/$$/fd/2` ;
	chomp($com);
	print "##filterBadCalls=<commandLine=\"$com\">\n";

	# remove samples that don't exist in $metadata (anymore) or are not of interest
	chomp($line);
	my @fields = split(/\t/,$line);
	foreach my $i (reverse(9..$#fields)) {
	    # reverse so we can splice bad columns out
	    if (! $samples{$fields[$i]}) {
		splice(@fields,$i,1);
		$skippedCols[$i] = 1;
	    }
	}
	print join("\t",@fields)."\n";
	last;
    }
    else {
	die "E: parsing header, found bad line:\n$line";
    }
}

# flush stdout before starting our eatTmpFiles job
STDOUT->flush();


#############################################
# parse data lines

# create fork manager
my $pm = new Parallel::ForkManager($numJobs);

# spawn a child process that waits for workers to finish producing batches,
# and prints the tmpfiles to stdout in correct order, cleaning up behind itself
if (! $pm->start) {
    &eatTmpFiles($tmpDir);
    $pm->finish;
}

# boolean flag, true iff current batch is the last
my $lastBatch = 0;
# number of the current batch
my $batchNum = 0;

while (!$lastBatch) {
    $batchNum++;
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

    # create tmp output filehandle for this batch
    my $tmpOut = "$tmpDir/$batchNum.vcf";
    open(my $tmpOutFH, "> $tmpOut") || die "cannot open $tmpOut for writing\n";

    # process this batch
    &processBatch(\@lines,$tmpOutFH,\%filterParams,\@skippedCols,$verbose);

    # done, close tmp FH and create flag-file
    close($tmpOutFH) || die "cannot close tmp outFH $tmpOutFH\n";
    my $tmpOutFlag = "$tmpDir/$batchNum.done";
    open(OUTFLAG, "> $tmpOutFlag") || die "cannot open flagfile $tmpOutFlag for writing\n";
    print OUTFLAG "$batchNum\n";
    close(OUTFLAG);
    $pm->finish;
}

# some children are still processing batches, but we know the last
# batchNum that will ever exist, tell &eatTmpFiles() so it can exit
# (of course if you change $tmpOutLast you have to edit &eatTmpFiles)
my $tmpOutLast = "$tmpDir/lastBatch";
open(OUTLAST, "> $tmpOutLast") || die "cannot open tmp-last-file $tmpOutLast for writing\n";
print OUTLAST "$batchNum\n";
close OUTLAST;

$pm->wait_all_children;

$now = strftime("%F %T", localtime);
rmdir($tmpDir) || 
    die "E: $now - all done but cannot rmdir tmpDir $tmpDir, why? $!\n";
warn "I: $now - DONE running: ".join(" ", $0, @ARGV)."\n";



#############################################
## subs


# process a batch of lines
# args:
# - ref to array of chomped lines
# - outFH open filehandle to print to
# - hashref with filter params
# - ref to array saying which columns to skip
# - $verbose, 0 is quiet, increase value for more verbosity
sub processBatch {
    (@_ == 5) || die "E: processBatch needs 5 args\n";
    my ($linesR,$outFH,$filterParamsR,$skippedColsR,$verbose) = @_;

    # counters for number of blatant errors fixed to HV or HET
    my $fixedToHV = 0;
    my $fixedToHET = 0;

    foreach my $line (@$linesR) {
	# $keepLine: boolean, true if at least one sample is not ./. and 0/0 after filtering
	my $keepLine = 0;
	my @data = split(/\t/, $line);
	(@data >= 10) || die "no sample data in line?\n$line\n";
	# if no ALT in line, skip immediately
	($data[4] eq '.') && next;
	# first 9 fields are copied except AF is added to FORMAT after GT
	my $lineToPrint = join("\t",@data[0..7]);
	my $format = $data[8];
	my $newFormat = $format;
	($newFormat =~ s/^GT:/GT:AF:/)  || 
	    die "E: cannot add AF after GT in format: $format\n";
	$lineToPrint .= "\t$newFormat";
	# %format: key is a FORMAT key (eg GQX), value is the index of that key in $format
	my %format;
	{
	    my @format = split(/:/, $format);
	    foreach my $i (0..$#format) {
		$format{$format[$i]} = $i ;
	    }
	}
	# sanity: make sure the fields we need are there
	# don't check DP since AD is checked
	(defined $format{"GQX"}) || die "no GQX key in FORMAT string for line:\n$line\n";
	(defined $format{"GT"}) || die "no GT key in FORMAT string for line:\n$line\n";
	(defined $format{"AD"}) || die "no AD key in FORMAT string for line:\n$line\n";

	# now deal with actual data fields
	foreach my $i (9..$#data) {
	    ($skippedColsR->[$i]) && next;
	    my $data = $data[$i];
	    # if genotype is already '.' or './.' == NOCALL, just use ./.
	    if (($data =~ m~^\.$~) || ($data =~ m~^\.:~) || ($data =~ m~^\./\.~)) {
		$lineToPrint .= "\t./." ;
		next;
	    }
	    # otherwise examine content and apply filters
	    my @thisData = split(/:/, $data) ;

	    if ((! $thisData[$format{"GQX"}]) || ($thisData[$format{"GQX"}] eq '.') ||
		($thisData[$format{"GQX"}] < $filterParamsR->{"minGQX"})) {
		# GQX undefined or too low, change to NOCALL
		$lineToPrint .= "\t./.";
		next;
	    }

	    # grab the depth (DP or sumOfADs, whichever is defined and higher)
	    my $thisDP = -1;
	    if ((defined $format{"DP"}) && ($thisData[$format{"DP"}]) && ($thisData[$format{"DP"}] ne '.')) {
		$thisDP = $thisData[$format{"DP"}];
	    }
	    if ((defined $format{"AD"}) && ($thisData[$format{"AD"}]) && ($thisData[$format{"AD"}] =~ /^[\d,]+$/)) {
		my $sumOfADs = 0;
		foreach my $ad (split(/,/,$thisData[$format{"AD"}])) {
		    $sumOfADs += $ad;
		}
		($thisDP < $sumOfADs) && ($thisDP = $sumOfADs);
	    }

	    # if depth too low or undefined for this sample, change to NOCALL
	    if ($thisDP < $filterParamsR->{"minDP"}) {
		$lineToPrint .= "\t./.";
		next;
	    }
	    # with GATK I had some issues with DP=0 calls, causing illegal divisions
	    # by zero when calculating AF, but that is now skipped above

	    # clean up Strelka GTs and calculate AF, for fracVarReads filter:
	    # Strelka makes some phased calls sometimes, homogenize as unphased
	    $thisData[$format{"GT"}] =~ s~\|~/~ ;
	    # Strelka also makes some hemizygous calls (eg when the position
	    # is in a HET deletion), makes sense but still, homogenize as HOMO
	    $thisData[$format{"GT"}] =~ s~^(\d+)$~$1/$1~;
	    # grab geno
	    my ($geno1,$geno2) = split(/\//, $thisData[$format{"GT"}]);
	    ((defined $geno1) && (defined $geno2)) ||
		die "E: a sample's genotype cannot be split: ".$thisData[$format{"GT"}]."in:\n$line\n";
	    # make sure alleles are in sorted order
	    if ($geno2 < $geno1) {
		my $genot = $geno1;
		$geno1 = $geno2;
		$geno2 = $genot;
		$thisData[$format{"GT"}] = "$geno1/$geno2";
	    }

	    my $af = '.';
	    if (($geno2 != 0) && (($geno1 == 0) || ($geno1 == $geno2))) {
		# 0/x HET or x/x HV, AD should always be there
		if ((! $thisData[$format{"AD"}]) || ($thisData[$format{"AD"}] !~ /^[\d,]+$/)) {
		    die "E: GT is HET or HV but we don't have AD or AD data is blank in:\n$line\nright after:\n$lineToPrint\n";
		}
		my @ads = split(/,/, $thisData[$format{"AD"}]);
		# $geno2 is always the index of the VAR (thanks to sorting above)
		my $fracVarReads = $ads[$geno2] / $thisDP ;
		if ($fracVarReads < $filterParamsR->{"minAF"}) {
		    # fracVarReads too low, change to NOCALL
		    $lineToPrint .= "\t./.";
		    next;
		}
		else {
		    # keeping, round AF to nearest float with 2 decimals
		    $af = sprintf("%.2f",$fracVarReads);
		}
	    }
	    # else this is HR or x/y, minAF doesn't apply, use default AF='.'

	    # we have $thisDP and $af , fix blatantly wrong calls
	    if (($thisDP >= $filterParamsR->{"minDP_HV"}) && ($geno1 == 0) &&
		($af ne '.') && ($af >= $filterParamsR->{"minAF_HV"})) {
		# change to HV
		$thisData[$format{"GT"}] = "$geno2/$geno2";
		$fixedToHV++;
		if ($verbose >= 2) {
		    # warn with chrom pos ref > alts sample dp af
		    warn "I: fix to HV, $data[0]:$data[1] $data[3] > $data[4] sample ".($i-9)." DP=$thisDP AF=$af\n";
		}
	    }
	    if (($thisDP >= $filterParamsR->{"minDP_HET"}) && ($geno1 != 0) && ($af ne '.') && 
		($af >= $filterParamsR->{"minAF_HET"}) && ($af <= $filterParamsR->{"maxAF_HET"})) {
		# change to HET
		$thisData[$format{"GT"}] = "0/$geno2";
		$fixedToHET++;
		if ($verbose >= 2) {
		    # warn with chrom pos ref > alts sample dp af
		    warn "I: fix to HET, $data[0]:$data[1] $data[3] > $data[4] sample ".($i-9)." DP=$thisDP AF=$af\n";
		}
	    }

	    # other filters (eg strandDisc) would go here

	    # OK data passed all filters but $thisData(GT) may have been changed
	    # -> fix GT in $data and set $keepLine if not HOMOREF
	    ($data =~ s/^([^:]+):/$thisData[$format{"GT"}]:$af:/) || 
		die "cannot fix GT to $thisData[$format{'GT'}] and add AF $af after the geno in: $data\n";
	    $lineToPrint .= "\t$data";
	    ($thisData[$format{"GT"}] ne '0/0') && ($keepLine = 1);
	}
	# done with $line, print if at least one sample is not NOCALL|HOMOREF
	($keepLine) && (print $outFH "$lineToPrint\n");
    }
    # INFO with number of fixed calls in this batch, we don't care that this
    # comes out of order to stderr
    ($verbose) && ($fixedToHV) && (warn "I: fixed $fixedToHV calls from HET to HV\n");
    ($verbose) && ($fixedToHET) && (warn "I: fixed $fixedToHET calls from HV to HET\n");
}


###############
# this function waits for flagfiles to be created and "eats" the 
# corresponding tmpFiles in order, starting at 1.
# "eating" means print lines to stdout and remove the tmpfiles.
# we also watch for $tmpOutLast, a file that will tell us
# the last batch number to wait for
sub eatTmpFiles {
    (@_ == 1) || die "eatTmpFiles needs 1 arg.\n";
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
		    die "E in eatTmpFiles, flagfile $tmpOutFlag exists but cant read tmpFile $tmpOut: $!\n";
	    while(<IN>) {
		print $_;
	    }
	    close(IN);
	    (unlink($tmpOut,$tmpOutFlag) == 2) ||
		die "E in eatTmpFiles, done with tmpFile $tmpOut and flagfile $tmpOutFlag but cannot unlink them: $!\n";
	    my $now = strftime("%F %T", localtime);
	    warn("I: $now - done processing batch $nextBatch\n");
	    $nextBatch++;
	    next;
	}

	elsif (-e $tmpOutLast) {
	    open (IN, $tmpOutLast) || 
		die "E in eatTmpFiles, cannot open tmpOutLast $tmpOutLast although it exists: $!\n";
	    $lastBatch = <IN>;
	    chomp($lastBatch);
	    close(IN);
	    unlink($tmpOutLast) || 
		die "E in eatTmpFiles, cannot unlink tmpOutLast $tmpOutLast: $!\n";
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

