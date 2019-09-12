#!/usr/bin/perl

# NTM
# 09/07/2019 (but starting from 1_filterBadCalls.pl which is older)


# Parses on stdin a Strelka GVCF file with one or more sample data columns;
# Args: see $USAGE.
# Prints to stdout a VCF file where:
# - grexomes that don't appear in $metadata "grexomeID" column are removed;
# - non-variant lines are removed;
# - the variant calls in data columns are replaced by ./. if a condition
#   is not met, or if previous call was '.' (the Strelka NOCALL);
# - lines where every sample is now ./. or 0/0 are skipped;
# - AF is added to FORMAT right after GT, and every 0/x or x/x call gets
#   for AF the fraction of variant reads (rounded to 2 decimals), HR
#   and x/y calls get '.'


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
my $batchSize = 100000;


#############################################
## options / params from the command-line

# filter cutoffs
my $minDP = 0;
my $minGQX = 0;
my $minFracVarReads = 0;

# metadata XLSX, no default
my $metadata;

# for multi-threading, need to create a tmpDir. It will
# be removed when we are done and must not pre-exist.
# To improve performance it should be on a ramdisk.
my $tmpDir = "tmpdir_filterBadCalls/";

# number of parallel jobs to run
my $numJobs = 16;

# help: if true just print $USAGE and exit
my $help = '';

GetOptions ("metadata=s" => \$metadata,
	    "minDP=i" => \$minDP,
            "minGQX=i"   => \$minGQX,
            "minFracVarReads=f"  => \$minFracVarReads,
	    "jobs=i" => \$numJobs,
	    "tmpdir=s" => \$tmpDir,
	    "help" => \$help)
    or die("Error in command line arguments\n");

my $USAGE = "Parse a Strelka GVCF on stdin, print to stdout a similar GVCF where calls that fail basic quality filters are changed to NOCALL and lines are only printed if at least one sample has a non-HR genotype call.\n
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--metadata string [no default] : patient metadata xlsx file, with path
--minDP int [$minDP] : must have DP or DPI >= minDP
--minGQX int [$minGQX] : must have GQX >= minGQX
--minFracVarReads float [$minFracVarReads] must have fraction of variant reads for the called variant AD[i]/DP >= minFracVarReads, filter is not applied if called genotype is HR or VAR1/VAR2 with 2 non-ref alleles
--tmpdir string [default = $tmpDir] : subdir where tmp files will be created (on a RAMDISK if possible), must not pre-exist and will be removed after execution
--jobs N [default = $numJobs] : number of parallel jobs=threads to run
--help : print this USAGE";

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
# parse patient metadata file

# key==existing grexome, value==1
my %grexomes = ();

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
    my ($grexCol) = (-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	# if column has no header just ignore it
	(defined $cell) || next;
	($cell->value() eq "grexomeID") &&
	    ($grexCol = $col);
    }
    ($grexCol >= 0) ||
	die "E parsing xlsx: no column title is grexomeID\n";

    foreach my $row ($rowMin+1..$rowMax) {
	my $grexome = $worksheet->get_cell($row, $grexCol)->value;
	# skip "none" lines
	($grexome eq "none") && next;
	(defined $grexomes{$grexome}) && 
	    die "E parsing xlsx: have 2 lines with grexome $grexome\n";
	$grexomes{$grexome} = 1;
    }
}

#############################################
# deal with headers

my $now = strftime("%F %T", localtime);
warn "I: $now - starting to run: ".join(" ", $0, @ARGV)."\n";

# array, same number of elements as there are columns in the #CHROM line
# (and hence in each data line), value is true iff column must be skipped
# (corresponding to grexomes that no longer exist in the metadata file, 
# eg they were dupes of other grexomes with better sequencing)
my @skippedCols = ();

# parse header, just copy it except we remove grexomes that don't exist anymore
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

	# remove grexomes that don't exist in $metadata (anymore)
	chomp($line);
	my @fields = split(/\t/,$line);
	foreach my $i (reverse(9..$#fields)) {
	    # reverse so we can splice bad columns out
	    if (! $grexomes{$fields[$i]}) {
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
    &processBatch(\@lines,$tmpOutFH,$minDP,$minGQX,$minFracVarReads,\@skippedCols);

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
# - filter cutoffs $minDP, $minGQX, $minFracVarReads
# - ref to array saying which columns to skip
sub processBatch {
    (@_ == 6) || die "E: processBatch needs 6 args\n";
    my ($linesR,$outFH,$minDP,$minGQX,$minFracVarReads,$skippedColsR) = @_;

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
	# we use either DP or DPI as long as one is present and not '.'
	(defined $format{"GQX"}) || die "no GQX key in FORMAT string for line:\n$line\n";
	(defined $format{"DP"}) || (defined $format{"DPI"}) ||
	    die "no DP or DPI key in FORMAT string for line:\n$line\n";
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
		if ((! $thisData[$format{"AD"}]) || ($thisData[$format{"AD"}] =~ /^[\.,]+$/)) {
		    die "E: GT is HET or HV but we don't have AD or AD data is blank in:\n$line\nright after:\n$lineToPrint\n";
		}
		my @ads = split(/,/, $thisData[$format{"AD"}]);
		# $geno2 is always the index of the VAR (thanks to sorting above)
		my $fracVarReads = $ads[$geno2] / $thisDP ;
		if ($fracVarReads < $minFracVarReads) {
		    # fracVarReads too low, change to NOCALL
		    $lineToPrint .= "\t./.";
		    next;
		}
		else {
		    # keeping, round AF to nearest float with 2 decimals
		    $af = sprintf("%.2f",$fracVarReads);
		}
	    }
	    # else this is HR or x/y, minFracVarReads doesn't apply, use default AF='.'
	    ($data =~ s/^([^:]+):/$thisData[$format{"GT"}]:$af:/) || 
		die "cannot add fixed GT $thisData[$format{'GT'}] and AF $af after the geno in: $data\n";

	    # other filters (eg strandDisc) would go here

	    # OK data passed all filters, print and set $keepLine if not HOMOREF
	    $lineToPrint .= "\t$data";
	    ($thisData[$format{"GT"}] ne '0/0') && ($keepLine = 1);
	}
	# done with $line, print if at least one sample is not NOCALL|HOMOREF
	($keepLine) && (print $outFH "$lineToPrint\n");
    }
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

