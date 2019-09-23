#!/usr/bin/perl

# 25/03/2018
# NTM

# Takes as arguments a $metadata xlsx file, a $candidatesFile
#  xlsx file, an $outDir and a $tmpDir that don't exist; 
# reads on stdin a fully annotated TSV file;
# makes $outDir and creates in it one gzipped TSV file per cohort.
# The cohorts are defined in $metadata.
# For each sample, any identified causal (mutation in a) gene is grabbed
# from $metadata.
# @notControls defined at the top of this script says which cohorts
# should NOT be used as negative controls for each other.
# For optimal performance $tmpDir should be on a RAMDISK (eg tmpfs).
#
# For each $cohort, the GENO columns HV/HET/OTHER/HR are modified as follows:
# - the HR GENO column is removed (but samples are COUNTed, see below).
# - we make new NEGCTRL_* columns placed immediately after the HV/HET/OTHER columns.
# - NEGCTRL_* columns list all samples falling in that GENO category and:
#   * having an identified $causalGene, whatever their cohort (except in lines where 
#     SYMBOL==$causalGene, see below), or
#   * belonging to another cohort that isn't defined in @notControls for $cohort.
# - for samples with an identified $causalGene and lines where SYMBOL==$causalGene,
#   the sample is dealt with as if it didn't have a $causalGene, ie it stays 
#   HV/HET/OTHER for his cohort, is ignored in his @notControls cohorts, and 
#   goes in NEGCTRL_* columns in other cohorts.
#
# A new KNOWN_CANDIDATE_GENE column is inserted right after SYMBOL:
# it holds the "Level" value parsed from  $candidatesFile if SYMBOL is a known 
# candidate gene for this cohort (as specified in $candidatesFile), 
# 0 otherwise. Any $causalGene from $metadata is considered a
# a known candidate gene with Level=5.
#
# New COUNT_$cohort_$geno and COUNT_NEGCTRL_$geno columns are created
# for each  GENO (HV, HET, OTHER, HR) in that order.
# These columns contain the total number of samples listed in the
# corresponding GENO column (except for HR, which has no GENO column).
# For HR we count all samples (ie don't care about @notControls or $causalGene).
# The COUNTs are inserted right after the new KNOWN_CANDIDATE_GENE column.
#
# Lines where no samples from the cohort are HV|HET (for this alt allele)
# are skipped. We rely on the fact that vcf2tsv.pl moved HV/HET genotypes
# concerning other alleles to OTHER (but we check it).


use strict;
use warnings;
use Spreadsheet::XLSX;
use POSIX qw(strftime);
use Parallel::ForkManager;


#############################################
## hard-coded stuff that shouldn't change much

# @notControls: array of arrayrefs, each arrayref holds cohorts that
# should NOT be used as neg controls for each other.
# The cohort names must match the "pathology" column of the $metadata xlsx
# (this is checked).
# NOTE: @notControls IS DUPLICATED IN extractTranscripts.pl, IF IT IS CHANGED HERE IT 
# MUST ALSO BE CHANGED THERE 
my @notControls = (["Flag","Astheno","Headless"],
		   ["Azoo","Ovo","Macro","IOP"],
		   ["Globo","Macro","Terato"]);


# max number of lines to read in a single batch. Each batch is then
# processed by a worker thread.
my $batchSize = 20000;

# number of jobs
my $numJobs = 16;


#########################################################

(@ARGV == 4) || die "needs 4 args: a patient metadata xlsx, a candidateGenes xlsx, and non-existing dirs outDir and tmpDir\n";
my ($metadata, $candidatesFile, $outDir, $tmpDir) = @ARGV;

(-e $outDir) && 
    die "found argument $outDir but it already exists, remove it or choose another name.\n";
mkdir($outDir) || die "cannot mkdir outDir $outDir\n";

(-e $tmpDir) && 
    die "found argument $tmpDir but it already exists, remove it or choose another name.\n";
mkdir($tmpDir) || die "cannot mkdir tmpDir $tmpDir\n";

my $now = strftime("%F %T", localtime);
warn "I: $now - starting to run: ".join(" ", $0, @ARGV)."\n";


#########################################################
# parse known candidate genes file

# %knownCandidateGenes: key==$cohort, value is a hashref whose keys 
# are gene names and values are the "Level" from $candidatesFile,
# or 5 if the gene is "Causal" for a $cohort patient in $metadata.
# I use %knownCandidatesSeen (defined below) to sanity-check the lists: any gene 
# name that is never seen will be reported to stderr (and probably a typo needs fixing).
my %knownCandidateGenes = ();

(-f $candidatesFile) ||
    die "E: the supplied candidates file $candidatesFile doesn't exist\n";
{
    my $workbook = Spreadsheet::XLSX->new("$candidatesFile");
    (defined $workbook) ||
	die "E when parsing xlsx\n";
    ($workbook->worksheet_count() == 1) ||
	die "E parsing xlsx: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($pathoCol, $geneCol,$levelCol) = (-1,-1,-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	($cell->value() eq "pathology") &&
	    ($pathoCol = $col);
	($cell->value() eq "Candidate gene") &&
	    ($geneCol = $col);
	($cell->value() eq "Level") &&
	    ($levelCol = $col);
     }
    ($pathoCol >= 0) ||
	die "E parsing xlsx: no col title is pathology\n";
    ($geneCol >= 0) ||
	die "E parsing xlsx: no col title is Candidate gene\n";
    ($levelCol >= 0) ||
	die "E parsing xlsx: no col title is Level\n";
    
    foreach my $row ($rowMin+1..$rowMax) {
	my $cohort = $worksheet->get_cell($row, $pathoCol)->unformatted();
	my $gene = $worksheet->get_cell($row, $geneCol)->unformatted();
	my $level = $worksheet->get_cell($row, $levelCol)->unformatted();

	# clean up $gene
	$gene =~ s/^\s+//;
	$gene =~ s/\s+$//;

	(defined $knownCandidateGenes{$cohort}) ||
	    ($knownCandidateGenes{$cohort} = {});
	(defined $knownCandidateGenes{$cohort}->{$gene}) && 
	    die "E parsing candidatesFile xlsx: have 2 lines with same gene $gene and pathology $cohort\n";
	$knownCandidateGenes{$cohort}->{$gene} = $level;
    }
}


#########################################################
# parse patient metadata file

# key==sample id, value is the $cohort this sample belongs to
my %sample2cohort = ();
# cohort names
my @cohorts = ();
# causal gene, key==sample id, value == HGNC gene name
my %sample2causal = ();

(-f $metadata) ||
    die "E: the supplied metadata file doesn't exist\n";
{
    # for cohort names we use a temp hash to avoid redundancy
    my %cohorts;
    my $workbook = Spreadsheet::XLSX->new("$metadata");
    (defined $workbook) ||
	die "E when parsing xlsx\n";
    ($workbook->worksheet_count() == 1) ||
	die "E parsing xlsx: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($grexCol, $cohortCol,$causalCol) = (-1,-1,-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	# if column has no header just ignore it
	(defined $cell) || next;
	($cell->value() eq "grexomeID") &&
	    ($grexCol = $col);
	($cell->value() eq "pathology") &&
	    ($cohortCol = $col);
	($cell->value() eq "Causal gene") &&
	    ($causalCol = $col);
     }
    ($grexCol >= 0) ||
	die "E parsing xlsx: no column title is grexomeID\n";
    ($cohortCol >= 0) ||
	die "E parsing xlsx: no col title is pathology\n";
    ($causalCol >= 0) ||
	die "E parsing xlsx: no col title is Causal gene\n";
    
    foreach my $row ($rowMin+1..$rowMax) {
	my $grexome = $worksheet->get_cell($row, $grexCol)->value;
	# skip "none" lines
	($grexome eq "none") && next;
	(defined $sample2cohort{$grexome}) && 
	    die "E parsing xlsx: have 2 lines with grexome $grexome\n";
	my $cohort = $worksheet->get_cell($row, $cohortCol)->value;
	$sample2cohort{$grexome} = $cohort;
	$cohorts{$cohort} = 1;
	if ($worksheet->get_cell($row, $causalCol)) {
	    my $causal = $worksheet->get_cell($row, $causalCol)->value;
	    # clean up a bit, remove leading or trailing whitespace
	    $causal =~ s/^\s+//;
	    $causal =~ s/\s+$//;
	    $sample2causal{$grexome} = $causal;
	    # add to knownCandidateGenes with level 5
	    (defined $knownCandidateGenes{$cohort}) || ($knownCandidateGenes{$cohort} = {});
	    $knownCandidateGenes{$cohort}->{$causal} = 5;
	}
    }
    @cohorts = sort(keys(%cohorts));
}

# for sanity-checking known candidate genes:
# fill this now so we also check the causal genes from $metadata
my %knownCandidatesSeen;
foreach my $c (keys(%knownCandidateGenes)) {
    foreach my $gene (keys(%{$knownCandidateGenes{$c}})) {
	$knownCandidatesSeen{$gene} = 0;
    }
}

#########################################################
# check @notControls cohort names and store in %notControls hash

# %notControls: key is a cohort name, value is a hashref
# with keys == cohorts that shouldn't be used as negative 
# controls for this cohort, value==1
my %notControls = ();

foreach my $notConR (@notControls) {
    foreach my $cohort (@$notConR) {
	(grep($cohort eq $_, @cohorts)) ||
	    die "E in extractCohorts: cohort $cohort from notControls is not in cohorts @cohorts\n";
	(defined $notControls{$cohort}) || ($notControls{$cohort} = {});
	foreach my $notC (@$notConR) {
	    ($notC eq $cohort) && next;
	    $notControls{$cohort}->{$notC} = 1;
	}
    }
}

#########################################################

# array of filehandles open for writing, one for each cohort, same indexes as @cohorts
my @outFHs;
foreach my $cohorti (0..$#cohorts) {
    my $outFile = "$outDir/$cohorts[$cohorti].csv.gz";
    open (my $FH, "| gzip -c > $outFile") || die "cannot gzip-open $outFile for writing";
    $outFHs[$cohorti] = $FH ;
}

#########################################################
# headers

my $header = <STDIN>;
chomp($header);
my @headers = split(/\t/, $header);

# the genotype categories (don't change this)
my @genoCategories = ("HV","HET","OTHER","HR");

# useful columns: SYMBOL and genos in @genoCategories order
my $symbolCol;
my @genoCols;
foreach my $i (0..$#headers) {
    ($headers[$i] eq "SYMBOL") && ($symbolCol = $i);
    foreach my $gi (0..$#genoCategories) {
	($headers[$i] eq $genoCategories[$gi]) && 
	    ($genoCols[$gi] = $i);
    }
}
($symbolCol) || die "could not find SYMBOL in headers\n";
foreach my $gi (0..$#genoCategories) {
    ($genoCols[$gi]) || die "cound not find $genoCategories[$gi] in headers\n";
    ($genoCols[$gi] == ($genoCols[0]+$gi)) || 
	die "E: GENO columns are not subsequent and in genoCategories order, the code relies on this\n";
    # actually maybe the code just requires that they be 
    # consecutive columns and start with HV (not sure)
}

# print new headers
foreach my $cohorti (0..$#cohorts) {
    # we always want to keep the first column and 
    # this simplifies things (no \t)
    my $toPrint = "$headers[0]";
    foreach my $i (1..$#headers) {
	if ($i == $symbolCol) {
	    $toPrint .= "\t$headers[$i]";
	    # KNOWN_CANDIDATE_GENE and COUNTs go right after SYMBOL
	    $toPrint .= "\tKNOWN_CANDIDATE_GENE";
	    foreach my $geno (@genoCategories) {
		$toPrint .= "\tCOUNT_".$cohorts[$cohorti]."_$geno";
	    }
	    # COUNT_NEGCTRL_* right after
	    foreach my $geno (@genoCategories) {
		$toPrint .= "\tCOUNT_NEGCTRL_$geno";
	    }
	}
	elsif (($i == $genoCols[0]) || ($i == $genoCols[1]) || ($i == $genoCols[2])) {
	    # HV/HET/OTHER
	    $toPrint .= "\t$headers[$i]\tNEGCTRL_$headers[$i]";
	}
	elsif ($i == $genoCols[3]) {
	    # HR is not printed
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
my $pm = new Parallel::ForkManager($numJobs);

# need a tmp file for listing the known candidates
my $tmpFileCandidatesSeen = "$tmpDir/allCandidates.seen";
open(my $knownCandidatesSeenFH, "> $tmpFileCandidatesSeen") ||
    die "E: cannot open tmpFileCandidatesSeen $tmpFileCandidatesSeen for writing\n";
# spawn a child process that waits for workers to finish producing batches,
# and prints the tmpfiles to @outFHs in correct order, cleaning up behind 
# itself. 
# Also prints seen candidate genes to $knownCandidatesSeenFH.
if (! $pm->start) {
    &eatTmpFiles($tmpDir,\@cohorts,\@outFHs,$knownCandidatesSeenFH);
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
    # NOTE: IF YOU CHANGE the tmp filenames below ($tmpOut, $tmpSeenFile, $tmpOutFlag),
    # you MUST EDIT &eatTmpFiles()

    # create tmp output filehandles for this batch
    my @tmpOutFHs = ();
    foreach my $i (0..$#cohorts) {
	my $tmpOut = "$tmpDir/$batchNum.$cohorts[$i].tsv";
	open(my $outFH, "> $tmpOut") || die "cannot open $tmpOut for writing\n";
	push(@tmpOutFHs, $outFH);
    }
    # create tmp output filehandle for candidatesSeen
    my $tmpSeenFile = "$tmpDir/$batchNum.seen";
    open(my $tmpSeenFH, "> $tmpSeenFile") || die "cannot open $tmpSeenFile for writing\n";

    # process this batch
    &processBatch(\@lines,\%knownCandidateGenes,\%sample2cohort,\@cohorts,
		  \%sample2causal,\%notControls,$symbolCol,\@genoCols,\@tmpOutFHs,$tmpSeenFH);

    # done, close tmp FHs and create flag-file
    foreach my $outFH (@tmpOutFHs,$tmpSeenFH) {
	close($outFH) || die "cannot close tmp outFH $outFH\n";
    }
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

foreach my $fh (@outFHs) {
    close($fh);
}
close($knownCandidatesSeenFH);

open(IN, "$tmpFileCandidatesSeen") ||
    die "cannot open tmpFileCandidatesSeen $tmpFileCandidatesSeen for reading";
while (my $gene = <IN>) {
    chomp($gene);
    $knownCandidatesSeen{$gene} = 1;
}
close(IN);
(unlink($tmpFileCandidatesSeen) == 1) ||
    die "E: cannot unlink tmpFileCandidatesSeen $tmpFileCandidatesSeen: $!\n";
foreach my $gene (keys(%knownCandidatesSeen)) {
    ($knownCandidatesSeen{$gene}) ||
	warn "W: \"known candidate gene\" $gene was never seen!! typo in metadata or candidates xlsx files?\n";
}

$now = strftime("%F %T", localtime);
rmdir($tmpDir) || 
    die "E: $now - all done but cannot rmdir tmpDir $tmpDir, why? $!\n";
warn "I: $now - DONE running $0\n";




#############################################
## subs

# process a batch of lines
# args:
# - ref to array of chomped lines
# - ref to %knownCandidateGenes
# - refs to %sample2cohort, to @cohorts, and to %sample2causal
# - ref to %notControls
# - $symbolCol, the column index of the SYMBOL column (in data lines)
# - ref to @genoCols holding column indexes of GENO columns (in @genoCategories order)
# - $tmpOutFilesR, ref to array of filehandles open for writing, one for each cohort,
#   same indexes as @cohorts
# - $tmpSeen, a filehandle open for writing, we will print one line per different
#   candidate gene seen in this batch of lines
sub processBatch {
    (@_ == 10) || die "E: processBatch needs 10 args\n";
    my ($linesR,$knownCandidateGenesR,$sample2cohortR,$cohortsR,$sample2causalR,
	$notControlsR,$symbolCol,$genoColsR,$tmpOutFilesR,$tmpSeen) = @_;

    # key == known candidate gene seen in this batch of lines, value==1
    my %candidatesSeen = ();

    foreach my $line (@$linesR) {
	my @fields = split(/\t/, $line, -1) ;

	# $symbol doesn't depend on cohorts
	my $symbol = $fields[$symbolCol];
	
      COHORT:
	foreach my $cohorti (0..$#$cohortsR) {
	    my $cohort = $cohortsR->[$cohorti];
	    # build array of 8 counts for $cohort: HV,HET,OTHER,HR and again for NEGCTRLs
	    my @counts = (0) x 8;
	    # also build array of 6 GENO columns: HV,NEGCTRL_HV,HET,NEGCTRL_HET,OTHER,NEGCTRL_OTHER
	    my @genos = ("") x 6;

	    # parse data
	    foreach my $gi (0..3) {
		my @genoData = split(/\|/,$fields[$genoColsR->[$gi]]);
		# sanity: at most one genotype except for OTHER column
		(@genoData <= 1) || ($gi==2) || 
		    die "E: more than one genoData for genotype $genoColsR->[$gi], impossible. Line:\n$line\n";
		foreach my $genoData (@genoData) {
		    ($genoData =~ /^(\d+\/\d+)~([^~\|]+)$/) ||
			die "E: cannot parse GENOS data $genoData in line:\n$line\n";
		    # $geno is the genotype (eg 1/1 or 0/2)
		    my $geno = $1;
		    my @samples = split(/,/,$2);
		    # @goodSamples will hold samples that should be counted for $cohort
		    # @badSamples will hold samples that should be counted as NEGCTRLs for $cohort
		    my @goodSamples = ();
		    my @badSamples = ();
		    foreach my $sample (@samples) {
			my $grexome = $sample;
			# remove trailing [DP:AF] if it's there (allowing AF > 1 for Strelka bug)
			$grexome =~ s/\[\d+:\d+\.\d\d\]$//;
			# sanity check
			($grexome =~ /^grexome\d+$/) || die "E grexome id $grexome illegal, sample was $sample\n";
			if ($sample2cohortR->{$grexome} eq $cohort) {
			    # $sample belongs to cohort
			    if (($gi == 3) || (! defined $sample2causalR->{$grexome}) || ($sample2causalR->{$grexome} eq $symbol)) {
				# we are HR or sample has no causal gene or it's the current gene
				push(@goodSamples,$sample);
			    }
			    else {
				# HV/HET/OTHER geno and this sample has a causal gene but not this gene:
				# this sample counts as NEGCTRL
				push(@badSamples,$sample);
			    }
			}
			elsif (($gi==3) || (! defined ${$notControlsR->{$cohort}}{$sample2cohortR->{$grexome}})) {
			    # sample is HR, or it is from another cohort that can be used as control 
			    # for $cohort
			    push(@badSamples,$sample);
			}
			elsif ((defined $sample2causalR->{$grexome}) && ($sample2causalR->{$grexome} ne $symbol)) {
			    # sample is from a notControls cohort but it has a causal gene (and it's not this gene)
			    push(@badSamples,$sample);
			}
			# else sample is from a different cohort but it's in @notControls, and it doesn't have
			# a causal gene (or it does but it is the current gene): ignore this sample, ie NOOP
		    }
		    
		    # OK, store counts and GENOs (careful with indexes)
		    if (@goodSamples) {
			$counts[$gi] += scalar(@goodSamples);
			if ($gi < 3) {
			    # don't store the GENO for HR ie gi==3
			    ($genos[$gi * 2]) && ($genos[$gi * 2] .= '|');
			    $genos[$gi * 2] .= "$geno~".join(',',@goodSamples);
			}
		    }
		    if (@badSamples) {
			$counts[$gi + 4] += scalar(@badSamples);
			if ($gi < 3) {
			    ($genos[1 + $gi * 2]) && ($genos[1 + $gi * 2] .= '|');
			    $genos[1 + $gi * 2] .= "$geno~".join(',',@badSamples);
			}
		    }
		}

		# if we just finished with HV and HET but there's no sample in either,
		# we can skip this line in this cohort
		if (($gi == 1) && ($counts[0] == 0) && ($counts[1] == 0)) {
		    next COHORT;
		}
		# otherwise: parse OTHER and HR data
	    }

	    # OK we have some data to print for $cohort
	    # always keep first field 
	    my $toPrint = "$fields[0]";
	    foreach my $i (1..$#fields) {
		if ($i == $symbolCol) {
		    # prepend apostrophe-space to gene names to avoid excel corrupting everything
		    $toPrint .= "\t\' $fields[$i]\t";
		    if (($knownCandidateGenesR->{$cohort}) && (my $level = $knownCandidateGenesR->{$cohort}->{$fields[$i]})) {
			$toPrint .= $level;
			$candidatesSeen{$fields[$i]} = 1;
		    }
		    else {
			# if not a known candidate use zero
			$toPrint .= "0";
		    }
		    # print all COUNTs
		    $toPrint .= "\t".join("\t",@counts);
		}
		elsif ($i == $genoColsR->[0]) {
		    # HV -> print all 6 GENO columns
		    # NOTE: we rely on the fact that the GENOs are consecutive and start with HV
		    $toPrint .= "\t".join("\t",@genos);
		}
		elsif (! grep(/^$i$/, @$genoColsR)) {
		    $toPrint .= "\t$fields[$i]";
		}
	    }
	    print { $tmpOutFilesR->[$cohorti] } "$toPrint\n";
	}
    }

    # done with this batch, print any seen candidate genes to $tmpSeen
    foreach my $gene (keys(%candidatesSeen)) {
	print $tmpSeen "$gene\n";
    }
}


###############
# this function waits for flagfiles to be created and "eats" the 
# corresponding tmpFiles in order, starting at 1.
# "eating" means print lines to the relevant $outFH or to $knownCandidatesSeenFH
# and remove the tmpfiles.
# we also watch for $tmpOutLast, a file that will tell us
# the last batch number to wait for
sub eatTmpFiles {
    (@_ == 4) || die "eatTmpFiles needs 4 args.\n";
    my ($tmpDir,$cohortsR,$outFHsR,$knownCandidatesSeenFH) = @_;

    # NOTE: all tmp filenames (eg $tmpOutLast) are hard-coded here and 
    # MUST MATCH those created by the worker threads.

    # when created, $tmpOutLast will contain the number of the last batch
    my $tmpOutLast = "$tmpDir/lastBatch";
    # $lastBatch remains undefined until we read it from $tmpOutLast
    my $lastBatch;

    # next batch to eat
    my $nextBatch = 1;

    while(1) {
	my $tmpSeenFile = "$tmpDir/$nextBatch.seen";
	my $tmpOutFlag = "$tmpDir/$nextBatch.done";

	if (-e $tmpOutFlag) {
	    # next batch of tmp files are ready
	    foreach my $i (0..$#$cohortsR) {
		my $tmpFile = "$tmpDir/$nextBatch.$cohortsR->[$i].tsv";
		open (IN, $tmpFile) || 
		    die "E in eatTmpFiles, flagfile $tmpOutFlag exists but cant read tmpFile $tmpFile: $!\n";
		while(<IN>) {
		    print {$outFHsR->[$i]} $_;
		}
		close(IN);
		(unlink($tmpFile) == 1) ||
		    die "E in eatTmpFiles, done with tmpFile $tmpFile but cannot unlink it: $!\n";
	    }
	    open(IN, $tmpSeenFile) ||
		die "E in eatTmpFiles, cannot open tmpSeenFile $tmpSeenFile: $!\n";
	    while(<IN>) {
		print $knownCandidatesSeenFH $_;
	    }
	    (unlink($tmpSeenFile,$tmpOutFlag) == 2) ||
		die "E in eatTmpFiles, done with files for batch $nextBatch but cannot unlink (both of) tmpSeen / tmpOutFlag: $!\n";

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

