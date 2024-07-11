#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2024
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


# 24/03/2018
# NTM

# Read on stdin a VCF file, write to stdout a similar VCF file
# with added VEP annotations.
# The VEP command-line args and VEP plugins used are hard-coded
# in &vepCommand(), edit that sub if you want different VEP fields
# (and then you should also edit @goodVeps at the top of 05_vcf2tsv.pl)
#
# 11/08/2019: adding a cache file (VEP is slow).
# CSQ are taken from the cachefile when the variant is in it,
# otherwise run VEP and add CSQ to cache.
# If VEP and/or the VEP-provided cache are updated, script dies and explains that
# cachefile must be manually removed. It will then be rebuilt from scratch
# on the next execution.

use strict;
use warnings;
use File::Basename qw(basename);
use Getopt::Long;
use POSIX qw(strftime);
# import LOCK_* constants
use Fcntl qw(:flock);
# Storable for the cache
use Storable qw(store retrieve lock_retrieve);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


### enclosing block for the "main", to enforce that the subs don't use any
### global/file-level-my variables except $0
{
    #########################################################################
    ## options / params from the command-line

    # our VEP cachefile, with path
    my $cacheFile;

    # human genome fasta, filename with path, for VEP --hgvs
    my $genome;

    # dir containing subdirs with the data required by the VEP plugins 
    # we use (eg dbNSFP/)
    my $dataDir;

    # tmp dir, preferably on a RAMdisk or at least SSD
    my $tmpDir;

    # debug: if true don't use any annotations from our cacheFile, just compare
    # the VEP annotations with those in the cache and report any diffs to stderr
    my $debug = '';

    # vep "binary", with path if it's not in PATH
    my $vepBin = "vep";

    # number of jobs that VEP can run in parallel (--fork)
    my $vepJobs = 6;

    # help: if true just print $USAGE and exit
    my $help = '';


    my $USAGE = "\nRead on stdin a VCF file, write to stdout a similar VCF file with added VEP annotations.\n
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--cacheFile string : filename with path where this script stores its cachefile
--genome string : ref genome fasta, with path
--dataDir string : dir containing subdirs with the data required by the VEP plugins we use (eg dbNSFP)
--tmpDir string : tmp dir, must not pre-exist, will be removed after running
--debug : don't use the cacheFile, just report discrepancies between it and VEP
--vep [default $vepBin] : VEP program, with path if it's not in PATH
--jobs [default $vepJobs] : number of jobs that VEP can run in parallel
--help : print this USAGE";

    GetOptions ("cacheFile=s" => \$cacheFile,
		"genome=s" => \$genome,
		"dataDir=s" => \$dataDir,
		"tmpDir=s" => \$tmpDir,
		"debug" => \$debug,
		"vep" => \$vepBin,
		"jobs" => \$vepJobs,
		"help" => \$help)
	or die("E $0: Error in command line arguments\n$USAGE\n");

    # make sure required options were provided and sanity check them
    ($help) && die "$USAGE\n\n";

    ($cacheFile) || die "E $0: you must provide a cacheFile\n$USAGE\n";
    # if cacheFile exists we need write-access to it
    if (-e $cacheFile) {
	(-w $cacheFile) ||
	    die "E $0: provided cacheFile $cacheFile exists but isn't writable\n";
    }
    # otherwise it's OK we'll create it, but we need write-access to the subdir
    else {
	open(TEST, ">>$cacheFile") ||
	    die "E $0: provided cacheFile $cacheFile doesn't exist and we can't create it\n";
	close(TEST);
	unlink($cacheFile);
    }

    ($genome) || die "E $0: you must provide a ref genome fasta file\n";
    (-f $genome) || die "E $0: provided genome fasta file doesn't exist\n";

    ($dataDir) || 
	die "E $0: you must provide a dataDir containing the data required by the VEP plugins we use\n";
    (-d $dataDir) ||
	die "E $0: the provided dataDir $dataDir doesn't exist or isn't a dir\n";

    ($tmpDir) || die "E $0: you must provide a non-existing tmpDir\n";
    (-e $tmpDir) &&
	die "E $0: tmpDir $tmpDir exists, please rm -r $tmpDir or use another tmpDir as arg\n";
    mkdir($tmpDir) || die "E $0: cannot create tmpDir $tmpDir\n";

    (`which $vepBin` =~ /$vepBin$/) ||
	die "E $0: the VEP executable $vepBin can't be found\n";

    ($vepJobs =~ /^\d+$/) || die "E $0: --jobs must be an int, you provided $vepJobs\n";

    
    my $now = strftime("%F %T", localtime);
    warn "I $now: $0 - starting to run\n";


    ##########################################################################

    # vep --no_stats results in buggy VEP output (VEP git issue 1034), so we produce
    # stats in $vepStats and remove the file when done
    my $vepStats = "$tmpDir/vepStats.html";

    # construct full VEP command-line
    my $vepCommand = &vepCommand($vepBin, $dataDir, $genome, $vepJobs, $vepStats);

    # read $cacheFile if provided
    # $cache is a hashref. key=="chr:pos:ref:alt" (for SNVs/indels) or 
    # "chr:pos:ref:alt:end" (for CNVs), value==CSQ, + 2 special keys
    # "VEPVERSION" and "INFOCSQ" used to store the VEP version and fields
    my $cache = {};
    # grab previously cached data except in debug mode
    if ((!$debug) && (-f $cacheFile) && (! -z $cacheFile)) {
	$cache = &lock_retrieve($cacheFile) ||
	    die "E $0: cachefile $cacheFile exists and is not empty but I can't retrieve hash from it.\n";
    }

    # $cacheUpdate is a hashref, same key->value pairs as $cache but it
    # will hold new VEP CSQs for keys that were missing in $cache
    my $cacheUpdate = {};

    #################################
    # headers:
    # - send header of VCF on stdin to VEP, to make sure the cache isn't stale;
    # - if AOK, print VEP-processed header immediately to stdout;
    # - store the #CHROM header line in $headerLine, we'll need it later for VEP.
    my $headerLine = &checkHeaders($cache, $cacheUpdate, $tmpDir, $vepCommand);
    if (!$headerLine) {
	# VEP version or CSQ mismatch detected
	unlink($vepStats) || die "E $0: cannot unlink vepStats file $vepStats\n";
	rmdir($tmpDir) ||
	    warn "W $0: VEP version or INFO-CSQ mismatch but can't rmdir tmpDir $tmpDir\n";
	die "E $0:  VEP version or INFO-CSQ mismatch. If you updated your VEP cache or VEP ".
	    "options/plugins this script's cachefile is now stale, you need to rm $cacheFile ".
	    "and rebuild a fresh cache.\n";
    }

    #################################
    # data
    # - $headerLine and data lines absent from cache are sent to VEP, and result
    #   is printed to $vcfFromVep;
    # - data lines found in cache are annotated using cache and printed
    #   to $vcfFromCache;
    # - whenever we change chromosomes, start a thread that merges the fromVep
    #   and fromCache files from previous chrom and prints results to stdout.
    # Both $vcf temp files are gzipped --fast, because they can get large and
    # $tmpDir can be a smallish ramdisk.

    # filenames and FILEHANDLES are per-chromosome
    my ($vcfFromVep, $vcfFromCache);
    my ($VCFVEP, $VCFCACHE);

    # chrom we were working on in previous stdin line
    my $prevChr;

    # parse data lines from stdin
    while (my $line = <STDIN>) {
	chomp($line);
	my @f = split(/\t/,$line,-1);
	(@f >= 8) || die "E $0: VCF line doesn't have >= 8 columns:\n$line\n";

	my $chr = $f[0];
	if (!defined $prevChr) {
	    # this is the first data line
	    $vcfFromVep = "$tmpDir/vcfFromVep_$chr.vcf.gz";
	    (-e $vcfFromVep) &&
		die "E $0: vcfFromVep $vcfFromVep already exists\n";
	    open($VCFVEP, "| $vepCommand | gzip -c --fast > $vcfFromVep") ||
		die "E $0: cannot run VEP on new data for chrom $chr\n";
	    $vcfFromCache = "$tmpDir/vcfFromCache_$chr.vcf.gz";
	    (-e $vcfFromCache) &&
		die "E $0: vcfFromCache $vcfFromCache already exists\n";
	    open($VCFCACHE, "| gzip -c --fast > $vcfFromCache") ||
		die "E $0: cannot open gzip pipe to vcfFromCache $vcfFromCache\n";
	    # send minimal header to VEP
	    print $VCFVEP $headerLine;
	    $prevChr = $chr;
	}
	elsif ($prevChr ne $chr) {
	    # we changed chrom, start merging and printing results from prevChr
	    close($VCFVEP);
	    close($VCFCACHE);
	    $now = strftime("%F %T", localtime);
	    warn "I $now: $0 - finished parsing/processing chrom $prevChr\n";
	    # merge prevChr files
	    &mergeAndPrint($vcfFromVep, $vcfFromCache, $cacheUpdate);
	    unlink($vcfFromVep, $vcfFromCache);
	    # set up new outputs
	    $vcfFromVep = "$tmpDir/vcfFromVep_$chr.vcf.gz";
	    (-e $vcfFromVep) &&
		die "E $0: vcfFromVep $vcfFromVep already exists\n";
	    open($VCFVEP, "| $vepCommand | gzip -c --fast > $vcfFromVep") ||
		die "E $0: cannot run VEP on new data for chrom $chr\n";
	    $vcfFromCache = "$tmpDir/vcfFromCache_$chr.vcf.gz";
	    (-e $vcfFromCache) &&
		die "E $0: vcfFromCache $vcfFromCache already exists\n";
	    open($VCFCACHE, "| gzip -c --fast > $vcfFromCache") ||
		die "E $0: cannot open gzip pipe to vcfFromCache $vcfFromCache\n";
	    # send minimal header to VEP
	    print $VCFVEP $headerLine;
	    $prevChr = $chr;
	}

	# no else: process $line whether we changed chrom or not
	# key: chrom:pos:ref:alt for SNVs/indels, chrom:pos:ref:alt:end for CNVs
	my $key = "$f[0]:$f[1]:$f[3]:$f[4]";
	if (($f[4] eq '<DUP>') || ($f[4] eq '<DEL>')) {
	    # CNV: grab END
	    ($f[7] =~ /END=(\d+)/) || 
		die "E $0: cannot grab END from CNV line:\n$line\n";
	    $key .= ":$1";
	}
	if (defined $cache->{$key}) {
	    # just copy CSQ from cache
	    ($f[7]) && ($f[7] ne '.') && ($f[7] .= ';');
	    ($f[7]) && ($f[7] eq '.') && ($f[7] = '');
	    $f[7] .= 'CSQ='.$cache->{$key};
	    print $VCFCACHE join("\t",@f)."\n";
	}
	else {
	    print $VCFVEP "$line\n";
	}
    }

    # merge VCFs for last chrom, unless stdin was completely empty
    if ($prevChr) {
	$now = strftime("%F %T", localtime);
	warn "I $now: $0 - finished parsing/processing chrom $prevChr\n";
	close($VCFVEP);
	close($VCFCACHE);
	# merge last chr files
	&mergeAndPrint($vcfFromVep, $vcfFromCache, $cacheUpdate);
	unlink($vcfFromVep, $vcfFromCache);

	# update cache with new CSQs
	&updateCache($cacheFile, $cacheUpdate, $debug);
    }

    # clean up
    unlink($vepStats) || die "E $0: cannot unlink vepStats file $vepStats\n";
    rmdir($tmpDir) || die "E $0: cannot rmdir tmpdir $tmpDir\n";

    $now = strftime("%F %T", localtime);
    warn "I $now: $0 - ALL DONE, completed successfully!\n";
}


##########################################################################
### SUBS
##########################################################################

#######################
# construct and return the full VEP command-line
sub vepCommand {
    (@_ == 5) || die "E $0: vepCommand needs 5 args";
    my ($vepBin, $dataDir, $genome, $vepJobs, $vepStats) = @_;

    # --no_stats results in buggy VEP output (VEP git issue 1034), so we produce
    # stats in $vepStats and the caller must remove the file when done

    # VEP command, reading on stdin and printing to stdout
    my $vepCommand = $vepBin;
    $vepCommand .= " --offline --format vcf --vcf" ;
    # cache to use: refseq, merged, or ensembl (default)
    # $vepCommand .= " --merged" ;
    $vepCommand .= " --force_overwrite";
    $vepCommand .= " --stats_file $vepStats";
    $vepCommand .= " --allele_number"; # for knowing which CSQ annotates which ALT
    $vepCommand .= " --canonical --biotype --xref_refseq --symbol --mane";
    $vepCommand .= " --numbers --total_length  --variant_class";
    # report where the variant lies in the miRNA secondary structure
    $vepCommand .= " --mirna";
    # reduce default distance (5000) for annotating [upstream|downstream]_gene_variant
    $vepCommand .= " --distance 1000";
    # right-align indels before consequence calculation, see:
    # https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#shifting
    $vepCommand .= " --shift_3prime 1";
    $vepCommand .= " --sift b --polyphen b";
    # --af is 1KG-phase3 global AF, --af_1kg is per continent AFs
    $vepCommand .= " --af --af_1kg --af_gnomade --af_gnomadg";
    $vepCommand .= " --check_existing";
    # Don't URI escape HGVS strings
    $vepCommand .= " --no_escape";
    # commenting out "--domains", it's a bit massive and in non-deterministic order
    # and we don't curently look at it
    # also removing --pubmed, don't think anyone looks at that either
    # I also tried out --regulatory but it's not really usable IMO, eg we obtain
    # ENSR and ENSM features but we can't know the target genes...
    ## other possibilities to consider: --tsl --appris 
    $vepCommand .= " --fasta $genome --hgvs";

    # plugins:
    # full --plugin string with all VEP plugins we want to use and
    # associated datafiles (with paths)
    my $vepPlugins = "";
    # CADD - dbNSFP provides it (among many other things) for coding and consensus
    # splice site variants, but not for deeper intronic variants
    my $caddPath = "$dataDir/CADD/";
    (-d $caddPath) || die "E $0: CADD datadir doesn't exist: $caddPath\n";
    $vepPlugins .= " --plugin CADD,$caddPath/whole_genome_SNVs.tsv.gz,$caddPath/gnomad.genomes.r3.0.indel.tsv.gz";
    # dbNSFP
    my $dbNsfpPath = "$dataDir/dbNSFP/";
    (-d $dbNsfpPath) || die "E $0: dbNSFP datadir doesn't exist: $dbNsfpPath\n";
    # comma-separated list of fields to retrieve from dbNSFP, there are MANY
    # possibilities, check the README in $dbNsfpPath
    my $dbNsfpFields = "MutationTaster_pred,REVEL_rankscore,CADD_raw_rankscore";
    # MetaRNN: both MetaRNN_rankscore and MetaRNN_pred: T(olerated) or D(amaging)
    $dbNsfpFields .= ",MetaRNN_rankscore,MetaRNN_pred";
    $vepPlugins .= " --plugin dbNSFP,$dbNsfpPath/dbNSFP4.3a.gz,transcript_match=1,$dbNsfpFields";
    # dbscSNV (splicing), data is with dbNSFP (same authors), specify 
    # assembly GRCh38 as second param because the plugin can't figure it out
    $vepPlugins .= " --plugin dbscSNV,$dbNsfpPath/dbscSNV1.1_GRCh38.txt.gz,GRCh38";

    # spliceAI - I installed the plugin but it also needs data, to DL that data you
    # have to create an account, provide your email and personal details... see:
    # https://github.com/Ensembl/VEP_plugins/blob/release/105/SpliceAI.pm
    my $spliceAIPath = "$dataDir/SpliceAI/";
    $vepPlugins .= " --plugin SpliceAI,".
	"snv=$spliceAIPath/spliceai_scores.raw.snv.hg38.vcf.gz,".
	"indel=$spliceAIPath/spliceai_scores.raw.indel.hg38.vcf.gz";

    # AlphaMissense, data file must be DL'd as $alphaMSfile and tabix-indexed, see:
    # https://github.com/Ensembl/VEP_plugins/blob/release/111/AlphaMissense.pm
    my $alphaMSfile = "$dataDir/AlphaMissense/AlphaMissense_hg38.tsv.gz";
    if (-f $alphaMSfile) {
	$vepPlugins .= " --plugin AlphaMissense,file=$alphaMSfile";
    }
    else {
	warn "W: $0 - AlphaMissense VEP plugin won't be used because file $alphaMSfile doesn't exist\n";
    }

    $vepCommand .= $vepPlugins;

    # --fork borks when $vepJobs==1
    ($vepJobs > 1) && ($vepCommand .= " --fork $vepJobs") ;
    # write output to stdout so we can pipe it to another program
    $vepCommand .= " -o STDOUT" ;

    return($vepCommand);
}

#######################
# checkHeaders:
# - read VCF headers from STDIN;
# - run them through VEP, make sure the cache isn't stale;
# - if AOK, print VEP-processed header to stdout;
# - if $cache was empty, populate VEPVERSION and INFOCSQ in $cacheUpdate;
# - return the #CHROM header line, or "" if the cache is stale
sub checkHeaders {
    (@_ == 4) || die "E $0: checkHeaders needs 4 args";
    my ($cache, $cacheUpdate, $tmpDir, $vepCommand) = @_;

    # #CHROM header line to return, will be set to "" if cache is stale
    my $headerLine;

    # make a small VCF containing all the headers from STDIN
    my $vcfHeader = "$tmpDir/vcfHeader.vcf";
    open(my $VCFHEAD, "> $vcfHeader") ||
	die "E $0: cannot open vcfHeader $vcfHeader for writing\n";

    while (my $line = <STDIN>) {
	($line =~ /^#/) ||
	    die "E $0 in checkHeaders: STDIN isn't a correct VCF, found line:\n$line";
	# header lines go to VCFHEAD
	print $VCFHEAD $line;
	# line #CHROM is always last header line
	($line =~ /^#CHROM/) && ($headerLine = $line) && last;
    }
    close($VCFHEAD);

    ($headerLine) ||
	die "E $0 in checkHeaders: STDIN is not a VCF, comment lines don't end with #CHROM line\n";

    # run VEP on header file
    open(my $VCFHEAD_OUT, "$vepCommand < $vcfHeader |") ||
	die "E $0: cannot run VEP on header file with:\n$vepCommand < $vcfHeader\n";
    # check that the cache matches the VEP version and has the correct VEP columns,
    # if AOK print processed headers to stdout
    my $headerToPrint = "";
    while (my $line = <$VCFHEAD_OUT>) {
	$headerToPrint .= $line;
	chomp($line);
	if ($line =~ /^##VEP=/) {
	    # make sure VEP version + DBs match the cache
	    # need to remove timestamp
	    my $lineClean = $line;
	    ($lineClean =~ s/time="[^"]+" //) ||
		die "E $0: cannot remove timestamp from ##VEP line:\n$line\n";
	    # also remove any path before /.vep/ in cache= so different users can use the
	    # same cacheFile if they have their own VEP install with same cache versions
	    $lineClean =~ s~(cache=")[^"]+(/.vep/)~$1$2~; # no "||die", eg if the user changed his VEPDIR
	    # also remove ensembl* fields, they appear in a random order that changes from run to run:
	    # https://github.com/Ensembl/ensembl-vep/issues/1540
	    $lineClean =~ s/\sensembl\S+//g ||
		die "E $0: cannot remove ensembl* version info from ##VEP line:\n$line\n";
	    if (defined $cache->{"VEPversion"}) {
		my $cacheLine = $cache->{"VEPversion"};
		if ($cacheLine ne $lineClean) {
		    # version mismatch
		    $headerLine = "";
		    warn "E: $0 - cached VEP version and ##VEP line from VCF are different:\n$cacheLine\n$lineClean\n";
		    last;
		}
	    }
	    else {
		$cacheUpdate->{"VEPversion"} = $lineClean;
	    }
	}
	elsif ($line =~ /^##INFO=<ID=CSQ/) {
	    # make sure the INFO->CSQ fields from file and cache match
	    if (defined $cache->{"INFOCSQ"}) {
		my $cacheLine = $cache->{"INFOCSQ"};
		if ($cacheLine ne $line) {
		    $headerLine = "";
		    warn "E: $0 - cacheLine and INFO-CSQ line from VCF are different:\n$cacheLine\n$line\n";
		    last;
		}
	    }
	    else {
		$cacheUpdate->{"INFOCSQ"} = $line;
	    }
	}
    }
    close($VCFHEAD_OUT);
    unlink($vcfHeader);

    # print headers unless there was a VEP version or INFO-CSQ mismatch
    ($headerLine) && (print $headerToPrint);

    # return $headerLine, whether empty or not
    return($headerLine);
}

#######################
# Merge $vcfFromVep and $vcfFromCache, printing resulting VCF to stdout.
# Header lines are ignored.
# $vcfFromVep and $vcfFromCache must both concern a single identical CHROM.
# All CSQs in $vcfFromVep are added to $cacheUpdate, a hashref similar
# to $cache: key=="chr:pos:ref:alt" for SNVs/indels and "chrom:pos:ref:alt:end"
# for CNVs, value==CSQ
sub mergeAndPrint {
    (@_ == 3) || die "E $0: mergeAndPrint needs 3 args";
    my ($vcfFromVep, $vcfFromCache, $cacheUpdate) = @_;

    open(my $FROMVEP, "gunzip -c $vcfFromVep |") || 
	die "E $0: cannot gunzip-open FROMVEP $vcfFromVep\n";
    open (my $FROMCACHE, "gunzip -c $vcfFromCache |") || 
	die "E $0: cannot gunzip-open FROMCACHE $vcfFromCache\n";

    # next data lines from each file
    my ($nextVep, $nextCache);
    # POS from each next line (undef if no more data)
    my ($nextVepPos, $nextCachePos);
    # also need ALT and END (set to -1 if absent), to obtain deterministic merged output
    my ($nextVepAlt, $nextCacheAlt, $nextVepEnd, $nextCacheEnd);
    # CHROM, must be the same for both files and all lines
    my $chr;

    #######################
    # prime with first lines, skipping headers
    do {
	$nextVep = <$FROMVEP>;
    } while(($nextVep) && ($nextVep =~ /^#/));
    if ($nextVep) {
	# grab next* values for deterministic merging of fromVep and fromCache, and
	# also populate $cacheUpdate for updating cache at the end
	my @f = split(/\t/,$nextVep);
	(@f >= 8) || die "E $0: nextVep line doesn't have >=8 columns:\n$nextVep";
	($chr,$nextVepPos,$nextVepAlt)=($f[0],$f[1],$f[4]);
	# key for $cacheUpdate: chrom:pos:ref:alt for SNVs/indels, chrom:pos:ref:alt:end for CNVs
	my $key = "$f[0]:$f[1]:$f[3]:$f[4]";
	$nextVepEnd = -1;
	if (($nextVepAlt eq '<DUP>') || ($nextVepAlt eq '<DEL>')) {
	    # CNV: grab END
	    ($f[7] =~ /END=(\d+)/) || 
		die "E $0: cannot grab END in nextVep line:\n$nextVep";
	    $nextVepEnd = $1;
	    $key .= ":$1";
	}
	($f[7] =~ /CSQ=([^;]+)/) || die "E $0: cannot grab CSQ in nextVep line:\n$nextVep";
	my $csq = $1;
	$cacheUpdate->{$key} = $csq;
    }

    $nextCache = <$FROMCACHE>;
    if ($nextCache) {
	my @f = split(/\t/,$nextCache);
	my $thisChr = $f[0];
	if (!defined $chr) {
	    $chr = $thisChr;
	}
	else {
	    ($chr eq $thisChr) ||
		die"E $0: vcfFromCache has line with bad chrom, expected $chr:\n$nextCache\n";
	}
	($nextCachePos,$nextCacheAlt) = ($f[1],$f[4]);
	$nextCacheEnd = -1;
	if (($nextCacheAlt eq '<DUP>') || ($nextCacheAlt eq '<DEL>')) {
	    # CNV: grab END
	    ($f[7] =~ /END=(\d+)/) || 
		die "E $0: cannot grab END in nextCache line:\n$nextCache";
	    $nextCacheEnd = $1;
	}
    }
    
    #######################
    # as long as there is data
    while ($nextVep || $nextCache) {
	if (($nextVep) && ((! $nextCache) ||
			   ($nextVepPos < $nextCachePos) ||
			   (($nextVepPos == $nextCachePos) && ($nextVepAlt lt $nextCacheAlt)) ||
			   (($nextVepPos == $nextCachePos) && ($nextVepAlt eq $nextCacheAlt) && ($nextVepEnd < $nextCacheEnd)) )) {
	    # still have VEP data and it must be printed now
	    print $nextVep;
	    # read next VEP line
	    $nextVep = <$FROMVEP>;
	    if ($nextVep) {
		# grab next* values for deterministic merging of fromVep and fromCache, and
		# also populate $cacheUpdate for updating cache at the end
		my @f = split(/\t/,$nextVep);
		(@f >= 8) || die "E $0: nextVep line doesn't have >=8 columns:\n$nextVep";
		my $thisChr = $f[0];
		($chr eq $thisChr) ||
		    die "E $0: vcfFromVep has line with bad chrom, expected $chr:\n$nextVep\n";
		($nextVepPos,$nextVepAlt)=($f[1],$f[4]);
		# key for $cacheUpdate: chrom:pos:ref:alt for SNVs/indels, chrom:pos:ref:alt:end for CNVs
		my $key = "$f[0]:$f[1]:$f[3]:$f[4]";
		$nextVepEnd = -1;
		if (($nextVepAlt eq '<DUP>') || ($nextVepAlt eq '<DEL>')) {
		    # CNV: grab END
		    ($f[7] =~ /END=(\d+)/) || 
			die "E $0: cannot grab END in nextVep line:\n$nextVep";
		    $nextVepEnd = $1;
		    $key .= ":$1";
		}
		($f[7] =~ /CSQ=([^;]+)/) || die "E $0: cannot grab CSQ in nextVep line:\n$nextVep";
		my $csq = $1;
		$cacheUpdate->{$key} = $csq;
	    }
	}
	else {
	    # cache data must be printed
	    print $nextCache;
	    $nextCache = <$FROMCACHE>;
	    if ($nextCache) {
		my @f = split(/\t/,$nextCache);
		my $thisChr = $f[0];
		($chr eq $thisChr) ||
		    die"E $0: vcfFromCache has line with bad chrom, expected $chr:\n$nextCache\n";
		($nextCachePos,$nextCacheAlt) = ($f[1],$f[4]);
		$nextCacheEnd = -1;
		if (($nextCacheAlt eq '<DUP>') || ($nextCacheAlt eq '<DEL>')) {
		    # CNV: grab END
		    ($f[7] =~ /END=(\d+)/) || 
			die "E $0: cannot grab END in nextCache line:\n$nextCache";
		    $nextCacheEnd = $1;
		}
	    }
	}
    }
    
    #######################
    # clean up
    close($FROMVEP);
    close($FROMCACHE);
}

#######################
# - Read cache from $cacheFile and compare with  entries from %$cacheupdate;
# - save new entries from %$cacheupdate to $cacheFile (unless $debug);
# - if hashrefs have common keys, warn if the values disagree.
# This sub is thread-safe (can be called by several jobs running in parallel).
sub updateCache {
    (@_ == 3) || die "E $0: updateCache needs 3 args";
    my ($cacheFile, $cacheUpdate, $debug) = @_;

    # return immediately if there's nothing to update
    (%$cacheUpdate) || return();
    
    # lock $cacheFile and retrieve a fresh copy (which will include any
    # updates made by other processes while this process was running), then
    # add any new entries from %$cacheUpdate and store
    # -> works correctly even if several jobs are running in parallel
    open(CACHELOCK, ">>$cacheFile") ||
	die "E $0: I want to update cacheFile but I can't open it for locking\n";
    flock(CACHELOCK, LOCK_EX) || die "E $0: cannot get lock on cachefile\n";

    my $cache = {};
    if (! -z $cacheFile) {
	$cache = &retrieve($cacheFile) ||
	    die "E $0: I can't retrieve fresh copy of cache from cachefile $cacheFile\n";
    }
    foreach my $k (keys(%$cacheUpdate)) {
	if (defined($cache->{$k})) {
	    # annotation for key $k has been added to $cacheFile while we were running
	    # (or was already there if $debug), make sure it's consistent
	    ($cache->{$k} eq $cacheUpdate->{$k}) ||
		warn "W $0: cacheFile entry for $k was added while we were running (or we are in debugVep mode)".
		" and disagrees with our entry, we keep the cached version:\n".
		$cache->{$k}."\n".$cacheUpdate->{$k}."\n";
	}
	else {
	    $cache->{$k} = $cacheUpdate->{$k};
	}
    }
    if (! $debug) {
	&store($cache, $cacheFile) || 
	    die "E $0: produced/updated cache but cannot store to cachefile $cacheFile\n";
    }

    # remove file if it exists but is empty (can happen in debug mode with 
    # no pre-existing cachefile)
    (-z $cacheFile) && unlink($cacheFile);

    # ok, release lock
    flock(CACHELOCK, LOCK_UN);
    close(CACHELOCK);
}


