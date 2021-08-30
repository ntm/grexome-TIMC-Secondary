#!/usr/bin/perl

# 24/03/2018
# NTM

# Read on stdin a VCF file, write to stdout a similar VCF file
# with added VEP annotations.
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


##########################################################################
## hard-coded stuff that shouldn't change much

# vep executable, with path if it's not in PATH
my $vepBin = "vep";

# number of jobs that VEP can run in parallel (--fork)
my $vepJobs = 6;

# this script also hard-codes the VEP command-line args, VEP plugins used,
# etc... this can be adjusted in the "VEP command-line" section below


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

# help: if true just print $USAGE and exit
my $help = '';


my $USAGE = "\nRead on stdin a VCF file, write to stdout a similar VCF file with added VEP annotations.\n
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--cacheFile string : filename with path where this script stores its cachefile
--genome string : ref genome fasta, with path
--dataDir string : dir containing subdirs with the data required by the VEP plugins we use (eg dbNSFP)
--tmpDir string : tmp dir, must not pre-exist, will be removed after running
--debug : don't use the cacheFile, just report discrepancies between it and VEP
--help : print this USAGE";

GetOptions ("cacheFile=s" => \$cacheFile,
	    "genome=s" => \$genome,
	    "dataDir=s" => \$dataDir,
	    "tmpDir=s" => \$tmpDir,
	    "debug" => \$debug,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($cacheFile) || die "E $0: you must provide a cacheFile\n";
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

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";

#########################################################################
# construct the full VEP command-line

# full --plugin string with all VEP plugins we want to use and
# associated datafiles (with paths)
my $vepPlugins = "";
{
    # CADD - no longer needed, dbNSFP provides it (among many other things)
    #my $caddPath = "$dataDir/CADD/";
    # (-d $caddPath) || die "E $0: CADD datadir doesn't exist: $caddPath\n";
    #$vepPlugins .= " --plugin CADD,$caddPath/whole_genome_SNVs.tsv.gz,$caddPath/InDels.tsv.gz";

    # dbNSFP
    my $dbNsfpPath = "$dataDir/dbNSFP/";
    (-d $dbNsfpPath) || die "E $0: dbNSFP datadir doesn't exist: $dbNsfpPath\n";
    # comma-separated list of fields to retrieve from dbNSFP, there are MANY
    # possibilities, check the README in $dbNsfpPath
    my $dbNsfpFields = "MutationTaster_pred,REVEL_rankscore,CADD_raw_rankscore";
    # MetaRNN: both MetaRNN_rankscore and MetaRNN_pred: T(olerated) or D(amaging)
    $dbNsfpFields .= ",MetaRNN_rankscore,MetaRNN_pred";
    $vepPlugins .= " --plugin dbNSFP,$dbNsfpPath/dbNSFP4.2a.gz,$dbNsfpFields";

    # dbscSNV (splicing), data is with dbNSFP (same authors), specify 
    # assembly GRCh38 as second param because the plugin can't figure it out
    $vepPlugins .= " --plugin dbscSNV,$dbNsfpPath/dbscSNV1.1_GRCh38.txt.gz,GRCh38";
}


# VEP command, reading on stdin and printing to stdout
my $vepCommand = $vepBin;
$vepCommand .= " --offline --format vcf --vcf" ;
# cache to use: refseq, merged, or ensembl (default)
# $vepCommand .= " --merged" ;
$vepCommand .= " --force_overwrite";# --no_stats" ;
$vepCommand .= " --allele_number"; # for knowing which CSQ annotates which ALT
$vepCommand .= " --canonical --biotype --xref_refseq --symbol --mane";
$vepCommand .= " --numbers --total_length  --variant_class" ;
$vepCommand .= " --sift b --polyphen b";
# --af is 1KG-phase3 global AF, --af_1kg is per continent AFs
$vepCommand .= " --af --af_1kg --af_gnomad";
$vepCommand .= " --pubmed --check_existing";
# Don't URI escape HGVS strings
$vepCommand .= " --no_escape";
# commenting out "--domains", it's a bit massive and in non-deterministic order
# and we don't curently look at it
## other possibilities to consider: --regulatory --tsl --appris 
$vepCommand .= " --fasta $genome --hgvs";
# plugins:
$vepCommand .= $vepPlugins;
# --fork borks when $vepJobs==1
($vepJobs > 1) && ($vepCommand .= " --fork $vepJobs") ;
# write output to stdout so we can pipe it to another program
$vepCommand .= " -o STDOUT" ;


##########################################################################
# parse VCF on stdin.
# - headers and data lines absent from cache are printed to 
#   $vcf4vep (gzipped --fast);
# - data lines found in cache are annotated using cache 
#   and printed to $vcfFromCache.
my $vcf4vep = "$tmpDir/vcf4vep.vcf.gz";
open(VCF4VEP, "| gzip -c --fast > $vcf4vep") ||
    die "E $0: cannot open gzip pipe to vcf4vep $vcf4vep\n";
my $vcfFromCache = "$tmpDir/vcfFromCache.vcf.gz";
open(VCFCACHE, "| gzip -c --fast > $vcfFromCache") ||
    die "E $0: cannot open gzip pipe to vcfFromCache $vcfFromCache\n";

# a small VCF containing only the headers is also made, for testing
# the VEP version etc...
my $vcf4vepTest = "$tmpDir/vcf4vepVersionTest.vcf";
open(VEPTEST, "> $vcf4vepTest") ||
    die "E $0: cannot open vcf4vepTest $vcf4vepTest for writing\n";

# $cache is a hashref. key=="chr:pos:ref:alt", value==CSQ
my $cache = {};
# grab previously cached data except in debug mode
if ((!$debug) && (-f $cacheFile) && (! -z $cacheFile)) {
    $cache = &lock_retrieve($cacheFile) ||
	die "E $0: cachefile $cacheFile exists and is not empty but I can't retrieve hash from it.\n";
}
# $cacheUpdate is a hashref, similar to $cache but holding info
# that wasn't found in $cache; key=="chr:pos:ref:alt", value==CSQ
my $cacheUpdate = {};

# header
while (my $line = <STDIN>) {
    # header lines go to VCF4VEP and also to VEPTEST
    print VCF4VEP $line;
    print VEPTEST $line;
    # line #CHROM is always last header line
    ($line =~ /^#CHROM/) && last;
}

# run VEP on the small test file
close(VEPTEST);
open(VEPTEST_OUT, "$vepCommand < $vcf4vepTest |") ||
    die "E $0: cannot run VEP on testfile with:\n$vepCommand < $vcf4vepTest\n";
# check that the cache matches the VEP and cache versions and has the correct VEP columns
while (my $line = <VEPTEST_OUT>) {
    chomp($line);
    if ($line =~ /^##VEP=/) {
	# make sure VEP version + DBs match the cache
	# need to remove timestamp
	my $lineClean = $line;
	($lineClean =~ s/time="[^"]+" //) ||
	    die "E $0: cannot remove timestamp from ##VEP line:\n$line\n";
	# also remove any path before /.vep/ in cache= so different users can use the
	# same $cacheFile if they have their own VEP install with same cache versions
	$lineClean =~ s~(cache=")[^"]+(/.vep/)~$1$2~; # no "||die", eg if the user changed his VEPDIR
	if (defined $cache->{"VEPversion"}) {
	    my $cacheLine = $cache->{"VEPversion"};
	    if ($cacheLine ne $lineClean) {
		# version mismatch, clean up and die
		close(VCF4VEP);
		close(VCFCACHE);
		unlink($vcf4vep,$vcfFromCache,$vcf4vepTest);
		rmdir($tmpDir) || warn "W $0: VEP version mismatch but can't rmdir tmpDir $tmpDir\n";
		die "E: $0 - cached VEP version and ##VEP line from VCF are different:\n$cacheLine\n$lineClean\n".
		    "if you updated your VEP cache this script's cachefile is now stale, you need to rm $cacheFile".
		    " (or change the cacheFile in $0)\n\n";
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
		close(VCF4VEP);
		close(VCFCACHE);
		unlink($vcf4vep,$vcfFromCache,$vcf4vepTest);
		rmdir($tmpDir) || warn "W $0: INFO-CSQ mismatch but can't rmdir tmpDir $tmpDir\n";
		die "E: $0 - cacheLine and INFO-CSQ line from VCF are different:\n$cacheLine\n$line\n".
		    "if you really want to use this vcf from STDIN you need to rm $cacheFile (or change the cacheFile in $0)\n\n";
	    }
	}
	else {
	    $cacheUpdate->{"INFOCSQ"} = $line;
	}
    }
}
close(VEPTEST_OUT);
unlink($vcf4vepTest);


# data
while (my $line = <STDIN>) {
    chomp($line);
    my @f = split(/\t/,$line,-1);
    (@f >= 8) || die "E $0: VCF line doesn't have >= 8 columns:\n$line\n";
    # key: chrom:pos:ref:alt
    my $key = "$f[0]:$f[1]:$f[3]:$f[4]";
    if (defined $cache->{$key}) {
	# just copy CSQ from cache
	($f[7]) && ($f[7] ne '.') && ($f[7] .= ';');
	($f[7]) && ($f[7] eq '.') && ($f[7] = '');
	$f[7] .= 'CSQ='.$cache->{$key};
	print VCFCACHE join("\t",@f)."\n";
    }
    else {
	print VCF4VEP "$line\n";
    }
}

close(VCF4VEP);
close(VCFCACHE);


$now = strftime("%F %T", localtime);
warn "I $now: $0 - finished parsing stdin and splitting it into $tmpDir files\n";


##########################################################################
# run VEP on $vcf4vep, producing $vcfFromVep
my $vcfFromVep = "$tmpDir/vcfFromVep.vcf.gz";

system("gunzip -c $vcf4vep | $vepCommand | gzip -c --fast > $vcfFromVep") ;

$now = strftime("%F %T", localtime);
warn "I $now: $0 - finished running VEP on the new variants\n";

##########################################################################
# merge $vcfFromVep and $vcfFromCache, printing resulting VCF to stdout;
# also fill $cacheUpdate with new CSQs from $vcfFromVep

open(VCFVEP,"gunzip -c $vcfFromVep |") || 
    die "E $0: cannot gunzip-open VCFVEP $vcfFromVep\n";
open (VCFCACHE, "gunzip -c $vcfFromCache |") || 
    die "E $0: cannot gunzip-open VCFCACHE $vcfFromCache\n";

# copy headers from VCFVEP
while (my $line = <VCFVEP>) {
    print $line;
    ($line =~ /^#CHROM/) && last;
}

# next data lines from each file
my $nextVep = <VCFVEP>;
my $nextCache = <VCFCACHE>;

# (chr,pos) from each next line (undef if no more data),
# for chr we use the number and replace X,Y,M with 23-25
my ($nextVepChr,$nextVepPos);
if ($nextVep && ($nextVep =~ /^chr(\w+)\t(\d+)\t/)) {
    ($nextVepChr,$nextVepPos)=($1,$2);
    if ($nextVepChr eq "X") {$nextVepChr = 23;}
    elsif ($nextVepChr eq "Y") {$nextVepChr = 24;}
    elsif ($nextVepChr eq "M") {$nextVepChr = 25;}
}
elsif ($nextVep) {
    die"E $0: vcfFromVep has a first data line but I can't parse it:\n$nextVep\n";
}

my ($nextCacheChr,$nextCachePos);
if ($nextCache && ($nextCache =~ /^chr(\w+)\t(\d+)\t/)) {
    ($nextCacheChr,$nextCachePos)=($1,$2);
    if ($nextCacheChr eq "X") {$nextCacheChr = 23;}
    elsif ($nextCacheChr eq "Y") {$nextCacheChr = 24;}
    elsif ($nextCacheChr eq "M") {$nextCacheChr = 25;}
}
elsif ($nextCache) {
    die"E $0: vcfFromCache has a first data line but I can't parse it:\n$nextCache\n";
}


# as long as there is data
while ($nextVep || $nextCache) {
    if ($nextVep) {
	# still have VEP data
	if ((! $nextCache) || ($nextVepChr < $nextCacheChr) ||
	    (($nextVepChr == $nextCacheChr) && ($nextVepPos <= $nextCachePos))) {
	    # and this VEP data must be printed now!
	    print $nextVep;

	    # also save for updating cache at the end
	    my @f = split(/\t/,$nextVep);
	    (@f >= 8) || die "E $0: nextVep line doesn't have >=8 columns:\n$nextVep";
	    # key: chrom:pos:ref:alt
	    my $key = "$f[0]:$f[1]:$f[3]:$f[4]";
	    ($f[7] =~ /CSQ=([^;]+)/) || die "E $0: cannot grab CSQ in nextVep line:\n$nextVep";
	    my $csq = $1;
	    $cacheUpdate->{$key} = $csq;

	    # read next VEP line
	    $nextVep = <VCFVEP>;
	    if ($nextVep && ($nextVep =~ /^chr(\w+)\t(\d+)\t/)) {
		($nextVepChr,$nextVepPos)=($1,$2);
		if ($nextVepChr eq "X") {$nextVepChr = 23;}
		elsif ($nextVepChr eq "Y") {$nextVepChr = 24;}
		elsif ($nextVepChr eq "M") {$nextVepChr = 25;}
	    }
	    elsif ($nextVep) {
		die"E $0: vcfFromVep has a data line but I can't parse it:\n$nextVep\n";
	    }
	    next;
	}
	# else print $nextCache and update it, but we do this also if ! $nextVep below
    }

    # no else because we want to use $nextCache as long as we didn't "next" above
    print $nextCache;
    $nextCache = <VCFCACHE>;
    if($nextCache && ($nextCache =~ /^chr(\w+)\t(\d+)\t/)) {
	($nextCacheChr,$nextCachePos)=($1,$2);
	if ($nextCacheChr eq "X") {$nextCacheChr = 23;}
	elsif ($nextCacheChr eq "Y") {$nextCacheChr = 24;}
	elsif ($nextCacheChr eq "M") {$nextCacheChr = 25;}
    }
    elsif ($nextCache) {
	die"E $0: vcfFromCache has a data line but I can't parse it:\n$nextCache\n";
    }
}

close(VCFVEP);
close(VCFCACHE);

##########################################################################
# save any new entries to $cacheFile: we lock $cacheFile and retrieve a fresh copy
# (which will include any updates made by other processes while this process was running),
# then add any new entries and store
# -> works correctly even if several jobs are running in parallel
open(CACHELOCK, ">>$cacheFile") ||
    die "E $0: I want to update cacheFile but I can't open it for locking\n";
flock(CACHELOCK, LOCK_EX) || die "E $0: cannot get lock on cachefile\n";

if (! -z $cacheFile) {
    $cache = &retrieve($cacheFile) ||
	die "E $0: I can't retrieve cache from fresh copy of cachefile $cacheFile\n";
}
foreach my $k (keys(%$cacheUpdate)) {
    if (defined($cache->{$k})) {
	# annotation for key $k has been added to $cacheFile while we were running
	# (or was already there if $debug), make sure it's consistent
	($cache->{$k} eq $cacheUpdate->{$k}) ||
	    warn "W $0: cacheFile entry for $k was added while we were running and disagrees with our entry, we keep the cached version:\n".
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


##########################################################################
# clean up
unlink($vcf4vep) || die "E $0: cannot unlink tmpfile vcf4vep $vcf4vep\n";
unlink($vcfFromCache) || die "E $0: cannot unlink tmpfile vcfFromCache $vcfFromCache\n";
unlink($vcfFromVep) || die "E $0: cannot unlink tmpfile vcfFromVep $vcfFromVep\n";
rmdir($tmpDir) || die "E $0: cannot rmdir tmpdir $tmpDir\n";


$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";
