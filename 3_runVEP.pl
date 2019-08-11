#!/usr/bin/perl

# 24/03/2018
# NTM

# Take a single arg: a non-existing tmp dir.
# Read on stdin a VCF file, write to stdout a similar VCF file
# with added VEP annotations.
#
# 11/08/2019: adding a cache file (VEP is slow).
# CSQ are taken from the cachefile when the variant is in it,
# otherwise run VEP and add CSQ to cache.
# The cachefile must be manually removed if VEP and/or the VEP
# cache are updated (I am not implementing a timestamp check),
# it will then be rebuilt from scratch.


use strict;
use warnings;
# Storable for the cache
use Storable;

#use lib '/home/nthierry/PierreRay/Grexome/SecondaryAnalyses/';
# export PERL5LIB=... before calling script, it's more flexible
use  grexome_sec_config qw(vep_bin vep_jobs genome_fasta vep_plugins vep_cacheFile);


(@ARGV == 1) || die "runVep.pl needs a non-existing tmpdir as arg\n";
my ($tmpDir) = @ARGV;
(-e $tmpDir) && die "E: tmpDir $tmpDir exists, please rm -r $tmpDir or use another tmpDir as arg\n";
mkdir($tmpDir) || die "E: cannot create tmpDir $tmpDir\n";


##########################################################################
# parse VCF on stdin.
# - headers and data lines absent from cache are printed to 
#   $vcf4vep (gzipped --fast);
# - data lines found in cache are annotated using cache 
#   and printed to $vcfFromCache.
my $vcf4vep = "$tmpDir/vcf4vep.vcf.gz";
open(VCF4VEP, "| gzip -c --fast > $vcf4vep") ||
    die "cannot open gzip pipe to vcf4vep $vcf4vep\n";
my $vcfFromCache = "$tmpDir/vcfFromCache.vcf.gz";
open(VCFCACHE, "| gzip -c --fast > $vcfFromCache") ||
    die "cannot open gzip pipe to vcfFromCache $vcfFromCache\n";


my $cacheFile = &vep_cacheFile();
# $cache is a hashref. key=="chr:pos:ref:alt", value==CSQ
my $cache = {};
# grab previously cached data
if (-f $cacheFile) {
    $cache = &retrieve($cacheFile) ||
	die "E: cachefile $cacheFile exists but I can't retrieve hash from it.\n";
}


# header
while (my $line = <STDIN>) {
    chomp($line);
    if ($line =~ /^##INFO=<ID=CSQ/) {
	# make sure the INFO->CSQ fields from file and cache match
	if (defined $cache->{"INFOCSQ"}) {
	    my $cacheLine = $cache->{"INFOCSQ"};
	    ($cacheLine eq $line) ||
		die "cacheLine and INFO-CSQ line from VCF are different:\n$cacheLine\n$line\n".
		"if you really want to use this vcf from STDIN you need to remove the cachefile $cacheFile\n";
	}
	else {
	    $cache->{"INFOCSQ"} = $line;
	}
    }
    # in any case print header lines to VCF4VEP
    print VCF4VEP "$line\n";
    # line #CHROM is always last header line
    ($line =~ /^#CHROM/) && last;
}
# data
while (my $line = <STDIN>) {
    chomp($line);
    my @f = split(/\t/,$line,-1);
    (@f >= 8) || die "VCF line doesn't have >= 8 columns:\n$line\n";
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


##########################################################################
# run VEP on $vcf4vep, producing $vcfFromVep
my $vcfFromVep = "$tmpDir/vcfFromVep.vcf.gz";

my $vepCommand = &vep_bin() ;
$vepCommand .= " --offline --format vcf --vcf" ;
# cache to use: refseq, merged, or ensembl (default)
# $vepCommand .= " --merged" ;
$vepCommand .= " --force_overwrite --no_stats" ;
$vepCommand .= " --allele_number"; # for knowing which CSQ annotates which ALT
$vepCommand .= " --canonical --biotype --xref_refseq --flag_pick_allele_gene";
# instead of --everything we select relevant options (eg not --regulatory)
$vepCommand .= " --sift b --polyphen b --symbol --numbers --total_length" ;
$vepCommand .= " --gene_phenotype --af --af_1kg --af_esp --af_gnomad";
$vepCommand .= " --pubmed --variant_class --check_existing ";
# commenting out "--domains", it's a bit massive and in non-deterministic order
# and I don't think anyone looks at it
$vepCommand .= " --fasta ".&genome_fasta()." --hgvs";
# plugins:
$vepCommand .= &vep_plugins();
# --fork borks when vep_jobs==1
(&vep_jobs() > 1) && ($vepCommand .= " --fork ".&vep_jobs()) ;
# write output to stdout so we can pipe it to another program
$vepCommand .= " -o STDOUT" ;

system("gunzip -c $vcf4vep | $vepCommand | gzip -c --fast > $vcfFromVep") ;


##########################################################################
# merge $vcfFromVep and $vcfFromCache, printing resulting VCF to stdout;
# also update cache with new CSQs from $vcfFromVep

open(VCFVEP,"gunzip -c $vcfFromVep |") || 
    die "cannot gunzip-open VCFVEP $vcfFromVep\n";
open (VCFCACHE, "gunzip -c $vcfFromCache |") || 
    die "cannot gunzip-open VCFCACHE $vcfFromCache\n";

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
    die"vcfFromVep has a first data line but I can't parse it:\n$nextVep\n";
}

my ($nextCacheChr,$nextCachePos);
if ($nextCache && ($nextCache =~ /^chr(\w+)\t(\d+)\t/)) {
    ($nextCacheChr,$nextCachePos)=($1,$2);
    if ($nextCacheChr eq "X") {$nextCacheChr = 23;}
    elsif ($nextCacheChr eq "Y") {$nextCacheChr = 24;}
    elsif ($nextCacheChr eq "M") {$nextCacheChr = 25;}
}
elsif ($nextCache) {
    die"vcfFromCache has a first data line but I can't parse it:\n$nextCache\n";
}


# as long as there is data
while ($nextVep || $nextCache) {
    if ($nextVep) {
	# still have VEP data
	if ((! $nextCache) || ($nextVepChr < $nextCacheChr) ||
	    (($nextVepChr == $nextCacheChr) && ($nextVepPos <= $nextCachePos))) {
	    # and this VEP data must be printed now!
	    print $nextVep;

	    # also update cache
	    my @f = split(/\t/,$nextVep);
	    (@f >= 8) || die "nextVep line doesn't have >=8 columns:\n$nextVep";
	    # key: chrom:pos:ref:alt
	    my $key = "$f[0]:$f[1]:$f[3]:$f[4]";
	    ($f[7] =~ /CSQ=([^;]+)/) || die "cannot grab CSQ in nextVep line:\n$nextVep";
	    my $csq = $1;
	    $cache->{$key} = $csq;

	    # read next VEP line
	    $nextVep = <VCFVEP>;
	    if ($nextVep && ($nextVep =~ /^chr(\w+)\t(\d+)\t/)) {
		($nextVepChr,$nextVepPos)=($1,$2);
		if ($nextVepChr eq "X") {$nextVepChr = 23;}
		elsif ($nextVepChr eq "Y") {$nextVepChr = 24;}
		elsif ($nextVepChr eq "M") {$nextVepChr = 25;}
	    }
	    elsif ($nextVep) {
		die"vcfFromVep has a data line but I can't parse it:\n$nextVep\n";
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
	die"vcfFromCache has a data line but I can't parse it:\n$nextCache\n";
    }
}

close(VCFVEP);
close(VCFCACHE);

# save cache
&store($cache, $cacheFile) || 
    die "E: produced/updated cache but cannot store to cachefile $cacheFile\n";

# clean up
unlink($vcf4vep) || die "cannot unlink tmpfile vcf4vep $vcf4vep\n";
unlink($vcfFromCache) || die "cannot unlink tmpfile vcfFromCache $vcfFromCache\n";
unlink($vcfFromVep) || die "cannot unlink tmpfile vcfFromVep $vcfFromVep\n";
rmdir($tmpDir) || die "cannot rmdir tmpdir $tmpDir\n";
