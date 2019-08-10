#!/usr/bin/perl

# 24/03/2018
# NTM

# Read on stdin a VCF file, write to stdout a similar VCF file
# with added VEP annotations.

use strict;
use warnings;

#use lib '/home/nthierry/PierreRay/Grexome/SecondaryAnalyses/';
# export PERL5LIB=... before calling script, it's more flexible
use  grexome_sec_config qw(vep_bin vep_jobs genome_fasta vep_plugins);

(@ARGV == 0) || die "runVep.pl doesn't take args\n";



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
# trying to add --synonyms because I noticed some 'ENST\d+\.\d:c\.\d+N>' ,
# but actually most of the HGVS are fine... Need to try on a full file though
# UPDATE 07/08/2019: found out that the fasta index on luxor was KO, missing
# every chrom beyond chr4... remade it, keep searching for 'N>' in HGVS though...
#$vepCommand .= " --synonyms /home/nthierry/.vep/homo_sapiens/97_GRCh38/chr_synonyms.txt ";
# plugins:
$vepCommand .= &vep_plugins();
# --fork borks when vep_jobs==1
(&vep_jobs() > 1) && ($vepCommand .= " --fork ".&vep_jobs()) ;
# write output to stdout so we can pipe it to another program
$vepCommand .= " -o STDOUT" ;

system($vepCommand) ;
