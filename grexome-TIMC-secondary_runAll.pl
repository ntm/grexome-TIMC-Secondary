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


# 26/08/2021
# NTM

# Wrapper for running grexome-TIMC-Secondary.pl in parallel on all variant-caller GVCFs.
# This script is install-specific and can't be re-used as-is, but it's a great
# full-blown example since this is how we use grexome-TIMC-Secondary routinely:
# When new samples are received:
# -> update GVCFs (with grexome-TIMC-Primary, https://github.com/ntm/grexome-TIMC-Primary )
# -> update CNV calls (with JACNEx, https://github.com/ntm/JACNEx )
# -> run this script, it auto-finds the most recent GVCFs ($infile) and CNVs ($cnvFile) and
#    produces analysis-ready files in a date-stamped ($secDate) outDir.

use strict;
use warnings;

my $secDate =  `/bin/date +%y%m%d`;
chomp($secDate);

# set to "--canonical" to limit to canonical transcripts
#my $canon = "";
my $canon = "--canonical";


($canon) || ($secDate .= "_AllTranscripts");
my $metaPath = "/home/nthierry/GrexomeZoufris/WES/";
my $binPath = "/home/nthierry/Software/Grexome-TIMC/grexome-TIMC-Secondary/";

foreach my $caller ("Strelka", "GATK") {
#foreach my $caller ("Strelka") {
#foreach my $caller ("GATK") {
    my $callerLow = lc($caller);
    my $com = "perl $binPath/grexome-TIMC-secondary.pl --samples=$metaPath/patient_summary.xlsx ";
    $com .= "--pathologies=$metaPath/MetaData/pathologies.xlsx ";
    $com .= "--candidateGenes=$metaPath/MetaData/candidateGenes_allPhenotypes.xlsx ";
#    $com .= "--candidateGenes=$metaPath/MetaData/PIWI_genes.xlsx ";
    my $infile = `/bin/ls -1rt ../GVCFs_grexome/GVCFs_${caller}_Filtered_Merged/grexomes_${callerLow}_merged_*.g.vcf.gz | tail -n 1`;
    chomp($infile);
    my $cnvFile = `/bin/ls -1rt ../RunJACNEx/JACNEx_canonical_*/VCFs/CNVs_*_NODUPs.vcf.gz | tail -n 1`;
    chomp($cnvFile);
    $com .= "--infile=$infile --cnvs=$cnvFile ";
    $com .= "--outDir=SecondaryAnalyses_${secDate}_$caller ";
    $com .= "$canon ";
    # compare VEP annotations for the same variants from one VEP run to another
    # $com .= "--debugVep ";
    $com .= "2> grexomeTIMCsec_${secDate}_$caller.log ";
    # move log into caller results subdir
    $com .= "; mv grexomeTIMCsec_${secDate}_$caller.log SecondaryAnalyses_${secDate}_$caller/ ";
    # running for all callers in parallel:
    $com = "( $com ) &";
    system($com);
}
