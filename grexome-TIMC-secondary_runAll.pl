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

# Wrapper for running grexome-TIMC-Secondary.pl in parallel on 
# all variant-caller GVCFs.
# This script is very much install-specific and can't be re-used as-is,
# but it serves as a full-blown example since this is how we use the pipeline
# routinely: when we obtain new samples we simply update $primDate, $secDate and
# $cnvFile , then run the script.

use strict;
use warnings;

my $primDate = "250511";
my $secDate =  "250515";

my $cnvFile = "/home/nthierry/PierreRay_DATA/RunJACNEx/JACNEx_canonical_114/VCFs/CNVs_2025-05-12_15-54-14_NODUPs.vcf.gz";

# set to "--canonical" to limit to canonical transcripts
#my $canon = "";
my $canon = "--canonical";


($canon) || ($secDate .= "_AllTranscripts");
my $metaPath = "/home/nthierry/GrexomeZoufris/WES/";
my $binPath = "/home/nthierry/Software/Grexome-TIMC/grexome-TIMC-Secondary/";

foreach my $caller ("Strelka", "GATK") {
    my $callerLow = lc($caller);
    my $com = "perl $binPath/grexome-TIMC-secondary.pl --samples=$metaPath/patient_summary.xlsx ";
    $com .= "--pathologies=$metaPath/MetaData/pathologies.xlsx ";
    $com .= "--candidateGenes=$metaPath/MetaData/candidateGenes_allPhenotypes.xlsx ";
    $com .= "--infile=../GVCFs_grexome/GVCFs_${caller}_Filtered_Merged/grexomes_${callerLow}_merged_$primDate.g.vcf.gz ";
    ($cnvFile) && ($com .= "--cnvs=$cnvFile ");
    $com .= "--outDir=SecondaryAnalyses_${secDate}_$caller ";
    $com .= "$canon ";
    $com .= "2> grexomeTIMCsec_${secDate}_$caller.log ";
    # move log into caller results subdir
    $com .= "; mv grexomeTIMCsec_${secDate}_$caller.log SecondaryAnalyses_${secDate}_$caller/ ";
    # running for all callers in parallel:
    $com = "( $com ) &";
    system($com);
}
