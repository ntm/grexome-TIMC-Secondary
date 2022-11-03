#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2022
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

my $primDate = "220517";
my $secDate =  "220518";

my $cnvFile = "/data/septiera/InfertilityCohort_Results/Calling_results_ExomeDepth_220427_574pat/CNVResults_ExomeDepth_574samples_220428.vcf";
my $canon = "";
# set to "--canonical" to limit to canonical transcripts
# $canon = "--canonical";


($canon) || ($secDate .= "_AllTranscripts");
my $metaPath = "/home/nthierry/VariantCalling/GrexomeFauve/Grexome_Metadata/";
my $binPath = "/home/nthierry/Software/Grexome-TIMC/grexome-TIMC-Secondary/";

foreach my $caller ("Strelka", "GATK") {
    my $callerLow = lc($caller);
    my $com = "perl $binPath/grexome-TIMC-secondary.pl --samples=$metaPath/patient_summary.xlsx ";
    $com .= "--pathologies=$metaPath/3-CandidateGenes/pathologies.xlsx ";

    $com .= "--candidateGenes=";
    foreach my $pathos ("NOA-OA-OMD-POF-DSD", "preprm2", "MMAF-PCD-AST-Necro",
			"Globo-Macro-Headless-FF-PreImpA-OATS-Terato", "OG") {
	$com .= "$metaPath/3-CandidateGenes/candidateGenes_$pathos.xlsx,";
    }
    ($com =~ s/,$/ /) ||
	die "E adding pathos to com: cannot replace trailing comma with space in:\n$com\n";

    $com .= "--infile=../GVCFs_grexome/GVCFs_${caller}_Filtered_Merged/grexomes_${callerLow}_merged_$primDate.g.vcf.gz ";
    $com .= "--cnvs=$cnvFile ";
    $com .= "--outDir=SecondaryAnalyses_${secDate}_$caller ";
    $com .= "$canon ";
    #$com .= "--debugVep ";
    $com .= "2> grexomeTIMCsec_${secDate}_$caller.log ";
    # move log into caller results subdir
    $com .= "; mv grexomeTIMCsec_${secDate}_$caller.log SecondaryAnalyses_${secDate}_$caller/ ";
    # running for all callers in parallel:
    $com = "( $com ) &";
    system($com);
}
