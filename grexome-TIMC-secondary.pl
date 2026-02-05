#!/usr/bin/env perl


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


# NTM
# 27/07/2020


# This is a wrapper script for the grexome-TIMC secondary analysis
# pipeline.
# Args: see $USAGE.


use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Copy qw(copy);
use File::Basename qw(basename);
use File::Path qw(remove_tree);
use File::Temp qw(tempdir);
use FindBin qw($RealBin);

use lib "$RealBin";
use grexome_metaParse qw(parsePathologies parseSamples parseCandidateGenes);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## hard-coded stuff that shouldn't change much

# number of parallel jobs to run for the initial "bgzip -d",
# and then for steps 1, 4 and 9.
# These defaults are good for us (dual Xeon 4114) but tweaking
# could improve performance depending on your hardware.
my $numJobsGunzip = 8;
my $numJobs1 = 20;
my $numJobs4 = 20;
my $numJobs9 = 16;

# name (+path if needed) of bash binary, needed for -o pipefail
my $bash = "bash";
system("which $bash &> /dev/null") &&
    die "E $0: the bash executable $bash can't be found\n";

# name (+path if needed) of gzip-like binary, we like bgzip (multi-threaded)
my $bgzip = "bgzip";
system("which $bgzip &> /dev/null") &&
    die "E $0: the bgzip executable $bgzip can't be found\n";
# full command. If you don't have bgzip you could use gzip but without --threads
$bgzip = "$bgzip -cd --threads $numJobsGunzip ";


#############################################
## options / params from the command-line

# metadata file with all samples
my $samples;

# pathologies metadata file
my $pathologies;

# comma-separated list of path+files holding known candidate genes
my $candidateGenes;

# input bgzipped multi-sample GVCF or VCF
my $inFile;

# input multi-sample VCF containing CNV calls, possibly (b)gzipped
my $cnvs;

# path+file with sampleIDs that belong to a sub-cohort (specific Cohort
# and Transcripts files will be produced for them), filename == [subcohortName].txt
# with one sampleID per line
my $subcohortFile;

# outDir must not exist, it will be created and populated with
# subdirs (containing the pipeline results), logfiles (in debug mode),
# and copies of all provided metadata files.
my $outDir;

# species, for 04_runVEP.pl
my $species = "homo_sapiens";

# path+file of the config file holding all install-specific params,
# defaults to the distribution-povided file that you can edit but
# you can also copy it elsewhere and customize it, then use --config
my $config = "$RealBin/grexomeTIMCsec_config.pm";

# if $canon is true, results are restricted to CANONICAL transcripts.
# Otherwise you can filter the CSVs using the CANONICAL column,
# and we also still produce canonical-only Cohorts files (because
# the unrestricted Cohorts files can be too heavy for oocalc/excel)
my $canon = '';

# debug: if true:
# - run each step one after the other (no pipes)
# - check result of previous step before starting next (die if prev step failed)
# - keep all intermediate files (no cleanup)
# - produce individual logfiles for each step
my $debug = '';

# $debugVep: activate --debug specifically for 04_runVep.pl
my $debugVep = '';

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "Parse a GVCF or VCF, run the complete grexome-TIMC secondary analysis pipeline, 
and produce results (+ logs and copies of the metadata) in the provided outDir (which must not exist).
Each step of the pipeline is a stand-alone self-documented script, this is just a wrapper.
Every install-specific param (eg paths to shared data) should be in grexomeTIMCsec_config.pm.
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--samples : samples metadata xlsx file, with path
--pathologies : [optional] pathologies metadata xlsx file, with path
--candidateGenes : [optional] known candidate genes in xlsx files, comma-separated, with paths
--infile : bgzipped multi-sample GVCF or VCF file to parse
--species string [default $species]: species, as expected by VEP (eg mus_musculus)
--cnvs : [optional] multi-sample VCF file containing CNV calls, possibly (b)gzipped, conventions are:
                ALT is <DEL> or <DUP>, INFO contains END=, FORMAT must start with GT:GQ:FR:BPR and possibly :BP
--subcohort : [optional] txt file holding sampleIDs (one per line), specific Cohort and Transcripts
            files will be produced for them, filename is used as subcohort name
--outdir : subdir where results will be created, must not pre-exist
--config [defaults to grexomeTIMCsec_config.pm alongside this script] : your customized copy (with path) of the distributed *config.pm
--canonical : restrict results to canonical transcripts
--debug : activate debug mode => slower, keeps all intermediate files, produce individual logfiles
--debugVep : activate debug mode specifically for VEP => compare VEP annotations between runs
--help : print this USAGE";

GetOptions ("samples=s" => \$samples,
            "pathologies=s" => \$pathologies,
            "candidateGenes=s" => \$candidateGenes,
            "infile=s" => \$inFile,
            "species=s" => \$species,
            "cnvs=s" => \$cnvs,
            "subcohort=s" => \$subcohortFile,
            "outdir=s" => \$outDir,
            "config=s" => \$config,
            "canonical" => \$canon,
            "debug" => \$debug,
            "debugVep" => \$debugVep,
            "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($samples) || die "E $0: you must provide a samples file\n\n$USAGE\n";
(-f $samples) || die "E $0: the supplied samples file doesn't exist:\n$samples\n";

($inFile) || die "E $0: you must provide an input bgzipped (G)VCF file\n";
(-f $inFile) || die "E $0: the supplied infile doesn't exist\n";
($inFile =~ /\.gz$/) || die "E $0: the supplied infile doesn't seem bgzipped\n";

($species =~ /^\w+$/) || die "E $0: VEP species must be alphanumeric+underscores\n";

if ($cnvs) {
    (-f $cnvs) || die "E $0: the supplied cnvs file doesn't exist\n";
}

# immediately import $config, so we die if file is broken
# if $config doesn't have a path component, prepend ./ to avoid loading the dist version
# (in case the dist was copied into current dir and customized but not renamed)
($config =~ m~/~) || ($config = "./$config");
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n";
require($config);
grexomeTIMCsec_config->import(qw(refGenome vepCacheFile vepPluginDataPath fastTmpPath), 
                              qw(coveragePath gtexDatafile gtexFavoriteTissues));

($outDir) || die "E $0: you must provide an outDir\n";
(-e $outDir) && 
    die "E $0: outDir $outDir already exists, remove it or choose another name.\n";
mkdir($outDir) || die "E $0: cannot mkdir outDir $outDir\n";


# copy all provided metadata files into $outDir
copy($samples, $outDir) ||
    die "E $0: cannot copy samples metadata to outDir: $!\n";
# use the copied versions in scripts (eg if original gets edited while analysis is running...)
$samples = "$outDir/".basename($samples);

if ($pathologies) {
    (-f $pathologies) || die "E $0: the supplied pathologies file doesn't exist\n";
    copy($pathologies, $outDir) ||
        die "E $0: cannot copy pathologies metadata to outDir: $!\n";
    $pathologies = "$outDir/".basename($pathologies);
}

my @candNew = ();
if ($candidateGenes) {
    foreach my $candFile (split(/,/, $candidateGenes)) {
        (-f $candFile) ||
            die "E $0: the supplied candidateGenes file $candFile doesn't exist\n";
        copy($candFile, $outDir) ||
            die "E $0: cannot copy candidateGenes file $candFile to outDir: $!\n";
        # use the copies in script
        push(@candNew, "$outDir/".basename($candFile));
    }
    $candidateGenes = join(',', @candNew);
}

if ($subcohortFile) {
    (-f $subcohortFile) || die "E $0: the supplied subcohort file doesn't exist\n";
    copy($subcohortFile, $outDir) ||
        die "E $0: cannot copy subcohort file to outDir: $!\n";
    # use copy from now on
    $subcohortFile = "$outDir/".basename($subcohortFile);
}


#############################################
my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";

# sanity-check all provided metadata files: just parse and ignore the results,
# we'll die if anything is wrong, but call parseCandidates in verbose mode
if ($pathologies) {
    &parsePathologies($pathologies);
    &parseSamples($samples, $pathologies);
    if ($candidateGenes) {
        &parseCandidateGenes($candidateGenes, $samples, 1, $pathologies);
    }
}
else {
    &parseSamples($samples);
    if ($candidateGenes) {
        &parseCandidateGenes($candidateGenes, $samples, 1);
    }
}

# pre-process subcohortFile (if provided):
# - make sure every subCohort sampleID exists
# - produce one subcohortFile per pathologyID, store results in hash with
#   key == subcFile, value == pathoID
my %subcFile2patho = ();
if ($subcohortFile) {
    my ($sample2pathoR) = &parseSamples($samples);

    # parse subcohortFile and separate sampleIDs by patho, save in hash:
    # key == patho, value == \n-separated list of samples
    my %patho2subcSamps = ();

    open(SUBC, "$subcohortFile") ||
        die "E $0: cannot open subcohortFile $subcohortFile for reading\n";
    while (my $line = <SUBC>) {
        chomp($line);
        # skip blank lines
        ($line =~ /^\s*$/) && next;
        # remove leading or trailing blank chars and make sure we have a reasonable ID
        # (no whitespace, no ( or [)
        ($line =~ /^\s*(\S+)\s*$/) ||
            die "E $0: cannot find reasonable sampleID in line $line from subcohort file\n";
        my $samp = $1;
        ($samp =~ /[(\[]/) && 
            die "E $0: sampleID $samp from subcohort file contains ( or [, illegal\n";
        (defined $sample2pathoR->{$samp}) ||
            die "E $0: sampleID $samp from subcohort file doesn't exist in samples metadata\n";
        my $patho = $sample2pathoR->{$samp};
        (defined $patho2subcSamps{$patho}) || ($patho2subcSamps{$patho} = "");
        $patho2subcSamps{$patho} .= "$samp\n";
    }
    close(SUBC);
    
    foreach my $patho (keys(%patho2subcSamps)) {
        # create per-patho subcohort files alongside the copied $subcohortFile
        my $subcThisPatho = $subcohortFile;
        ($subcThisPatho =~ s/\.txt$/$patho.txt/) ||
            die "E $0: subcohort filename must end with .txt";
        (-e $subcThisPatho) && die "E $0: subcThispatho $subcThisPatho exists, impossible!\n";
        open(SUBC_OUT, "> $subcThisPatho") ||
            die "E $0: cannot open subcThisPatho for $patho: $!";
        print SUBC_OUT $patho2subcSamps{$patho};
        close(SUBC_OUT);
        $subcFile2patho{$subcThisPatho} = $patho;
    }
}


# sanity-check: any sample listed in $samplesFile MUST be present in $inFile,
# otherwise we can get errors in later steps (eg step 7).
# We also need the number of samples in $samples, to set $min_hr (for filtering).
my $numSamples = 0;
{
    # grab the list of samples present in $inFile:
    my @samplesFromInFile = split(/\s+/, `zgrep --max-count=1 '#CHROM' $inFile`);
    # first 9 columns aren't sampleIDs
    splice(@samplesFromInFile,0,9);
    # key == sampleID present in $inFile, value==1
    my %samplesFromIn = ();
    foreach my $sample (@samplesFromInFile) {
        $samplesFromIn{$sample} = 1;
    }

    # grab the samples listed in $samplesFile: just use the first hashref from
    # parseSamples, keys are samples
    my ($samplesFromMetadataR) = &parseSamples($samples);
    $numSamples = scalar(keys(%$samplesFromMetadataR));
    foreach my $s (sort keys(%$samplesFromMetadataR)) {
        (defined $samplesFromIn{$s}) ||
            die "E $0: every sample defined in samples file $samples MUST be present in inFile $inFile, sample $s isn't.\n";
    }
}

# variant-caller string/name, will be added to all final filenames
my $caller;
# try to find the caller name in $inFile
foreach my $c ("Strelka", "GATK", "ElPrep", "DeepVariant") {
    if ($inFile =~ /$c/i) {
        $caller = $c;
        warn "I $0: variant-caller id $caller will be appended to all filenames\n";
        last;
    }
}
($caller) ||
    warn "W $0: variant-caller could not be auto-detected, filenames won't be tagged\n";


#############################################

# all intermediate results / tmp working files are created in $tmpdir,
# a randomly-named subdir of &fastTmpPath() to avoid clashes
# $tmpdir is removed afterwards except in debug mode
my $tmpdir = tempdir(DIR => &fastTmpPath());

######################
# STEPS 1-9, piped into each other except in debug mode

# decompress infile and step 1
my $com = "$bgzip $inFile | perl $RealBin/01_filterBadCalls.pl --samplesFile=$samples --tmpdir=$tmpdir/FilterTmp/ --jobs $numJobs1 ";
if ($debug) {
    # specific logfile from step and save its output
    $com .= "2> $outDir/step1.err > $outDir/step1.out";
    system($com) && die "E $0: debug mode on, step1 failed, examine step1.err\n";
    # next step will read this step's output
    $com = "cat $outDir/step1.out ";
}

# step 2 if we have a CNVs file
if ($cnvs) {
    $com .= " | perl $RealBin/02_collateVCFs.pl --vcf=$cnvs ";
    if ($debug) {
        $com .= "2> $outDir/step2.err > $outDir/step2.out";
        system($com) && die "E $0: debug mode on, step2 failed, examine step2.err\n";
        $com = "cat $outDir/step2.out ";
    }
}
else {
    warn "I $0: CNVs not provided, step 02-collateVCFs skipped\n";
}

# step 3
$com .= " | perl $RealBin/03_sampleData2genotypes.pl ";
if ($debug) {
    $com .= "2> $outDir/step3.err > $outDir/step3.out";
    system($com) && die "E $0: debug mode on, step3 failed, examine step3.err\n";
    $com = "cat $outDir/step3.out ";
}

# step 4
$com .= " | perl $RealBin/04_runVEP.pl --cacheFile=".&vepCacheFile()." --genome=".&refGenome().
    " --dataDir=".&vepPluginDataPath()." --tmpDir=$tmpdir/runVepTmpDir/ --species=$species --jobs=$numJobs4";
($debugVep) && ($com .= "--debug ");
if ($debug) {
    $com .= "2> $outDir/step4.err > $outDir/step4.out";
    system($com) && die "E $0: debug mode on, step4 failed, examine step4.err\n";
    $com = "cat $outDir/step4.out ";
}

# step 5
$com .= " | perl $RealBin/05_vcf2tsv.pl ";
if ($debug) {
    $com .= "2> $outDir/step5.err > $outDir/step5.out";
    system($com) && die "E $0: debug mode on, step5 failed, examine step5.err\n";
    $com = "cat $outDir/step5.out ";
}

# step 6
$com .= " | perl $RealBin/06_checkCandidatesExist.pl --samples=$samples ";
($candidateGenes) && ($com .= "--candidateGenes=$candidateGenes ");
if ($debug) {
    $com .= "2> $outDir/step6.err > $outDir/step6.out";
    system($com) && die "E $0: debug mode on, step6 failed, examine step6.err\n";
    $com = "cat $outDir/step6.out ";
}

# step 7 - immediately filter variants on IMPACT, BIOTYPE, CANONICAL and AFs (human-only)
$com .= " | perl $RealBin/07_filterVariants.pl --logtime --no_mod --no_pseudo --no_nmd ";
($canon) && ($com .= "--canonical ");
if (($species eq 'homo_sapiens') || ($species eq 'human')) {
    # global freq <= 1% in each dataset, and per-population freq <= 5%
    $com .= "--max_af_global 0.01 --max_af_perPop 0.05 ";
}
if ($debug) {
    $com .= "2> $outDir/step7.err > $outDir/step7.out";
    system($com) && die "E $0: debug mode on, step7 failed, examine step7.err\n";
    $com = "cat $outDir/step7.out ";
}

# step 8
$com .= " | perl $RealBin/08_addGTEX.pl --favoriteTissues=".&gtexFavoriteTissues()." --gtex=".&gtexDatafile($RealBin)." ";
if ($debug) {
    $com .= "2> $outDir/step8.err > $outDir/step8.out";
    system($com) && die "E $0: debug mode on, step8 failed, examine step8.err\n";
    $com = "cat $outDir/step8.out ";
}

# step 9
$com .= " | perl $RealBin/09_extractCohorts.pl --samples=$samples ";
($pathologies) && ($com .= "--pathologies=$pathologies ");
($candidateGenes) && ($com .= "--candidateGenes=$candidateGenes ");
$com .= "--outDir=$tmpdir/Cohorts/ --tmpDir=$tmpdir/TmpExtract/ --jobs=$numJobs9 ";
if ($debug) {
    $com .= "2> $outDir/step9.err";
    system($com) && die "E $0: debug mode on, step9 failed, examine step9.err\n";
}
else {
    # after step9 we have several files (one per cohort) so no more piping,
    # run steps 1-9, all logs go to stderr, fail if any component of the pipe fails
    $com =~ s/"/\\"/g;
    $com = "$bash -o pipefail -c \" $com \"";
    system($com) && die "E $0: steps 1 to 9 have failed, check for previous errors in this logfile\n";
}

######################
# subsequent steps work on the individual CohortFiles

# STEP 10 - TRANSCRIPTS , BEFORE FILTERING on max_ctrl*, adding patientIDs
$com = "perl $RealBin/11_extractTranscripts.pl --indir $tmpdir/Cohorts/ --outdir $tmpdir/Transcripts_noIDs/ ";
($pathologies) && ($com .= "--pathologies=$pathologies ");
($debug) && ($com .= "2> $outDir/step10t.err");
system($com) && die "E $0: step10-transcripts failed\n";

$com = "perl $RealBin/12_addPatientIDs.pl $samples $tmpdir/Transcripts_noIDs/ $outDir/Transcripts/ ";
($debug) && ($com .= "2>> $outDir/step10t.err");
system($com) && die "E $0: step10-transcripts-addPatientIDs failed\n";
(! $debug) && remove_tree("$tmpdir/Transcripts_noIDs/");


# STEP 10b - filter variants on COUNTs and reorder columns, clean up unfiltered CohortFiles
$com = "perl $RealBin/10_filterAndReorderAll.pl --indir $tmpdir/Cohorts/ --outdir $tmpdir/Cohorts_Filtered/ ";
$com .= "--reorder --favoriteTissues=".&gtexFavoriteTissues()." ";
# set min_hr to 20% of $numSamples
my $min_hr = int($numSamples * 0.2);
$com .= "--min_hr=$min_hr ";
# hard-coded max_ctrl_hv:
my $max_ctrl_hv = 3;
# set max_ctrl_het to 2% of $numSamples (but always at least 2 * $max_ctrl_hv)
my $max_ctrl_het = int($numSamples * 0.02);
($max_ctrl_het < 2 * $max_ctrl_hv) && ($max_ctrl_het = 2 * $max_ctrl_hv);
$com .= "--max_ctrl_hv=$max_ctrl_hv --max_ctrl_het=$max_ctrl_het ";
($debug) && ($com .= "2> $outDir/step10f.err");
system($com) && die "E $0: step10-filter failed\n";
# remove unfiltered results in non-debug mode
(! $debug) && remove_tree("$tmpdir/Cohorts/");


# STEP 11 - SAMPLES
$com = "perl $RealBin/11_extractSamples.pl --samples $samples ";
$com .= "--indir $tmpdir/Cohorts_Filtered/ --outdir $outDir/Samples/ ";
(&coveragePath()) && ($com .= "--covdir ".&coveragePath());
($debug) && ($com .= " 2> $outDir/step11s.err");
system($com) && die "E $0: step11-samples failed\n";


# STEP 11 - FINAL COHORTFILES , require at least one HV or HET sample and add patientIDs
$com = "perl $RealBin/11_requireUndiagnosed.pl $tmpdir/Cohorts_Filtered/ $tmpdir/Cohorts_FINAL/ ";
($debug) && ($com .= "2> $outDir/step11-finalCohorts.err");
system($com) && die "E $0: step11-finalCohorts failed\n";

$com = "perl $RealBin/12_addPatientIDs.pl $samples $tmpdir/Cohorts_FINAL/ $outDir/Cohorts/ ";
($debug) && ($com .= "2>> $outDir/step11-finalCohorts.err");
system($com) && die "E $0: step11-finalCohorts-addPatientIDs failed\n";

(! $debug) && remove_tree("$tmpdir/Cohorts_Filtered/", "$tmpdir/Cohorts_FINAL/");


# STEP 12 - COHORTFILES_CANONICAL , if called without --canon we still
# produce CohortFiles restricted to canonical transcripts (in case
# the AllTranscripts CohortFiles are too heavy for oocalc/excel)
if (! $canon) {
    $com = "perl $RealBin/10_filterAndReorderAll.pl --indir $outDir/Cohorts/ --outdir $outDir/Cohorts_Canonical/ --canon ";
    ($debug) && ($com .= "2> $outDir/step12-cohortsCanonical.err");
    system($com) && die "E $0: step12-cohortsCanonical failed\n";
}


######################
# STEP13 - rename and clean up

# APPENDVC: append $caller (if it was auto-detected) to all final filenames
if ($caller) {
    open (FILES, "find $outDir/ -name \'*csv\' |") ||
        die "E $0: step13-appendVC cannot find final csv files with find\n";
    while (my $f = <FILES>) {
        chomp($f);
        my $new = $f; 
        ($new =~ s/\.csv$/.$caller.csv/) || 
            die "E $0: step13-appendVC cannot add $caller as suffix to $new\n";
        (-e $new) && 
            die "E $0: step13-appendVC want to rename to new $new but it already exists?!\n";
        rename($f,$new) ||
            die "E $0: step13-appendVC cannot rename $f $new\n";
    }
    close(FILES);
}

# REMOVEEMPTY: remove files with no data lines (eg if infile concerned only some samples)
open (FILES, "find $outDir/ -name \'*csv\' |") ||
    die "E $0: step13-removeEmpty cannot find final csv files with find\n";
while (my $f = <FILES>) {
    chomp($f);
    my $wc = `head -n 2 $f | wc -l`;
    # there's always 1 header line
    ($wc > 1) || unlink($f) ||
        die "E $0: step13-removeEmpty cannot unlink $f: $!\n";
}
close(FILES);


######################
# STEP 14 - SUBCOHORT, must run after APPENDVC and REMOVEEMPTY
#
# Only runs if called with --subcohort :
# the idea is to produce Cohorts and Transcripts files corresponding to subsets of
# samples. Typically the subsets are samples that were provided by collaborators,
# and this allows us to send them the results concerning their patients.
if ($subcohortFile) {
    # 13_extractSubcohort.pl gets called many times but we prefer a single pair of log
    # messages -> warn here, pretending to be 13_extractSubcohort.pl
    $now = strftime("%F %T", localtime);
    warn "I $now: 13_extractSubcohort.pl - starting to run\n";

    mkdir("$outDir/SubCohort") || die "E $0: cannot mkdir $outDir/SubCohort\n";
    mkdir("$outDir/SubCohort/Samples") || die "E $0: cannot mkdir $outDir/SubCohort/Samples\n";
    # output filenames: input filename prepended with $subcName
    my $subcName = basename($subcohortFile);
    ($subcName =~ s/.txt$//); # sanity already checked
    my $outFileRoot = "$outDir/SubCohort/$subcName";

    foreach my $subcFile (keys(%subcFile2patho)) {
        my $patho = $subcFile2patho{$subcFile};
        my $com = "( perl $RealBin/13_extractSubcohort.pl $subcFile ";
        $com .= "< $outDir/Cohorts/$patho.final.patientIDs.csv > $outFileRoot.$patho.cohort.csv ) && ";
        if (! $canon) {
            $com .= "( perl $RealBin/13_extractSubcohort.pl $subcFile ";
            $com .= "< $outDir/Cohorts_Canonical/$patho.final.patientIDs.canon.csv > $outFileRoot.$patho.cohort.canon.csv ) && ";
        }
        $com .= "( perl $RealBin/13_extractSubcohort.pl $subcFile ";
        $com .= "< $outDir/Transcripts/$patho.Transcripts.patientIDs.csv > $outFileRoot.$patho.transcripts.csv ) ";
        ($debug) && ($com .= "2>> $outDir/step14-subCohort.err ");
        system($com) && die "E $0: step14-subCohort failed\n";

        # also symlink Samples files
        open(SUBC, "$subcFile") || die "E $0: cannot open subcFile $subcFile for reading\n";
        while (my $samp = <SUBC>) {
            chomp($samp);
            my @infile = glob("$outDir/Samples/$patho.$samp.*");
            (@infile == 1) || die "E $0: extractSubcohort.pl finds several Samples files for $samp\n";
            my $sampFile = basename($infile[0]);
            symlink("../../Samples/$sampFile", "$outDir/SubCohort/Samples/$sampFile") ||
                die "E $0: cannot symlink $sampFile for subcohort\n";
        }
        close(SUBC);
    }
    
    $now = strftime("%F %T", localtime);
    warn "I $now: 13_extractSubcohort.pl - ALL DONE, completed successfully!\n";
}
else {
    warn "I $0: no provided subcohort, step14-subCohort skipped\n";
}


######################
# STEP15 - QC
#
# QC_CHECKCAUSAL: report hits of known causal genes by (severe) variants
# QC report will be printed to $qc_causal file
my $qc_causal = "$outDir/qc_causal.csv";

$com = "perl $RealBin/14_qc_checkCausal.pl --samplesFile=$samples --indir=$outDir/Samples/ ";
$com .= "> $qc_causal ";
($debug) && ($com .= "2> $outDir/step14-qc_causal.err");
system($com) && die "E $0: step15-qc_causal failed\n";


######################
# all done, clean up tmpdir
($debug) || rmdir($tmpdir) || die "E $0: all done but can't rmdir(tmpdir), why?\n";

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";
