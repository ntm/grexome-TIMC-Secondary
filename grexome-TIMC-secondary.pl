#!/usr/bin/env perl

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
use File::Temp qw(tempdir);
use FindBin qw($RealBin);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## hard-coded stuff that shouldn't change much

# number of parallel jobs to run for the initial "bgzip -d",
# and then for steps steps 1 and 6.
# These defaults are good for us (dual Xeon 4114) but tweaking
# could improve performance depending on your hardware.
my $numJobsGunzip = 6;
my $numJobs1 = 20;
my $numJobs6 = 16;

# name of the single logfile that will be produced in $outDir
# (in debug mode this is replaced by one stepX.err logfile per step)
my $logfile = "grexomeTIMCsec.err";

#############################################
## options / params from the command-line

# metadata file with all samples
my $metadata;

# path+file holding known candidate genes
my $candidateGenes;

# input bgzipped multi-sample GVCF or VCF
my $inFile;

# outDir must not exist, it will be created and populated with
# subdirs (containing the pipeline results), logfiles, and copies
# of the provided metadata and candidateGenes.
my $outDir;

# path+file of the config file holding all install-specific params,
# defaults to the distribution-povided file that you can edit but
# you can also copy it elsewhere and customize it, then use --config
my $config = "$RealBin/grexomeTIMCsec_config.pm";

# debug: if true:
# - run each step one after the other (no pipes)
# - check result of previous step before starting next (die if prev step failed)
# - keep all intermediate files (no cleanup)
# - produce individual logfiles for each step
my $debug = '';

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "Parse a GVCF or VCF, run the complete grexome-TIMC secondary analysis pipeline, and produce results (+ logs and copies of the metadata) in the provided outDir (which must not exist).
Each step of the pipeline is a stand-alone self-documented script, this is just a wrapper.
Every install-specific param (eg paths to required data) should be in grexomeTIMCsec_config.pm.
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--metadata string [no default] : patient metadata xlsx file, with path
--candidateGenes string [no default] : known candidate genes in xlsx file, with path
--infile string [no default] : bgzipped multi-sample GVCF or VCF file to parse
--outdir string [no default] : subdir where resulting cohort files will be created, must not pre-exist
--config string [$config] : your customized copy (with path) of the distributed *config.pm
--debug : activate debug mode => slower, keeps all intermediate files, produce individual logfiles
--help : print this USAGE";

GetOptions ("metadata=s" => \$metadata,
	    "candidateGenes=s" => \$candidateGenes,
	    "infile=s" => \$inFile,
	    "outdir=s" => \$outDir,
	    "config=s" => \$config,
 	    "debug" => \$debug,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($metadata) || die "E $0: you must provide a metadata file\n";
(-f $metadata) || die "E $0: the supplied metadata file doesn't exist:\n$metadata\n";
($candidateGenes) || die "E $0: you must provide a candidateGenes file\n";
(-f $candidateGenes) || die "E $0: the supplied candidateGenes file $candidateGenes doesn't exist:\n$candidateGenes\n";

($inFile) || die "E $0: you must provide an input bgzipped (G)VCF file\n";
(-f $inFile) || die "E $0: the supplied infile doesn't exist\n";
($inFile =~ /\.gz$/) || die "E $0: the supplied infiledoesn't seem bgzipped\n";

# immediately import $config, so we die if file is broken
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n";
require($config);
grexomeTIMCsec_config->import();

($outDir) || die "E $0: you must provide an outDir\n";
(-e $outDir) && 
    die "E $0: outDir $outDir already exists, remove it or choose another name.\n";
mkdir($outDir) || die "E $0: cannot mkdir outDir $outDir\n";

# copy the provided metadata and candidateGenes files into $outDir
copy($metadata, $outDir) ||
    die "E $0: cannot copy metadata to outDir: $!\n";
copy($candidateGenes, $outDir) ||
    die "E $0: cannot copy candidateGenes to outDir: $!\n";
# use the copied versions in scripts (eg if original gets edited while analysis is running...)
$metadata = "$outDir/".basename($metadata);
$candidateGenes = "$outDir/".basename($candidateGenes);

#############################################

my $now = strftime("%F %T", localtime);
warn "I $0: $now - starting to run: ".join(" ", $0, @ARGV)."\n";

# all intermediate results / tmp working files are created in $tmpdir,
# a randomly-named subdir of &fastTmpPath() to avoid clashes
# $tmpdir is removed afterwards except in debug mode
my $tmpdir = tempdir(DIR => &fastTmpPath());

######################
# STEPS 1-6, piped into each other except in debug mode

# decompress infile and step 1
my $com = "bgzip -cd -@".$numJobsGunzip." $inFile | perl $RealBin/1_filterBadCalls.pl --metadata=$metadata --tmpdir=$tmpdir/FilterTmp/ --jobs $numJobs1 ";
if ($debug) {
    # specific logfile from step and save its output
    $com .= "2> $outDir/step1.err > $outDir/step1.out";
    system($com) && die "E $0: debug mode on, step1 failed: $?";
    # next step will read this step's output
    $com = "cat $outDir/step1.out ";
}

# step 2
$com .= " | perl $RealBin/2_sampleData2genotypes.pl ";
if ($debug) {
    $com .= "2> $outDir/step2.err > $outDir/step2.out";
    system($com) && die "E $0: debug mode on, step2 failed: $?";
    $com = "cat $outDir/step2.out ";
}

# step 3
$com .= " | perl $RealBin/3_runVEP.pl --cacheFile=".&vepCacheFile()." --genome=".&refGenome()." --dataDir=".&vepPluginDataPath()." --tmpDir=$tmpdir/runVepTmpDir/ ";
if ($debug) {
    $com .= "2> $outDir/step3.err > $outDir/step3.out";
    system($com) && die "E $0: debug mode on, step3 failed: $?";
    $com = "cat $outDir/step3.out ";
}

# step 4
$com .= " | perl $RealBin/4_vcf2tsv.pl ";
if ($debug) {
    $com .= "2> $outDir/step4.err > $outDir/step4.out";
    system($com) && die "E $0: debug mode on, step4 failed: $?";
    $com = "cat $outDir/step4.out ";
}

# step 5
$com .= " | perl $RealBin/5_addGTEX.pl --favoriteTissues=".&gtexFavoriteTissues()." --gtex=".&gtexDatafile($RealBin)." ";
if ($debug) {
    $com .= "2> $outDir/step5.err > $outDir/step5.out";
    system($com) && die "E $0: debug mode on, step5 failed: $?";
    $com = "cat $outDir/step5.out ";
}

# step 6
$com .= " | perl $RealBin/6_extractCohorts.pl --metadata=$metadata --candidateGenes=$candidateGenes --outDir=$tmpdir/Cohorts/ --tmpDir=$tmpdir/TmpExtract/ --jobs=$numJobs6 ";
if ($debug) {
    $com .= "2> $outDir/step6.err";
    system($com) && die "E $0: debug mode on, step6 failed: $?";
}
else {
    # after step6 we have several files (one per cohort) so no more piping,
    # run steps 1-6 with a single logfile
    $com = "( $com ) 2> $outDir/$logfile";
    system($com) && die "E $0: steps 1 to 6 seem to have failed? $?";
}

######################
  
# STEP 7: filter variants and reorder columns, clean up unfiltered CohortFiles
$com = "perl $RealBin/7_filterAndReorderAll.pl $tmpdir/Cohorts/ $tmpdir/Cohorts_Filtered/ ";
if ($debug) {
    $com .= "2> $outDir/step7.err";
}
else {
    $com .= "2>> $outDir/$logfile";
    # remove unfiltered results in non-debug mode
    $com .= " ; rm -r $tmpdir/Cohorts/";
}
system($com) && die "E $0: step7 failed: $?";

# STEP 8 - SAMPLES
$com = "perl $RealBin/8_extractSamples.pl $metadata $tmpdir/Cohorts_Filtered/ $outDir/Samples/ ".&coveragePath()." ";
if ($debug) {
    $com .= "2> $outDir/step8s.err";
}
else {
    $com .= "2>> $outDir/$logfile";
}
system($com) && die "E $0: step8-samples failed: $?";

# STEP 8 - TRANSCRIPTS , adding patientIDs
$com = "( perl $RealBin/8_extractTranscripts.pl $tmpdir/Cohorts_Filtered/ $tmpdir/Transcripts_noIDs/ ; ";
$com .= " perl $RealBin/8_addPatientIDs.pl $metadata $tmpdir/Transcripts_noIDs/ $outDir/Transcripts/ )";
if ($debug) {
    $com .= "2> $outDir/step8t.err";
}
else {
    $com .= "2>> $outDir/$logfile";
    $com .= " ; rm -r $tmpdir/Transcripts_noIDs/ ";
}
system($com) && die "E $0: step8-transcripts failed: $?";

# STEP 9 - FINAL COHORTFILES , require at least one HV or HET sample and add patientIDs
$com = "( perl $RealBin/9_requireUndiagnosed.pl $tmpdir/Cohorts_Filtered/ $tmpdir/Cohorts_FINAL/ ; ";
$com .= " perl $RealBin/8_addPatientIDs.pl $metadata $tmpdir/Cohorts_FINAL/ $outDir/Cohorts/ )";
if ($debug) {
    $com .= "2> $outDir/step9-finalCohorts.err";
}
else {
    $com .= "2>> $outDir/$logfile";
    $com .= " ; rm -r $tmpdir/Cohorts_Filtered/ $tmpdir/Cohorts_FINAL/ ";
}
system($com) && die "E $0: step9-finalCohorts failed: $?";

######################
# STEP 9 - SUBCOHORTS (can run after requireUndiagnosed and addPatientIDs)
#
# this step is not generic yet...
# The idea is to produce Cohorts and Transcripts files corresponding to subsets of
# samples, all affected with the same pathology. Typically the subsets are samples
# that were provided by collaborators, and this allows us to send them the results
# concerning their patients.

# key==path+file defining a subCohort, value==pathology
my %subCohorts = ("~/VariantCalling/GrexomeFauve/Grexome_Metadata/4-SubCohorts/subCohort_FV.txt" => "Azoo",
		  "~/VariantCalling/GrexomeFauve/Grexome_Metadata/4-SubCohorts/subCohort_London.txt" => "Azoo",
		  "~/VariantCalling/GrexomeFauve/Grexome_Metadata/4-SubCohorts/subCohort_AzooZouari.txt" => "Azoo");

mkdir("$outDir/SubCohorts") || die "E $0: cannot mkdir $outDir/SubCohorts\n";
$com = "";
foreach my $subC (keys(%subCohorts)) {
    my $patho = $subCohorts{$subC};
    # grab filename from $subC and remove leading "subCohort_" and trailing .txt
    my $outFileRoot = basename($subC);
    ($outFileRoot =~ s/^subCohort_//) || die "E $0: cannot remove leading subCohort_ from subCohortFile $outFileRoot\n";
    ($outFileRoot =~ s/\.txt$//) || die "E $0: cannot remove .txt from subCohortFile $outFileRoot\n";
    $outFileRoot = "$outDir/SubCohorts/$outFileRoot";
   
    $com .= "perl $RealBin/9_extractSubcohort.pl $subC < $outDir/Cohorts/$patho.final.patientIDs.csv > $outFileRoot.cohort.csv ; ";
    $com .= "perl $RealBin/9_extractSubcohort.pl $subC < $outDir/Transcripts/$patho.Transcripts.patientIDs.csv > $outFileRoot.transcripts.csv ; ";
}
system($com) && die "E $0: step9-subCohorts failed: $?";

######################
# all done, clean up tmpdir
($debug) || rmdir($tmpdir) || die "E $0: all done but can't rmdir(tmpdir), why?\n";

$now = strftime("%F %T", localtime);
warn "I $0: $now - ALL DONE, completed successfully!\n";
