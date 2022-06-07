# NTM
# 27/07/2020

# Define the paths and filenames that should be install-specific and
# are needed by grexome-TIMC-secondary.pl .
# Every hard-coded path/filename in the pipeline should be here.
# Also define some data-specific or behavioral config (eg "gtexFavoriteTissues").
# Other hard-coded stuff you may want to tweak to adjust the pipeline's
# behavior include:
# - some stuff in 3_runVEP.pl (eg VEP plugins and params to use)
# - %filterParams in 1_filterBadCalls.pl
# - @goodVeps in 4_vcf2tsv.pl
# - @newOrder in 7_reorderColumns.pl
# -...?


package grexomeTIMCsec_config;

use strict;
use warnings;
use Exporter;
our @ISA = ('Exporter');
# NOTE: this file can and should be copied somewhere and customized, 
# therefore we should never "use grexomeTIMCsec_config" but instead
# provide the customized *config.pm as an argument, see --config in
# grexome-TIMC-secondary.pl for an example.
our @EXPORT_OK = qw(refGenome vepCacheFile vepPluginDataPath fastTmpPath 
                    coveragePath gtexDatafile gtexFavoriteTissues subCohorts);


#################################################################
# files and paths needed by the pipeline 


# Return the reference human genome fasta file, with path
sub refGenome {
    # return the first file that exists, so this works on all our servers
    foreach my $genome ("/home/nthierry/HumanGenome/hs38DH.fa",
			"/data/HumanGenome/hs38DH.fa") {
	(-f $genome) && return($genome);
    }
    # if we get here no file was found...
    die "E: no refGenome found, you need to edit *config.pm";
}


# Return the name (with path) of a cachefile that 3_runVEP.pl will use for
# storing VEP annotations.
# If the file doesn't exist it will be created, however the path must
# pre-exist.
# The cachefile can be shared by all users of the pipeline, and this would
# provide the best speedups; however all users need read-write access to it.
# NOTE that this cachefile is only used internally by 3_runVEP.pl
# and has nothing to do with the VEP-provided cache.
sub vepCacheFile {
    # you'll need to customize the array of possible values for $cachedir,
    # we use the first dir that exists (works on all our servers)
    foreach my $cachedir ("/data/nthierry/PierreRay/RunSecondaryAnalyses/",
			  "/home/nthierry/sshMounts/luxor/data/nthierry/PierreRay/RunSecondaryAnalyses/") {
	if (-d $cachedir) {
	    # the name of the file itself, can probably stay as-is
	    my $cacheFile = "VEP_cache";
	    $cacheFile = "$cachedir/$cacheFile";
	    # if it exists it must be writable
	    (! -e $cacheFile) || (-w $cacheFile) ||
		die "E: vepCacheFile $cacheFile exists but you can't write to it, fix the permissions or use another cachedir";
	    return($cacheFile);
	}
    }
    # no dir was found
    die "E: no vepCachePath found, you need to edit *config.pm";
}

# Return a dir containing subdirs with the data required by 
# the VEP plugins we use (seach for "plugin" in 3_runVEP.pl).
sub vepPluginDataPath {
    # return first existing subdir
    foreach my $dir ("/data/nthierry/", "/data/") {
	(-d $dir) && return($dir);
    }
    # no dir found
    die "E: no vepPluginDataPath found, you need to edit *config.pm";
}


# Return a tmp dir with fast RW access, ideally a ramdisk, otherwise
# hopefully at least an SSD.
# Example: on linux you can make a 96GB ramdisk (assuming you have enough RAM)
# accessible in /mnt/RamDisk/ with:
### sudo mkdir /mnt/RamDisk
### sudo mount -t tmpfs -o size=96g myramdisk /mnt/RamDisk/
# You can make the mount automatic on boot by adding to /etc/fstab:
### tmpfs /mnt/RamDisk tmpfs size=96g 0 0
sub fastTmpPath {
    my $ramdisk = "/mnt/RamDisk/";
    (-d $ramdisk) && return($ramdisk);
    die "E: no fastTmpPath found, you need to edit *config.pm";
}


# return a dir containing the per-sample coverage statistics, as produced
# by 0_coverage.pl ; or "" of you don't want to use this feature.
# If non-empty this allows 8_extractSamples.pl to append coverage summary
# statistics to each sample's header line.
sub coveragePath {
    # specify $covDir="" to disable this feature
    my $covDir = "/data/nthierry/PierreRay/CoverageResults/Coverage_Rolling/";
    if (($covDir eq "") || (-d $covDir)) {
	return($covDir);
    }
    else {
	die "E: covDir was provided but doesn't exist, you need to edit *config.pm";
    }
}

#################################################################
# Some more hard-coded stuff that affects the behavior of the pipeline


# Return the GTEX datafile, with path.
# Probably shouldn't be touched except when we update our GTEX data.
# This sub takes as arg the path to the grexome-TIMC-secondary codebase,
# because the current GTEX datafile is included in the released codebase
# but this *config.pm file could have been copied to a different subdir
# for easier customization.
sub gtexDatafile {
    (@_ == 1) || die "E: grexDatafile needs one arg";
    my ($secPath) = @_;
    # actual file name should only change when we update our GTEX data
    my $gtex = "$secPath/GTEX_Data/E-MTAB-5214-query-results.tpms.tsv";
    (-f $gtex) && return($gtex);
    die "E: no gtexDatafile found! unexpected since it should be part of the git repo...";
}

# Return a comma-separated list of favorite tissues, these must match
# column headers from the gtexDatafile
sub gtexFavoriteTissues {
    # default: working on infertility here...
    my $favoriteTissues = "testis,ovary";
    return($favoriteTissues);
}

# Return a ref to a hash:
# key==path+file defining a subCohort (filename must start with subCohort_ and end with .txt),
# value==pathology (must match the pathologyID column of the pathologies metadata xlsx).
sub subCohorts {
    my %subCohorts = ("/home/nthierry/VariantCalling/GrexomeFauve/Grexome_Metadata/4-SubCohorts/subCohort_FV.txt" => "NOA",
		      "/home/nthierry/VariantCalling/GrexomeFauve/Grexome_Metadata/4-SubCohorts/subCohort_London.txt" => "NOA",
		      "/home/nthierry/VariantCalling/GrexomeFauve/Grexome_Metadata/4-SubCohorts/subCohort_AzooZouari.txt" => "NOA");
    return(\%subCohorts);
}


# module loaded ok
1;
