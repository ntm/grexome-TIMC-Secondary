

package grexome_sec_config;
require Exporter;
our @ISA = qw(Exporter) ;
our @EXPORT_OK = qw(vep_bin vep_jobs genome_fasta vep_plugins vep_cacheFile gtexFile);

use strict;
use warnings;



########################
# runVEP.pl

# full path to vep
sub vep_bin {
    return("/usr/local/bin/vep");
}

# number of jobs that VEP can run in parallel (--fork)
sub vep_jobs {
    return(4);
}

# full path to human genome fasta, for VEP --hgvs
sub genome_fasta {
    # works on fauve or luxor, add other possibilities here if needed
    my @possibleGenomes = ("/data/HumanGenome/hs38DH.fa", "/home/nthierry/HumanGenome/hs38DH.fa");
    foreach my $genome (@possibleGenomes) {
	(-e $genome) && (return($genome));
    }
    # none of the @possibleGenomes exist
    die "E: genome_fasta in grexome_sec_config.pm can't find a genome fasta file\n";
}

# full path to the VEP cache file
sub vep_cacheFile {
    # works on fauve and luxor, file is actually currently on luxor
    my @possiblePaths = ("/data/nthierry/PierreRay/RunSecondaryAnalyses/",
			 "/home/nthierry/sshMounts/luxor/data/nthierry/PierreRay/RunSecondaryAnalyses/");
    foreach my $path (@possiblePaths) {
	(-d $path) && return("$path/VEP_cache");
    }
    die "E: vep_cacheFile can't find an existing dir in possiblePaths\n";
}

# full --plugin string with all VEP plugins we want to use and
# associated datafiles (with paths)
sub vep_plugins {
    # look for a data dir that works both on fauve and luxor
    my $dataDir = "/data/";
    (-d "$dataDir/nthierry/") && ($dataDir .= "nthierry/");

    my $plugins = "";

    # CADD - might not be needed if dbNSFP is good (it provides CADD 
    # among many other things), commenting out for now
    #my $caddPath = "$dataDir/CADD/";
    #$plugins .= " --plugin CADD,$caddPath/whole_genome_SNVs.tsv.gz,$caddPath/InDels.tsv.gz ";

    # dbNSFP
    my $dbNsfpPath = "$dataDir/dbNSFP/";
    # comma-separated list of fields to retrieve from dbNSFP, there are MANY
    # possibilities, check the README in $dbNsfpPath
    my $dbNsfpFields = "MutationTaster_pred,REVEL_rankscore,CADD_raw_rankscore";
    $plugins .= " --plugin dbNSFP,$dbNsfpPath/dbNSFP4.0a.gz,$dbNsfpFields ";

    # dbscSNV (splicing), data is with dbNSFP (same authors), specify 
    # assembly GRCh38 as second param because the plugin can't figure it out
    $plugins .= " --plugin dbscSNV,$dbNsfpPath/dbscSNV1.1_GRCh38.txt.gz,GRCh38 ";

    return($plugins);
}


########################
# addGTEX.pl

# full path to the GTEX datafile
sub gtexFile {
    # works on fauve and luxor
    my @possiblePaths = ("/home/nthierry/PierreRay/Grexome/SecondaryAnalyses/", 
			 "/home/nthierry/VariantCalling/GrexomeFauve/SecondaryAnalyses/");
    foreach my $path (@possiblePaths) {
	my $gtex = "$path/GTEX_Data/E-MTAB-5214-query-results.tpms.tsv";
	(-e $gtex) && (return($gtex)) ;
    }
    die "E: gtexFile in grexome_sec_config.pm can't find a GTEX datafile\n";
}

1;
