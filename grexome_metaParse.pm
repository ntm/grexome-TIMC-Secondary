
############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2025
#
# This file is part of grexome-TIMC-Primary and grexome-TIMC-Secondary, written by
# Nicolas Thierry-Mieg (CNRS, France) Nicolas.Thierry-Mieg@univ-grenoble-alpes.fr
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
# 01/04/2021

# Subs for parsing the grexome-TIMC metadata files.
# All subs die with an explicit error message if files are not as expected.
# Example metadata files are provided in the Documentation/ subdir.
# This module is shared by the grexome-TIMC-Primary and
# grexome-TIMC-Secondary pipelines.


package grexome_metaParse;

use strict;
use warnings;
use Spreadsheet::XLSX;
use Exporter;
our @ISA = ('Exporter');
our @EXPORT_OK = qw(parsePathologies parseSamples parseCandidateGenes);


#################################################################

# Parse pathologies metadata XLSX file: required columns are pathologyID,
# no_files, is_a and "compatibility groups" (can be in any order but they MUST exist).
#
# Return 2 hashrefs: ($ancestorsR, $compatR), keys are a pathologyID,
# values are hashrefs whose keys are the pathoIDs defined by:
# - ancestorsR: all no_files==0 ancestors (ontology-wise) of this patho;
# - compatR: all no_files==0 pathos that samples of this patho are "compatible" with;
# all values are 1.
#
# The pathologyIDs are structured as a DAG, is_a holds the (possibly empty)
# comma-separated list of parent pathologyIDs.
# Each pathologyID can belong to zero or more compatibility groups: members of
# a group may share common causal genes, so they are not used as negative controls
# for one another.
# A pathologyID with no_files==1 will never appear in the returned hashrefs,
# and must not appear in the samples metadata XLSX; such no_files==1 pathoIDs
# are only useful for human-friendly structuring of terms.
# A sample belonging to pathoID implicitly also belongs to all ancestors of pathoID,
# and is compatible with all descendants of pathoID (in addition to any pathologyID
# belonging to the same compatibility group as pathoID).
# Each sample serves as negative control for every pathologyID that is neither an
# ancestor nor compatible with its own pathoID.
#
# If the metadata file has errors, log as many as possible and die.
sub parsePathologies {
    my $subName = (caller(0))[3];
    (@_ == 1) || die "E: $subName needs one arg";
    my ($pathosFile) = @_;

    (-f $pathosFile) || die "E in $subName: provided file $pathosFile doesn't exist\n";

    # hashrefs to be returned
    my ($ancestorsR, $compatR) = ({}, {});

    # we will populate temp hash %compatGroups as we go: key==compatGroup id,
    # value==arrayref of pathologyIDs that belong to this compatGroup
    my %compatGroups = ();
    # also as we go: %noFiles, key==pathoID, value==1 (not defined if no_files==0)
    my %noFiles = ();
    # also as we go: %parents, key==pathoID, value==arrayref of parent pathoIDs
    my %parents = ();
    
    my $workbook = Spreadsheet::XLSX->new("$pathosFile");
    (defined $workbook) ||
        die "E in $subName: no workbook\n";
    ($workbook->worksheet_count() == 1) ||
        warn "W in $subName: expecting a single worksheet, got ".$workbook->worksheet_count().
        ", ignoring all worksheets except the first\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($pathoCol,$noFilesCol,$isaCol,$compatCol) = (-1,-1,-1,-1);
    foreach my $col ($colMin..$colMax) {
        my $cell = $worksheet->get_cell($rowMin, $col);
        # if column has no header just ignore it
        (defined $cell) || next;
        my $val = $cell->unformatted();
        if ($val eq "pathologyID") { $pathoCol = $col; }
        elsif ($val eq "no_files") { $noFilesCol = $col; }
        elsif ($val eq "is_a") { $isaCol = $col; }
        elsif ($val eq "compatibility groups") { $compatCol = $col; }
    }
    ($pathoCol >= 0) ||
        die "E in $subName: missing required column title: 'pathologyID'\n";
    ($noFilesCol >= 0) ||
        die "E in $subName: missing required column title: 'no_files'\n";
    ($isaCol >= 0) ||
        die "E in $subName: missing required column title: 'is_a'\n";
    ($compatCol >= 0) ||
        die "E in $subName: missing required column title: 'compatibility groups'\n";

    # when parsing data rows, log any errors and only die at the end
    my $errorsFound = 0;
    
    foreach my $row ($rowMin+1..$rowMax) {
        my $patho = $worksheet->get_cell($row, $pathoCol);
        # skip lines without a pathoID
        ($patho) || next;
        $patho = $patho->unformatted();
        # require alphanumeric+underscore strings
        if ($patho !~ /^\w+$/) {
            warn "E in $subName: pathologyIDs must be alphanumeric+underscore strings, found \"$patho\" in row ",$row+1,"\n";
            $errorsFound++;
            next;
        }
        if (defined $parents{$patho}) {
            warn "E in $subName: there are 2 lines with the same pathologyID $patho\n";
            $errorsFound++;
            next;
        }
        $parents{$patho} = [];

        # is_a must be a comma-separated list of pathologyIDs
        my $isas = $worksheet->get_cell($row, $isaCol);
        if (defined $isas) {
            $isas = $isas->unformatted();
            foreach my $isa (split(/,/, $isas)) {
                # just save the strings for now, when all parsed we will check that
                # they are valid pathologyIDs
                if ( grep( /^$isa$/, @{$parents{$patho}})) {
                    warn "E in $subName: is_a contains the same pathologyID $isa twice for $patho\n";
                    $errorsFound++;
                    next;
                }
                push(@{$parents{$patho}}, $isa);
            }
        }
        
        # no_files can be empty or 0 (false), or 1
        my $nf = $worksheet->get_cell($row, $noFilesCol);
        if (defined $nf) {
            $nf = $nf->unformatted();
            if ($nf eq '1') {
                $noFiles{$patho} = 1;
            }
            elsif ($nf ne '0') {
                warn "E in $subName: no_files must be empty or 0 or 1, found $nf for $patho\n";
                $errorsFound++;
            }
        }
        
        # 'Compatibility group' must be a comma-separated list of group identifiers (alphanum strings)
        my $compats = $worksheet->get_cell($row, $compatCol);
        (defined $compats) || next;
        # ignore if no_files==1
        ($noFiles{$patho}) && next;
        $compats = $compats->unformatted();
        foreach my $cg (split(/,/, $compats)) {
            if ($cg !~ /^\w+$/) {
                warn "E in $subName: compat group identifiers must be alphanumeric strings, found $cg for $patho\n";
                $errorsFound++;
                next;
            }
            (defined $compatGroups{$cg}) || ($compatGroups{$cg} = []);
            push(@{$compatGroups{$cg}}, $patho);
        }
    }

    ################################
    # DONE PARSING FILE

    # now check the %parents strings
    foreach my $patho (keys %parents) {
        foreach my $parent (@{$parents{$patho}}) {
            if (! defined $parents{$parent}) {
                warn "E in $subName: is_a must be a comma-separated list of pathologyIDs, found $parent for $patho\n";
                $errorsFound++;
            }
        }
    }
    
    if ($errorsFound) {
        die "E in $subName: encountered $errorsFound errors while parsing $pathosFile, please fix the file.\n";
    }

    # populate $ancestorsR: naive algorithm but should be fine, we don't expect a huge onthology
    foreach my $patho (keys %parents) {
        (defined $noFiles{$patho}) && next;
        $ancestorsR->{$patho} = {};
        my @toExamine = @{$parents{$patho}};
        while (my $anc = shift(@toExamine)) {
            (defined $noFiles{$anc}) || ($ancestorsR->{$patho}->{$anc} = 1);
            push(@toExamine, @{$parents{$anc}});
        }
    }

    # initialize $compatR
    foreach my $patho (keys %$ancestorsR) {
        $compatR->{$patho} = {};
    }
    # populate $compatR: start with all descendants of each pathoID
    foreach my $patho (keys %$ancestorsR) {
        foreach my $anc (keys %{$ancestorsR->{$patho}}) {
            $compatR->{$anc}->{$patho} = 1;
        }
    }
    # add compat info from %compatGroups
    foreach my $cgs (values(%compatGroups)) {
        foreach my $patho (@$cgs) {
            foreach my $compatPatho (@$cgs) {
                ($compatPatho eq $patho) && next;
                # a patho is never "compatible" with one of its ancestors (it's part of it)
                (defined $ancestorsR->{$patho}->{$compatPatho}) && next;
                $compatR->{$patho}->{$compatPatho} = 1;
            }
        }
    }
    
    return($ancestorsR, $compatR);
}


#################################################################

# Parse samples metadata XLSX file: required columns are "sampleID",
# "specimenID", "patientID", "pathologyID", and "Causal gene" (they
# can be in any order but they MUST exist). An optional column "Sex"
# is parsed if it exists.
# The optional second argument, if provided and non-empty, is the
# pathologies metadata XLSX file; it is used to make sure that every
# pathologyID in samples.xlsx is defined in pathologies.xlsx (ie
# sanity-check for typoes), without no_files=1 (if we have a sample for
# a pathologyID, we MUST make output files for this pathologyID).
#
# Return a list of 4 (or 5 if "Sex" column exists) hashrefs, the caller can use
# whichever it needs:
# - sample2patho: key is sampleID, value is pathologyID
# - sample2specimen: key is sampleID, value is specimenID
# - sample2patient: key is sampleID, value is patientID if it exists, specimenID otherwise
# - sample2causal: key is sampleID, value is causal gene if it exists (undef otherwise)
# - sample2sex (only if the sex column exists): key is sampleID, value is "M" or "F"
#
# If the metadata file has errors, log as many as possible and die.
sub parseSamples {
    my $subName = (caller(0))[3];
    (@_ == 1) || (@_ == 2) || die "E: $subName needs one or two args";
    my ($samplesFile, $pathosFile) = @_;

    (-f $samplesFile) ||
        die "E in $subName: provided samples file $samplesFile doesn't exist\n";
    (! $pathosFile) || (-f $pathosFile) ||
        die "E in $subName: optional pathologies file $pathosFile was provided but doesn't exist\n";

    # sample2patho: key is sampleID, value is pathologyID
    my %sample2patho;
    # sample2specimen: key is sampleID, value is specimenID
    my %sample2specimen;
    # sample2patient: key is sampleID, value is patientID if it exists, specimenID otherwise
    my %sample2patient;
    # sample2causal: key is sampleID, value is causal gene if it exists (undef otherwise)
    my %sample2causal;
    # sample2sex: key is sampleID, value is the sample's sex
    my %sample2sex;

    ################
    # parse $pathosFile immediately if it was provided
    my $ancestorsR;
    if ($pathosFile) {
        # grab just the first returned scalar
        $ancestorsR = &parsePathologies($pathosFile);
    }

    ################
    # parse headers
    my $workbook = Spreadsheet::XLSX->new("$samplesFile");
    (defined $workbook) ||
        die "E in $subName: no workbook\n";
    ($workbook->worksheet_count() == 1) ||
        warn "W in $subName: expecting a single worksheet, got ".$workbook->worksheet_count().
        ", ignoring all worksheets except the first\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($sampleCol, $pathoCol, $specimenCol, $patientCol, $causalCol, $sexCol) = (-1,-1,-1,-1,-1,-1);
    foreach my $col ($colMin..$colMax) {
        my $cell = $worksheet->get_cell($rowMin, $col);
        # if column has no header just ignore it
        (defined $cell) || next;
        my $val = $cell->unformatted();
        if ($val eq "sampleID") { $sampleCol = $col; }
        elsif ($val eq "pathologyID") { $pathoCol = $col; }
        elsif ($val eq "specimenID") { $specimenCol = $col; }
        elsif ($val eq "patientID") { $patientCol = $col; }
        elsif ($val eq "Causal gene") { $causalCol = $col; }
        elsif ($val eq "Sex") { $sexCol = $col; }
    }
    ($sampleCol >= 0) ||
        die "E in $subName: missing required column title: 'sampleID'\n";
    ($pathoCol >= 0) ||
        die "E in $subName: missing required column title: 'pathologyID'\n";
    ($specimenCol >= 0) ||
        die "E in $subName: missing required column title: 'specimenID'\n";
    ($patientCol >= 0) ||
        die "E in $subName: missing required column title: 'patientID'\n";
    ($causalCol >= 0) ||
        die "E in $subName: missing required column title: 'Causal gene'\n";
    # Sex is optional, don't test

    ################
    # parse data rows
    
    # when parsing data rows, log any errors and only die at the end
    my $errorsFound = 0;

    foreach my $row ($rowMin+1..$rowMax) {
        my $sample = $worksheet->get_cell($row, $sampleCol);
        if (! $sample) {
            warn "E in $subName, row ",$row+1,": every row MUST have a sampleID (use 0 for obsolete samples)\n";
            $errorsFound++;
            next;
        }
        $sample = $sample->unformatted();
        # skip "0" lines == obsolete samples
        ($sample eq "0") && next;
        if (defined $sample2patho{$sample}) {
            warn "E in $subName: found 2 lines with same sampleID $sample\n";
            $errorsFound++;
            next;
        }
        
        ################ sample2patho
        my $patho = $worksheet->get_cell($row, $pathoCol);
        if (! $patho) {
            warn "E in $subName, row ",$row+1,": every row MUST have a pathologyID\n";
            $errorsFound++;
            next;
        }
        $patho = $patho->unformatted();
        if ($pathosFile) {
            if (! defined $ancestorsR->{$patho}) {
                warn "E in $subName, row ",$row+1,": pathologyID $patho is not defined (or has no_files=1) in $pathosFile\n";
                $errorsFound++;
                next;
            }
        }
        else {
            # require alphanumeric strings
            if ($patho !~ /^\w+$/) {
                warn "E in $subName, row ",$row+1,": pathologyIDs must be alphanumeric strings, found \"$patho\"\n";
                $errorsFound++;
                next;
            }
        }
        $sample2patho{$sample} = $patho;

        ################ sample2specimen
        my $specimen = $worksheet->get_cell($row, $specimenCol);
        if (! $specimen) {
            warn "E in $subName, row ",$row+1,": every row MUST have a specimenID\n";
            $errorsFound++;
            next;
        }
        $specimen = $specimen->unformatted();
        # require alphanumeric or dashes (really don't want spaces or shell metachars)
        if ($specimen !~ /^[\w-]+$/) {
            warn "E in $subName, row ",$row+1,": specimenIDs must be alphanumeric (dashes allowed), found \"$specimen\"\n";
            $errorsFound++;
            next;
        }
        $sample2specimen{$sample} = $specimen;

        ################ sample2patient
        my $patient = $worksheet->get_cell($row, $patientCol);
        # patientID value is optional, we use specimenID if it's missing
        if ($patient) {
            $patient = $patient->unformatted();
            # require alphanum (or underscores as always) or dashes, this could probably
            # be relaxed
            if ($patient !~ /^[\w-]+$/) {
                warn "E in $subName, row ",$row+1,": patientIDs must be alphanumeric (dashes allowed), found \"$patient\"\n";
                $errorsFound++;
                next;
            }
        }
        else {
            $patient = $specimen;
        }
        $sample2patient{$sample} = $patient;

        ################ sample2causal
        my $causal = $worksheet->get_cell($row, $causalCol);
        if ($causal) {
            $causal = $causal->unformatted();
            # these are HUGO gene names -> must be alphanum+dashes
            if ($causal !~ /^[\w-]+$/) {
                warn "E in $subName, row ",$row+1,": causalGene must be alphanumeric (dashes allowed), found \"$causal\"\n";
                $errorsFound++;
                next;
            }
            $sample2causal{$sample} = $causal;
        }

        ################ sample2sex
        if ($sexCol >= 0) {
            my $sex = $worksheet->get_cell($row, $sexCol);
            if (! $sex) {
                warn "E in $subName, row ",$row+1,": since we have the optional \"Sex\" column, every row MUST have a sex\n";
                $errorsFound++;
                next;
            }
            $sex = $sex->unformatted();
            if (($sex ne "M") && ($sex ne "F")) {
                warn "E in $subName, row ",$row+1,": sex must be M or F, found \"$sex\"\n";
                $errorsFound++;
                next;
            }
            $sample2sex{$sample} = $sex;
        }
    }
    
    ################
    if ($errorsFound) {
        die "E in $subName: encountered $errorsFound errors while parsing $samplesFile, please fix the file.\n";
    }
    elsif ($sexCol >= 0) {
        return(\%sample2patho, \%sample2specimen, \%sample2patient, \%sample2causal, \%sample2sex);
    }
    else {
        return(\%sample2patho, \%sample2specimen, \%sample2patient, \%sample2causal);
    }
}


#################################################################

# Parse candidateGenes metadata XLSX file(s): required columns are
# "Gene", "pathologyID", and "Confidence score" (can be in any order
# but they MUST exist).
# There can be several candidateGenes files, comma-separated.
# The second argument is the samples metadata file: causalGenes in
# this file are added as candidate genes with score=5.
# The third argument $verbose is a Boolean: causalGenes that are absent
# from the candidateGenes file are reported only if $verbose (useful
# to avoid redundant stderr output).
# The optional fourth argument, if povided and non-empty, is the
# pathologies metadata XLSX file; it is used to make sure every
# pathologyID is defined in pathologies.xlsx (ie sanity-checking).
#
# Return a hashref:
# - key is a pathologyID
# - value is a hashref, with keys==Genes and value==confidenceScore
#
# If the metadata files have errors, log as many as possible and die.
sub parseCandidateGenes {
    my $subName = (caller(0))[3];
    (@_ == 3) || (@_ == 4) || die "E: $subName needs three or four args";
    my ($candidatesFiles, $samplesFile, $verbose, $pathosFile) = @_;

    my @candidateFiles = split(/,/, $candidatesFiles);

    foreach my $cf (@candidateFiles) {
        (-f $cf) || die "E in $subName: provided candidateGenes file $cf doesn't exist\n";
    }
    (-f $samplesFile) ||
        die "E in $subName: provided samples file $samplesFile doesn't exist\n";
    (! $pathosFile) || (-f $pathosFile) ||
        die "E in $subName: optional pathologies file $pathosFile was provided but doesn't exist\n";


    # %knownCandidateGenes: key==$pathoID, value is a hashref whose keys 
    # are gene names and values are the "Confidence score" from a $candidatesFile,
    # or 5 if the gene is "Causal" for a $pathoID patient in $samplesFile.
    my %knownCandidateGenes = ();


    ################
    # parse $samplesFile and $pathosFile immediately
    
    # $sample2pathoR: hashref, key==sample id, value is this sample's pathologyID
    my $sample2pathoR;
    # $samples2causalR: hashref, key==sample id, value == causal gene (HGNC gene name)
    my $sample2causalR;
    # $pathologiesR: hashref, keys are valid pathologies (just ignore the hash values)
    my $pathologiesR;
    {
        my @parsed;
        if ($pathosFile) {
            @parsed = &parseSamples($samplesFile, $pathosFile);
            # grab just the first returned scalar (this is ancestorsR, fine)
            $pathologiesR = &parsePathologies($pathosFile);
        }
        else {
            @parsed = &parseSamples($samplesFile);
        }
        $sample2pathoR = $parsed[0];
        $sample2causalR = $parsed[3];
    }

    ################
    # parse candidateGenes files
    
    # we want to log all errors and only die at the end
    my $errorsFound = 0;

    foreach my $candidatesFile (@candidateFiles) {
        my $workbook = Spreadsheet::XLSX->new("$candidatesFile");
        (defined $workbook) ||
            die "E in $subName: no workbook\n";
        ($workbook->worksheet_count() == 1) ||
            warn "W in $subName: expecting a single worksheet, got ".$workbook->worksheet_count().
            ", ignoring all worksheets except the first\n";
        my $worksheet = $workbook->worksheet(0);
        my ($colMin, $colMax) = $worksheet->col_range();
        my ($rowMin, $rowMax) = $worksheet->row_range();
        # check the column titles and grab indexes of our columns of interest
        my ($pathoCol, $geneCol,$scoreCol,$transmiCol) = (-1,-1,-1);
        foreach my $col ($colMin..$colMax) {
            my $cell = $worksheet->get_cell($rowMin, $col);
            # if column has no header just ignore it
            (defined $cell) || next;
            my $val = $cell->unformatted();
            if ($val eq "pathologyID") { $pathoCol = $col; }
            elsif ($val eq "Gene") { $geneCol = $col; }
            elsif ($val eq "Confidence score") { $scoreCol = $col; }
            elsif ($val =~ /^Transmission/) { $transmiCol = $col; }
        }
        ($pathoCol >= 0) ||
            die "E in $subName, parsing $candidatesFile: missing required column title: 'pathologyID'\n";
        ($geneCol >= 0) ||
            die "E in $subName, parsing $candidatesFile: missing required column title: 'Gene'\n";
        ($scoreCol >= 0) ||
            die "E in $subName, parsing $candidatesFile: missing required column title: 'Confidence score'\n";
        ($transmiCol >= 0) ||
            die "E in $subName, parsing $candidatesFile: missing required column title: 'Transmission...'\n";
        
        ################
        # parse data rows
        foreach my $row ($rowMin+1..$rowMax) {
            my $patho = $worksheet->get_cell($row, $pathoCol);
            my $gene = $worksheet->get_cell($row, $geneCol);
            my $score = $worksheet->get_cell($row, $scoreCol);
            my $transmission = $worksheet->get_cell($row, $transmiCol);

            # blank lines are allowed, silently skip them
            if ((! $patho) && (! $gene) && (! $score) && (! $transmission)) {
                next;
            }
            
            ################ pathology
            if (! $patho) {
                warn "E in $subName, parsing $candidatesFile row ", $row+1,
                    ": every row MUST have a pathologyID\n";
                $errorsFound++;
                next;
            }
            $patho = $patho->unformatted();
            if ($pathosFile) {
                if (! defined $pathologiesR->{$patho}) {
                    warn "E in $subName, parsing $candidatesFile row ", $row+1,
                        ": pathologyID $patho is not defined (or no_files=1) in the provided $pathosFile\n";
                    $errorsFound++;
                    next;
                }
            }
            else {
                # require alphanumeric strings
                if ($patho !~ /^\w+$/) {
                    warn "E in $subName, parsing $candidatesFile row ", $row+1,
                        ": pathologyIDs must be alphanumeric strings, found \"$patho\"\n";
                    $errorsFound++;
                    next;
                }
            }
            
            ################ gene
            if (! $gene) {
                warn "E in $subName, parsing $candidatesFile row ",$row+1,": every row MUST have a Gene\n";
                $errorsFound++;
                next;
            }
            $gene = $gene->unformatted();
            # these are HUGO gene names -> must be alphanum+dashes but allow (and ignore) trailing spaces
            if ($gene !~ /^([\w-]+)\s*$/) {
                warn "E in $subName, parsing $candidatesFile row ", $row+1,
                    ": Gene must be alphanumeric (dashes allowed), found \"$gene\"\n";
                $errorsFound++;
                next;
            }
            $gene = $1;

            ################ score
            if (! $score) {
                warn "E in $subName, parsing $candidatesFile row ", $row+1,
                    ": every row MUST have a Confidence score\n";
                $errorsFound++;
                next;
            }
            $score = $score->unformatted();
            # require alphanumeric strings
            if ($score !~ /^\w+$/) {
                warn "E in $subName, parsing $candidatesFile row ", $row+1,
                    ": Confidence score must be alphanumeric, found \"$score\"\n";
                $errorsFound++;
                next;
            }

            ################ transmission, can be empty
            if ($transmission) {
                $transmission = $transmission->unformatted();
                # require alphanumerics or -, ignore trailing spaces
                if ($transmission !~ /^([\w\-]+)\s*$/) {
                    warn "E in $subName, parsing $candidatesFile row ", $row+1,
                        ": Transmission must be alphanumeric ('-' allowed), found \"$transmission\"\n";
                    $errorsFound++;
                    next;
                }
            }
            
            ################ looks AOK, save data
            (defined $knownCandidateGenes{$patho}) ||
                ($knownCandidateGenes{$patho} = {});
            if (defined $knownCandidateGenes{$patho}->{$gene}) {
                warn "E in $subName, parsing $candidatesFile row ", $row+1,
                    ": this $gene - $patho association was already seen elsewhere\n";
                $errorsFound++;
                next;
            }
            $knownCandidateGenes{$patho}->{$gene} = $score;
            ($transmission) && ($knownCandidateGenes{$patho}->{$gene} .= ":$transmission");
        }
    }
    
    ################
    # add causalGenes to knownCandidateGenes with score 5
    foreach my $sample (keys %$sample2causalR) {
        my $patho = $sample2pathoR->{$sample};
        my $causal = $sample2causalR->{$sample};
        (defined $knownCandidateGenes{$patho}) || ($knownCandidateGenes{$patho} = {});
        # $causal could/should be in candidatesFiles
        ($knownCandidateGenes{$patho}->{$causal}) || (! $verbose) ||
            warn "I in $subName: gene $causal is marked causal of $patho for $sample, you may want to add it to a candidateGenes file\n";
        # in any case set score to 5
        $knownCandidateGenes{$patho}->{$causal} = 5;
    }
    
    ################
    if ($errorsFound) {
        die "E in $subName: encountered $errorsFound errors while parsing $candidatesFiles, please fix the file(s).\n";
    }
    else {
        return(\%knownCandidateGenes);
    }
}


# module loaded ok
1;
