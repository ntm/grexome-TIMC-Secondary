#!/bin/sh

# 18/09/2019
# NTM


# Build the list of canonical transcripts from Ensembl, by 
# connecting directly to the ensembl MYSQL database.
# Output to STDOUT a single column with ENST identifiers for 
# canonical transcripts.
#
# This is functionally equivalent to listCanonicalTranscripts_API.pl
# but the latter took 16h to run, this takes seconds...
# It can't be run at TIMC because it needs access to TCP/3306,
# which is filtered at TIMC. Just run it from eg nicofree.

DATABASE=homo_sapiens_core_97_38

echo "USE $DATABASE ; select t.stable_id from transcript t, gene g where t.transcript_id = g.canonical_transcript_id ;" | mysql -u anonymous -h ensembldb.ensembl.org | tail -n +2 | sort -u


