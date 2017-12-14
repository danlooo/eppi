#!/bin/bash

#
# get statistics about analysis
#

echo reads_file,$reads_file >> out/stat.csv
echo reads_count,$[$(wc -l $reads_file | cut -f 1 -d " ") / 4] >> out/stat.csv
echo contigs_count,$(grep ^">" out/contigs.fasta | wc -l) >> out/stat.csv

echo hmm_lca_contigs_count,$(wc -l out/hmm/lca_df.csv | cut -f 1 -d " ") >> out/stat.csv
