#!/usr/bin/env python3

#
# get all contig information and put them into one table
#

import pandas as pd
from glob import glob

# vrap basics
contigs_df = pd.read_csv("out/vrap/vrap_summary.csv", sep=";", usecols=[0,1,2,5], index_col=0)
contigs_df.columns = ["length", "orf_dens", "virus"]

# Read coverage
coverage_df = pd.read_csv("out/read-cov/coverage.tsv", sep="\t", header=None)
contigs_df["read_cov"] = coverage_df.groupby([0])[2].sum() / contigs_df.length

#HMM lca
hmm_lca_df = pd.read_csv("out/hmm/lca_df.csv", index_col=0)
contigs_df = contigs_df.join(hmm_lca_df, how="left")

# blast lca
blastx_lca_df = pd.read_csv("out/diamond/blastx-contigs-hits.tsv-lca.csv", index_col=0)
blastx_lca_df.columns = ["blastx_lca", "blastx_lca_level", "blastx_support"]
contigs_df = contigs_df.join(blastx_lca_df, how="left")

blastn_lca_df = pd.read_csv("out/blast/blastn-contigs-nt.tsv-lca.csv", index_col=0)
blastn_lca_df.columns = ["blastn_lca", "blastn_lca_level", "blastn_support"]
contigs_df = contigs_df.join(blastn_lca_df, how="left")

# deep learning
deep_learning_df = pd.read_csv("out/deep-learning/dnn-hits.tsv", sep="\t", index_col=0)
deep_learning_df.columns = ["deep_learning_species", "deep_learning_score"]
contigs_df = contigs_df.join(deep_learning_df, how="left")

contigs_df.to_csv("out/contigs.csv")
