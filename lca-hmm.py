#!/usr/bin/env python3

#
# HMM taxon hit table to lca
#

import pandas as pd
import sqlite3
import sys
import glob
import matplotlib.pyplot as plt
from collections import Counter

tax_db = sys.argv[1]
taxhit_glob = sys.argv[2]
lca_df_file = sys.argv[3]

# sequence must have at least min_support of all annotated hits assigned to this taxon
# 0.5= half of hits of the sequence  must assign the same taxon to be lca
min_support = 0.51

taxon_df = pd.DataFrame()
for tax_hit_file in glob.glob(taxhit_glob):
	cur_df = pd.read_csv(tax_hit_file, index_col=0).fillna(0)
	# ignore frame from sequence name
	cur_df.index = cur_df.index.map(lambda x: "_".join(x.split("_")[:-1]))
	taxon_df = taxon_df.add(cur_df, fill_value=0).fillna(0)

cur = sqlite3.connect(tax_db).cursor()

lca_df = pd.DataFrame(columns=["hmm_lca", "hmm_lca_level", "hmm_support"])

lineage_df = pd.DataFrame()
lineage_df["tax_id"] = taxon_df.columns
for level in ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]:
	lineage_df[level] = lineage_df.tax_id.apply(lambda x:
		cur.execute("select [%s] from lineage where tax_id == '%s' limit 1;" % (level, x)).fetchone()[0])

for level in ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]:
	taxa = set(lineage_df[level]) - {""} # ignore unnamed ones
	level_df = pd.DataFrame(columns=taxa, index=taxon_df.index).fillna(0)

	# get all species tax ids for the taxon level 
	for taxon in taxa:
		tax_ids = [x.tax_id for ix, x in lineage_df.iterrows() if x[level] == taxon]
		level_df[taxon] = taxon_df[tax_ids].T.sum() # sum of all hits from species in the taxon
	best_taxons = level_df.T.apply(lambda x: x.argmax())
	best_taxons_counts = bestTaxons = level_df.T.apply(lambda x: x.max())	
	best_taxons_fractions = best_taxons_counts / level_df.T.sum()
	is_lca = best_taxons_fractions > min_support

	# for each sequence with lca hit
	for seq in is_lca[is_lca == True].index:
		# do not override e.g. species lca with family lca
		if seq not in lca_df.index:
			lca_df.loc[seq] = [best_taxons[seq], level, best_taxons_fractions[seq]]

lca_df.to_csv(lca_df_file)

try:
	# plot lca level pie chart
	c = Counter(lca_df.lca_level)
	plt.pie(list(c.values()), labels=c.keys())
	plt.axis("equal")
	plt.title("HMM lowest common ancestor")
	plt.savefig(lca_df_file+".pdf")
except:
	print("lca plot could not be plotted. No X Server?")

