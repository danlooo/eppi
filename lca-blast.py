#!/usr/bin/env python3

#
# blast tsv lca
#

from io import StringIO
from collections import Counter
import re
import pandas as pd
import sqlite3
import sys

def annotate_lineage(df, accCol, dbpath):
	"""
	@param accCol: name of colukn in df containing accessions
	@return: annotated data frame
	"""
	cur = sqlite3.connect(dbpath).cursor()
	taxInfo = {}
	accs = set([x.split(".")[0] for x in set(df[accCol])])
	accs = [x.split("|")[-1] for x in accs] # convert ncbi accs like gi|323371659|gb|CP002521.1| too acc

	for a in accs:
		res = cur.execute("""
			SELECT *
			FROM acc_lineage
			WHERE acc = \"%s\"
			LIMIT 1;""" % a).fetchone()
		taxInfo[a] = res

	def query_dict(acc, item):
		try:
			if "|" in acc:
				# convert ncbi accs like gi|323371659|gb|CP002521.1| to acc
				acc = acc.split("|")[-2]
			return taxInfo[acc.split(".")[0]][item]
		except:
			return None

	df["tax_id"] = df[accCol].apply(lambda x: query_dict(x,1))
	df["species"] = df[accCol].apply(lambda x: query_dict(x,2))
	df["genus"] = df[accCol].apply(lambda x: query_dict(x,3))
	df["family"] = df[accCol].apply(lambda x: query_dict(x,4))
	df["order"] = df[accCol].apply(lambda x: query_dict(x,5))
	df["class"] = df[accCol].apply(lambda x: query_dict(x,6))
	df["phylum"] = df[accCol].apply(lambda x: query_dict(x,7))
	df["superkingdom"] = df[accCol].apply(lambda x: query_dict(x,8))
	return df

def get_lca(df, query_col, max_evalue, min_support):
	"""
	returns lowest common ancestor for each sequence
	@param df: data frame containing all hits annotated with their lineage
	@param query_col: name of the column in df containing query sequences
	@max_evalue: hits with greater evalues than this will be excluded
	@min_support: a taxon must have at least this fraction to be a lca
	"""
	lca_df = pd.DataFrame(columns = ["lca", "lca_level"])

	for seq in set(df[query_col]):
		hits = df.query("%s == '%s' and evalue <= %s" % (query_col, seq, max_evalue))

		for l in ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]:
			counts = Counter(hits[l])
			num_counts = sum(counts.values())
			counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)

			# skip empty assignments
			if len(counts) == 0: continue

			if counts[0][1] / num_counts > min_support:
				# most abundant taxa at level l has enough support and is therefore lca
				lca_df.loc[seq, "lca_level"] = l
				lca_df.loc[seq, "lca"] = counts[0][0] # first hit is lca
				lca_df.loc[seq, "support"] = counts[0][1] / num_counts
				break
	# filter failed results	
	lca_df = lca_df.query("lca != ''")
	return lca_df

blast_tsv_file = sys.argv[1]
tax_db_path = sys.argv[2]

df = pd.read_csv(blast_tsv_file, sep="\t", header=None)
df.columns = ["qseqid", "sseqid", "sstart", "send", "evalue", "btop"]
df = df.sort_values("qseqid")

df = annotate_lineage(df, "sseqid", tax_db_path)
lca_df = get_lca(df, "qseqid", 0.001, 0.51)
lca_df.to_csv(blast_tsv_file+"-lca.csv")
