#!/usr/bin/env python3

#
# eggnog hmm domain based search out to tax hit count table
#

import pandas as pd
import re
from io import StringIO
import sqlite3
import numpy as np
import sys

hmm_out_file = sys.argv[1]
tax_hit_out_file = sys.argv[2]

hmm_out_txt = open(hmm_out_file, "r").read()
hmm_out_txt = re.sub(" +", "\t", hmm_out_txt) # at least one space to tab
hmm_out_txt = re.sub("#.*\n", "", hmm_out_txt) # remove comment lines

df = pd.read_csv(StringIO(hmm_out_txt), sep="\t", header=None, usecols=[0,3,6])
df.columns = ["query", "target", "evalue"]

con = sqlite3.connect(sys.argv[3])
cur = con.cursor()

taxDict = {}
for t in set(df["target"]):
	taxDict[t] = cur.execute("select * from acc_lineage where acc == '%s'"
		% t).fetchall()

# filter significant hits
df = df.query("evalue <= 0.001")

#get taxcount matrix
taxs = []
for k,v in taxDict.items():
	# for each cluster
	for c in v:
		# add tax_id
		taxs += [c[1]]
taxs = list(set(taxs)) # get uniques
querys = list(set(df["query"]))
taxCountDf = pd.DataFrame(columns=taxs, index=querys)
taxCountDf = taxCountDf.fillna(0)

# for each hmm hit
for ri, r in df.iterrows():
	# get tax id hits for cluster
	taxIds = [x[1] for x in taxDict[r.target]]
	for tax in taxIds:
		# count hit
		taxCountDf.loc[r.query, tax] += 1
taxCountDf.to_csv(tax_hit_out_file)
