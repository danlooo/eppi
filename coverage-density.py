#!/usr/bin/env python3

#
# Calculates contig coverage density plots
# USAGE: coverage-density.py reads-against-contigs.bam out/path/
#

import os
import sys
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

sam_file = sys.argv[1]
out_path = sys.argv[2]

print("calculate raw coverage")
os.system("samtools depth %s > %scoverage.tsv" % (sam_file, out_path))

coverage_df = pd.read_csv(out_path+"coverage.tsv", header=None, sep="\t")
coverage_df.columns = ["sequence", "position", "coverage"]
seqs = sorted(set(coverage_df.sequence), key=lambda x: int(x.split("_")[1]))

# get only top 100 nodes
seqs = [x for x in seqs if "NODE" in x][:100]

# calculate position coverage heatmap
def scale(x, bins=100):
        res = []
        width = int(len(x) / bins)
        for i in range(bins):
                res += [np.mean(x[i*width:(i*width)+width])]
        return res

cov_heatmap_df = pd.DataFrame(columns=range(100))
print("calculate scaled position coverage")
for s in seqs:
        coverage = list(coverage_df.query("sequence == '%s'" % s).coverage)
        cov_heatmap_df.loc[s] = scale(coverage)
cov_heatmap_df.to_csv(out_path+"position-coverage-heatmap.csv")

# calculate coverage distribution heatmap
maxCov = 200
coverage_density_df = pd.DataFrame(columns = range(maxCov))
print("calculate coverage distribution")
for s in seqs:
        cs = coverage_df.query("sequence == '%s'" % s).coverage
        d = Counter(cs)
        coverage_density_df.loc[s] = [d[x] for x in range(maxCov)]
coverage_density_df.to_csv(out_path+"density-coverage-heatmap.csv")

if "DISPLAY" in os.environ:
	plt.figure(figsize=(10,5))
	sns.heatmap(coverage_density_df.fillna(0), cbar_kws={"label": "Anzahl Positionen (nt)"})
	plt.xlabel("Coverage")
	steps = 50
	plt.xticks(np.arange(0,200+steps,steps), np.arange(0,200+steps,steps))
	plt.ylabel("Contig")
	plt.yticks(np.arange(1,11)*10, reversed(range(0,100,10)))
	plt.savefig(out_path+"coverage_density.pdf")
	plt.close()
	# analyze position coverage heatmap
	plt.figure(figsize=(10,5))
	sns.heatmap(cov_heatmap_df, cbar_kws={"label": "Coverage"})
	plt.xlabel("Position in Contig (%)")
	plt.xticks(np.arange(10)*10, ["%s0%%" % x for x in range(10)])
	plt.yticks(np.arange(1,11)*10, reversed(range(0,100,10)))
	plt.ylabel("Contig")
	plt.savefig(out_path+"coverage_position.pdf")
	plt.close()
