#!/usr/bin/env python3

#
# Classify multifasta sequences based on DNN
#

import os
import sys
from io import StringIO
import pandas as pd
import numpy as np

import operator
from collections import Counter
from sklearn.preprocessing import LabelEncoder

import keras
from keras.models import load_model
from keras import utils

def classify(X, theta):
	"""
	if maximum of predictied probability (softmax) >= theta: argmax, else: Unknown
	"""
	maxima = np.apply_along_axis(np.max,1,X)
	arg_maxima = np.apply_along_axis(np.argmax,1,X)
	classes = [label_encoder.inverse_transform(arg_maxima[i])
		if maxima[i] > theta else "Unknown" for i in range(len(maxima))]
	# score is equal to activation of the best node, 1 if classified as unknown
	scores = [maxima[i] if maxima[i] > theta else 1 for i in range(len(maxima))]
	return (classes, scores)

model = load_model("/home/gi54cop/ma/dat/deep-learning/model20-sub.h5")
features = [x.strip() for x in open("/home/gi54cop/ma/dat/deep-learning/features.txt")]
targets = [x.strip() for x in open("/home/gi54cop/ma/dat/deep-learning/targets.txt")]
label_encoder = LabelEncoder().fit(targets)

# multifasta file containing sequences to be classified
fasta_path = sys.argv[1]
theta = float(sys.argv[2])

# prep samples
kmer_txt = os.popen("fasta2matrix.py -normalize frequency 6  %s" % fasta_path).read()
sample_df = pd.read_csv(StringIO(kmer_txt), index_col=0, sep="\t")
X = np.array(sample_df[features])

# prediction
classes, scores = classify(model.predict(X), theta)
res_df = pd.DataFrame()
res_df["sequence"] = sample_df.index
res_df["predicted"] = classes
res_df["score"] = scores
for ri, r in res_df.iterrows():
	if r.predicted != "Unknown": # skip Unknown samples
		print("%s\t%s\t%s" % (r.sequence, r.predicted, r.score))

