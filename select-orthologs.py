#!/usr/bin/python

'''
This script is used to select ortholog combinations of specific
interest, for further detailed study.

Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 30 March 2016
'''

import pandas as pd
import config
from collections import Counter

detcompcsv = pd.read_csv("%s/%s/%s-sppall_orthcomp-detail.csv" %(config.mainfolder,config.resfolder,config.idr),header=None)

#columns in the detcomp database
col1 = list(detcompcsv[detcompcsv.columns[0]]) #name of protein
col2 = list(detcompcsv[detcompcsv.columns[1]]) #sequence of protein (in motifs)
col3 = list(detcompcsv[detcompcsv.columns[2]]) #name of orthologous protein
col4 = list(detcompcsv[detcompcsv.columns[3]]) #sequence of orthologous protein

#generate a dictionary that contains the protein name as key, and the sequence as value
prots = col1+col3
seqs = col2+col4
mseqs = {a:b for a,b in zip(prots,seqs)}

#make a dictionary collecting all orthologs of a specific protein (key) as a list (value)
orthdict = {}
for a,b in zip(col1,col3):
	if a in orthdict:
		orthdict[a] += [b]
	else:
		orthdict[a] = [b]

#counts the iterations of each protein in the first column
countorth = Counter(col1)

#print to terminal those proteins (and their sequences) that have > 5 orthologs
for o in countorth:
	if countorth[o] > 5:
		print "\n",o
		for p in orthdict[o]:
			print p, "\t", mseqs[p]


#print config.translationdict_inv2
