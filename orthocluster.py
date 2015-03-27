#!/usr/bin/python

'''
This script can be used to calculate how many identified Dmel
orthologs in other genomes are clustered in the same group 
as their Dmel ortholog, and how their motif compositions differ.
It requires clustering results from cluster.py and the orthologs csv
file resulting from recblast.py.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 17 March 2015
'''

import csv,itertools,re
from jellyfish import levenshtein_distance as jld

#CUSTOMIZE INPUT FILES
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data"
clusterfile = "results/clustering-string_allz-average.csv"
orthologs = "sequences/dmel-orthologs.csv"
#END CUSTOMIZATION

def orthos(gene,k,ol):
	'''
	Look up per gene what their Dmel ortholog is. NB this may be multiple, as it could be multiple isoforms.
	'''
	isoforms = [ol[0][i] for i,j in enumerate(ol[k]) if j == gene]
	return isoforms #list of all isoforms
 
def clusters(genes,cdict): #genes is a list
	'''
	Look up the cluster they are assigned (both the Dmel and the Spp gene)
	'''
	clusters = [cdict[g.split('|')[2]] for g in genes]
	return clusters #list of clusters, same order as genes

def wordcomp(i,j):
	'''
	Calculate pairwise distances for all the strings collected. As strings
	are regular expressions, first expand the re. and then calculate all distances
	pairwise. Return the minimal distance.
	'''
	# translate strings[i] and strings[j] to all possible expressions
	sisplit = [part.split('|') for part in re.split(r'\{(.*?)\}',i)]
	sjsplit = [part.split('|') for part in re.split(r'\{(.*?)\}',j)]
	si,sj = [],[]
	for x in itertools.product(*sisplit):
		si.append(''.join(x))
	for x in itertools.product(*sjsplit):
		sj.append(''.join(x))
	# check jld of all strings[i] options against all strings[j] options
	distance = []
	for sin in si:
		for sjn in sj:
			distance.append(jld(sin,sjn))
	# return the minimum distance
	return min(distance)

def distance(gene,dmorths,sdict): #genes is a list
	'''
	Calculate the minimum distance between the gene and its dm orthologs
	where distances are calculated between strings that summarize protein
	architecture.
	'''
	dists = []
	for dmo in dmorths:
		g = sdict[gene[0].split('|')[2]]
		d = sdict[dmo.split('|')[2]]
		dist = wordcomp(g,d)
		dists.append(dist)
	return min(dists)

#Read ortholog file
of = csv.reader(open("%s/%s" %(dbfolder,orthologs)))
ol1 = [line for line in of] # list from the file
ol2 = [[x[i] for x in ol1] for i in range(len(ol1[0]))] # transpose the list

#Read cluster file and make dictionary
cf = csv.reader(open("%s/%s" %(dbfolder,clusterfile)))
cdict,sdict = {},{}
for line in cf:
	name = line[2]
	cdict[name] = line[3]
	sdict[name] = line[4]

#make a list of all the unique orthologs per species, and collect the isoforms
for k in range(1,len(ol2)):
	allorth = list(set(ol2[k][1:]))
	allorth.remove('')
	print 'species:', ol2[k][0]
	print 'total no. detected orthologs:', len(allorth)
	clcount = 0
	dists = []
	for gene in allorth:
		dmorths = orthos(gene,k,ol2)
		gene = [gene] #the clusters function requires a list input
		cl_gene = clusters(gene,cdict)
		cl_orth = clusters(dmorths,cdict)
		if cl_gene[0] in cl_orth: #If one of the isoforms has the same cluster hit as the spp ortholog, score as 1
			clcount += 1
		else:
			dist = distance(gene,dmorths,sdict)
			dists.append(dist)
	print 'total no. orthologs in same cluster:', clcount
	print 'percentage:', clcount/float(len(allorth)) * 100
	print 'average distance:', sum(dists)/float(len(dists))
