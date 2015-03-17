#!/usr/bin/python

'''
This script can be used to calculate how many identified Dmel
orthologs in other genomes are clustered in the same group 
as their Dmel ortholog?
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 17 March 2015
'''

import csv

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

#Read ortholog file
of = csv.reader(open("%s/%s" %(dbfolder,orthologs)))
ol1 = [line for line in of] # list from the file
ol2 = [[x[i] for x in ol1] for i in range(len(ol1[0]))] # transpose the list

#Read cluster file and make dictionary
cf = csv.reader(open("%s/%s" %(dbfolder,clusterfile)))
cdict = {}
for line in cf:
	name = line[2]
	cdict[name] = line[3]

#make a list of all the unique orthologs per species, and collect the isoforms
for k in range(1,len(ol2)):
	allorth = list(set(ol2[k][1:]))
	allorth.remove('')
	print 'species:', ol2[k][0]
	print 'total no. detected orthologs:', len(allorth)
	clcount = 0
	for gene in allorth:
		dmorths = orthos(gene,k,ol2)
		gene = [gene] #the clusters function requires a list input
		cl_gene = clusters(gene,cdict)
		cl_orth = clusters(dmorths,cdict)
		if cl_gene[0] in cl_orth: #If one of the isoforms has the same cluster hit as the spp ortholog, score as 1
			clcount += 1
	print 'total no. orthologs in same cluster:', clcount
	print 'percentage: ', clcount/float(len(allorth)) * 100
