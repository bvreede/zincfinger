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
from collections import Counter
import pylab as pl
import numpy as np

#CUSTOMIZE INPUT FILES
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data"
clusterfile = "results/clustering-string_all2-average.csv"
orthologs = "sequences/dmel-orthologs.csv" # ["sequences/%s-orthologs.csv" %spp for spp in ['dmel','dpul','smar','isca','turt','tcas']]
seqvsseq = "results/orthologs/dmeltoother.csv"
resultsummary = "results/orthologs/orthocomp_summary.csv"
hmimage = "images/ortholog_combinations.svg"
#END CUSTOMIZATION

#open resultfiles and make resultlists
seqcomp = open("%s/%s" %(dbfolder,seqvsseq), "w")
ressum = open("%s/%s" %(dbfolder,resultsummary), "w")
m2m,letters = [],[] #lists to collect combinations of letters (orthologs where two different domain classes are found in the same location) and the total appearance of those letters, respectively. 

#TAKE THIS FROM FINDMOTIF.PY: the alphabet and corresponding motiflist
motiflist = ['2_7_4','2_8_3','2_9_3','2_10_5','2_11_3','2_11_4','2_12_2','2_12_3','2_12_4','2_12_5','2_12_6','2_13_3',
'2_13_4','2_14_3','2_14_4','2_15_4','3_8_3','4_12_3','4_12_4','4_15_3']
alphabet = """ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`ABCDEFGHIJKLMNPQRSTUVWXYabcdefghijklmnopqrstuvwxy1234567890!@%^&*()_-+={}[]:;"'|\<,>.?/~`"""
translationdict = {motif: alphabet[a] for a,motif in enumerate(motiflist)} # dictionary of motifs and the corresponding string element


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
	pairwise. Return the minimal distance (with and without spaces).
	Additionally: find out which domain classes 'replace' each other in orthologs.
	Save these results in global lists.
	'''
	global m2m
	global letters
	# translate strings[i] and strings[j] to all possible expressions
	sisplit = [part.split('|') for part in re.split(r'\{(.*?)\}',i)]
	sjsplit = [part.split('|') for part in re.split(r'\{(.*?)\}',j)]
	si,sj = [],[]
	for x in itertools.product(*sisplit):
		si.append(''.join(x))
	for x in itertools.product(*sjsplit):
		sj.append(''.join(x))
	# check jld of all strings[i] options against all strings[j] options
	distance,distance2,sequencecomp = [],[],[]
	for sin in si:
		for sjn in sj:
			distance.append(jld(sin,sjn))
			#remove Zs and do it again
			sin = sin.replace('Z','')
			sjn = sjn.replace('Z','')
			distance2.append(jld(sin,sjn))
			sequencecomp.append("%s,%s" %(sin,sjn))
	# if there are differences: which motif turned into which?
	if min(distance2) > 0: #use the sequences without space/Z
		for k,d in enumerate(distance2):
			if d == min(distance2):
				#distance2 and sequencecomp have the same order, so sequencecomp[k] has sequences corresponding to distance = d
				comps = sequencecomp[k].split(',')
				sik,sjk = comps[0],comps[1]
				if len(sik) == len(sjk):
					for ki,si in enumerate(sik):
						letters.append(sjk[ki])
						letters.append(si)
						if si != sjk[ki]:
							fs = frozenset([si,sjk[ki]])
							m2m.append(fs)
					seqcomp.write("%s,%s\n" %(sik,sjk))
	# return the minimum distance
	return min(distance),min(distance2)

def distance(gene,dmorths,sdict): #NB: genes is a list
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

#Header line of results file
ressum.write("Species,Orthologs,Orth in same cluster,percentage,Identical orth,percentage,Levenshtein total,Levenshtein average,Identical orth (no space),Levenshtein total (no space),Levenshtein average (no space)\n")

for k in range(1,len(ol2)):
	#make a list of all the unique orthologs per species
	allorth = list(set(ol2[k][1:]))
	allorth.remove('')
	print 'species:', ol2[k][0]
	print 'total no. detected orthologs:', len(allorth)
	clcount,idcount,idcountnoZ = 0,0,0
	dists,distsnoZ = [],[]
	for gene in allorth:
		dmorths = orthos(gene,k,ol2) #collect all Dmel isoforms that are mentioned as orthologs for this gene
		gene = [gene] #the clusters function requires a list input
		cl_gene = clusters(gene,cdict)
		cl_orth = clusters(dmorths,cdict)
		if cl_gene[0] in cl_orth: #If one of the isoforms has the same cluster hit as the spp ortholog, score as 1
			clcount += 1
		#calculate the minimum distance between this gene and the dmel ortholog(s)
		dist,distnoZ = distance(gene,dmorths,sdict)
		dists.append(dist)
		distsnoZ.append(distnoZ)
		if dist == 0: # if the domain structures are identical, score as 1
			idcount += 1
		if distnoZ == 0:
			idcountnoZ += 1
	print 'total no. orthologs in same cluster:', clcount
	print 'percentage of orthologs in same cluster:', clcount/float(len(allorth)) * 100
	print 'number of identical orthologs:', idcount
	print 'percentage of identical orthologs:', idcount/float(len(allorth)) * 100
	print 'total levenshtein distance:', sum(dists)
	print 'average distance:', sum(dists)/float(len(dists))
	print 'number of identical orthologs (without space):', idcountnoZ
	print 'total levenshtein distance (without space):', sum(distsnoZ)
	print 'average distance (without space):', sum(distsnoZ)/float(len(distsnoZ))
	print '\n'
	ressum.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(ol2[k][0],len(allorth),clcount,clcount/float(len(allorth)),idcount,idcount/float(len(allorth)),sum(dists),sum(dists)/float(len(dists)),idcountnoZ,sum(distsnoZ),sum(distsnoZ)/float(len(distsnoZ))))

seqcomp.close()
ressum.close()

'''
Part II:
Which domain classes are parallel in non-identical orthologs?
Make a heatmap of this data.
'''
combocount = Counter(m2m) #dictionary with frozenset-motifcombinations, and their frequency
letterscount = Counter(letters) #dictionary with individual letters, and their frequency

#collect data for the heatmap in a 2d list
table = []
for m1 in motiflist:
	l1 = translationdict[m1] #corresponding string element of main motif
	total = letterscount[l1] #total frequency of this motif in the dataset
	row = [] #empty row that will collect relative frequency data
	for m2 in motiflist:
		l2 = translationdict[m2] #corresponding string element of comparing motif
		fs = frozenset([l1,l2]) #combined frozenset of main and comparing motif
		if fs in combocount: #if this combination is found:
			freq = float(combocount[fs]) #give total frequency of combination (float to enable float result)
		else:
			freq = 0
		row.append(freq/total)
	table.append(row)


###HEATMAP###
data = pl.array(table)
colourformap = "YlOrBr"
fig,ax = pl.subplots()
heatmap = pl.pcolor(data, cmap=colourformap)
cbar = pl.colorbar(heatmap)
	
# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)
ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)
pl.axis('tight') #remove the white bar
ax.invert_yaxis() #make sure it starts counting from the top
	
#make the labels
ax.set_xticklabels(motiflist, minor=False, rotation=90)
ax.set_yticklabels(motiflist, minor=False)
	
# save the figure
pl.savefig("%s/%s" %(dbfolder,hmimage), dpi = 300)
