#!/usr/bin/python
'''
This script can be used to create clusters of a fasta file with
strings (without regard for the exact meaning of the individual elements!).
It works through calculating Levenshtein distances of string pairs.
It returns a csv file with several cluster options based on the cluster
threshold.
The script requires several additional modules (e.g. numpy, scipy, ete2) that need
to be installed.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 15 October 2014
'''

import scipy, pylab, re, itertools, csv, os, math, config, sys
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
from jellyfish import levenshtein_distance as jld
from numpy import triu_indices as nti
from numpy import apply_along_axis as naaa
from ete2 import Tree
from random import shuffle


#input from command line: the motif sequence file (translated string) that needs to be clustered.
if len(sys.argv) <= 1:
	sys.exit("USAGE: cluster.py path/to/inputfile (input needs to be a motif sequence fasta file).\nOutputfolders are indicated in the script; edit the script if you want to alter them.")

infile = sys.argv[1]
infilebrev = infile.split('/')[-1].split('_')[0]
#NB! filenames should always start with an input ID specifier (separate elements in dashes) and end with output ID specifiers.
#e.g.: 150525-dmel_seq.fa or 141212-tcas_heatmap.svg


##### INPUT SPECIFICATIONS: CUSTOMIZE HERE! #####

#options for clustering:
clustermeth = "average"
threshold = 1.1547
clustercrit = "inconsistent"

#output files and folders:
#dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results"
clusterorder = "%s/%s/%s_clusters" %(config.mainfolder,config.resfolder,infilebrev) #database that contains group data
clusterevolv = "%s/%s/%s_cluster-newick" %(config.mainfolder,config.evfolder,infilebrev) #evolview input file in newick format
clustercolour = "%s/%s/%s_cluster-colours" %(config.mainfolder,config.evfolder,infilebrev) #evolview input file labeling clusters


#### END OF CUSTOMIZATION! PLEASE DON'T EDIT BELOW #####
# set up the measurement of levenshtein distance of a combination of strings
def wordcomp(coord):
	'''
	Calculate pairwise distances for all the strings collected. As strings
	are regular expressions, first expand the re. and then calculate all distances
	pairwise. Return the minimal distance.
	'''
	i,j = coord
	# translate strings[i] and strings[j] to all possible expressions
	sisplit = [part.split('|') for part in re.split(r'\{(.*?)\}',strings[i])]
	sjsplit = [part.split('|') for part in re.split(r'\{(.*?)\}',strings[j])]
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

def sof_tree2newick(T):
	'''
	Solution taken directly from stackoverflow: http://stackoverflow.com/questions/9364609/converting-ndarray-generated-by-	hcluster-into-a-newick-string-for-use-with-ete2
	Takes a tree object and returns the tree in newick format.
	Requires ete2.
	'''
	#ete2 section
	root = Tree()
	root.dist = 0
	root.name = "root"
	item2node = {T: root}
	
	to_visit = [T]
	while to_visit:
	    	node = to_visit.pop()
		cl_dist = node.dist /2.0
		for ch_node in [node.left, node.right]:
			if ch_node:
				ch = Tree()
				ch.dist = cl_dist
				ch.name = str(ch_node.id)
				item2node[node].add_child(ch)
				item2node[ch_node] = ch
				to_visit.append(ch_node)
	# This is your ETE tree structure
	return root.write() #translate the tree to a newick string


# Read the input fasta file; extract:
# string [with motif data]
# ID labels [customized to separate on | character] + genenames and protein numbers separately
inputfile = open(infile)
strings,gID,genenames,protnumbers = [],[],[],[] #will collect the clusterable data (1) + header information (2-4)
### SECTION BELOW FLAGGED FOR FIXING: USE FASTA READER TO DICTIONARY ###
for line in inputfile:
	if line[0] == ">": # indicates fasta header: collect header info
		header = line.strip()[1:].split('|')
		genenames.append(header[1])
		protnumbers.append(header[2])
		gID.append(header)
		flag = 1 # turn on data collector
	elif flag == 1:
		strings.append(line.strip())
		flag = 0 # turn off data collector

# make an array of coordinates for reciprocal comparisons
ar = nti(len(strings),1)
# actually run the comparisons over the entire length of the array
cr = naaa(wordcomp,0,ar)
# calculate the hierarchy given the pairwise distances provided.
C = sch.linkage(cr, method=clustermeth)
# turn the hierarchy into a tree object
T = sch.to_tree(C)
# translate the tree object to newick format
tree = sof_tree2newick(T)
# apply the labels
for n,s in enumerate(protnumbers): #use any desired label here.
	nwstr1 = '(' + str(n) + ':'
	nwstr2 = ',' + str(n) + ':'
	str1 = '(' + genenames[n] + '|' + s + ':'
	str2 = ',' + genenames[n] + '|' + s + ':'
	tree = tree.replace(nwstr1,str1)
	tree = tree.replace(nwstr2,str2)
# save the newick text to a file
newick = open("%s.txt" %(clusterevolv), "w")
newick.write(tree)
newick.close()


clustcoll = []
# define clusters given the threshold
L = sch.fcluster(C,t,criterion=clustercrit)
### I DO NOT UNDERSTAND THE FOLLOWING OPERATIONS ###
cID = list(L)
clustcoll.append(cID)
cldict = {}
for n,cluster in enumerate(cID):
	cldict[protnumbers[n]] = cluster

# save clusterdata as a database
orderfile = open("%.csv" %clusterorder, "w")
orderfile.write("Gene_stable_ID,Gene_name,Protein_stable_ID,cluster,sequence\n") #write header

for n,gene in enumerate(gID):
	orderfile.write("%s,%s,%s,%s,%s\n" %(gene[0],gene[1],gene[2],clustcoll[n],strings[n]))
orderfile.close()

# save clusterdata as EvolView-readable data
evolview = open("%s.txt" %clustercolours,"w")
evolview.write(" ## leaf background color\n\n")
colours = config.getColour(max(clustcoll[t]))
countclust = 0
whichclust = []
# for each gene
for n,gene in enumerate(gID):
	# get the cluster and the assigned colour
	cluster = clustcoll[n]
	clr = colours[cluster-1]
	# get the gene name
	gn = genenames[n] + '|' + protnumbers[n]
	# write to file
	evolview.write("%s\t%s\tprefix\n" %(gn,clr))
	if cluster not in whichclust:
		whichclust.append(cluster)
		countclust += 1
evolview.close()
print "Found %s (%s) clusters." %(max(clustcoll),countclust)

