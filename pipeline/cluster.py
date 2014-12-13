'''
This script can be used to create clusters of a fasta file with
strings (without regard for the exact meaning of the individual elements!).
It works through calculating Levenshtein distances, and therefore
clusters longer strings more efficiently.
It returns a csv file with several cluster options based on the cluster
threshold.
The script requires the following python modules:
- scipy
- pylab
- jellyfish
- numpy
- matplotlib
- itertools
- re
- csv
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 15 October 2014
'''

import scipy, pylab, re, itertools, csv, os
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
from jellyfish import levenshtein_distance as jld
from numpy import triu_indices as nti
from numpy import apply_along_axis as naaa


##### INPUT SPECIFICATIONS: CUSTOMIZE HERE! #####

#options for clustering:
clustermeth = "weighted"
threshold = [5,10,20,50,100,200,500]
clustercrit = "maxclust"

#input/output files and folders:
species = "dmel"
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data/results"
motiffile = "%s/motifseq_%s.fasta" %(dbfolder,species) #the file used for the clustering
infile = "%s/motifhits_%s.csv" %(dbfolder,species) #the file to apply the clustering to, and split into new files
clusterorder = "%s/clustering-string_%s-%s.csv" %(dbfolder,species,clustermeth) #this file will be made: contains group data

#### END OF CUSTOMIZATION! PLEASE DON'T EDIT BELOW #####


'''
Calculate pairwise distances for all the strings collected. As strings
are regular expressions, first expand the re. and then calculate all distances
pairwise. Return the minimal distance.
'''
# set up the measurement of levenshtein distance of a combination of strings
def wordcomp(coord):
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

'''
Read the input fasta file; extract:
- string [with motif data]
- ID labels [customized to separate on | character] + genenames and protein numbers separately
'''
inputfile = open(motiffile)
strings,gID,genenames,protnumbers = [],[],[],[] #will collect the clusterable data (1) + header information (2-4)
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
sch.dendrogram(C,labels=genenames)
plt.savefig("%s/dendrogram.png" %dbfolder)

'''
Collect data from the file that needs to be sorted into clusters.
'''
# collect infile to memory
infileread = csv.reader(open(infile))
ftc = []
for n,line in enumerate(infileread):
	if n == 0:
		ftchead = line
	else:
		ftc.append(line)
# find the column with protein ID (= unique identifying information)
pID = ftchead.index("Protein_stable_ID")


'''
interpret the clustering and apply it to a file with data (e.g. visualization,
GO terms, etc), so that these are clustered similarly.
Turned into a function so it can be repeated with different thresholds.
'''

def sortdata(thresh,nclust,cldict,outfolder):
	global pID,infile,ftchead
	for n in range(1,nclust+1):
		clusterfile = open("%s/%s_cluster%s.csv" %(outfolder,infile.split('/')[-1][:-4],n), "w")
		#write header
		lcollect = ""
		for item in ftchead:
			lcollect += item + ','
		clusterfile.write("%s\n" %lcollect[:-1])
		# write content
		for line in ftc:
			if cldict[line[pID]] == n:
				lcollect = ""
				for item in line:
					lcollect += item + ','
				clusterfile.write("%s\n" %lcollect[:-1])
		clusterfile.close()

clustcoll = []
for t in threshold:
	# define clusters given the threshold
	L = sch.fcluster(C,t,criterion=clustercrit)
	# translate the cluster into a dictionary
	cID = list(L)
	clustcoll.append(cID)
	cldict = {}
	for n,cluster in enumerate(cID):
		cldict[protnumbers[n]] = cluster
	# make a folder for these clusters
	outfolder = "%s/%s-clusters" %(dbfolder,t)
	if not os.path.exists(outfolder):
		os.system("mkdir %s" %(outfolder))
	sortdata(t,max(cID),cldict,outfolder)

# save clusterdata separately
orderfile = open(clusterorder, "w")
orderfile.write("Gene_stable_ID,Gene_name,Protein_stable_ID,")
tline = ""
for t in threshold:
	tline += str(t) + '_clusters,'
orderfile.write("%s\n" %tline[:-1])

for n,gene in enumerate(gID):
	orderfile.write("%s,%s,%s," %(gene[0],gene[1],gene[2]))
	ccline = ""
	for t,thresh in enumerate(threshold):
		ccline += str(clustcoll[t][n]) + ','
	orderfile.write("%s\n" %ccline[:-1])
orderfile.close()
