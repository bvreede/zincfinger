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
- collections
- jellyfish
- numpy
- matplotlib
- itertools
- re
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 15 October 2014
'''

import scipy, pylab, re, itertools
import scipy.cluster.hierarchy as sch
from collections import Counter
from jellyfish import levenshtein_distance as jld
from numpy import triu_indices as nti
from numpy import apply_along_axis as naaa
from matplotlib.pyplot import show

##### INPUT SPECIFICATIONS: CUSTOMIZE HERE! #####

#options for clustering:
clustermeth = "weighted"
threshold = 50
clustercrit = "maxclust"

#input/output files and folders:
species = "dmel"
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data"
motiffile = "%s/results/motifseq_%s.fasta" %(dbfolder,species) #the file used for the clustering
infile = "%s/results/motifhits_%s.csv" %(dbfolder,species) #the file to apply the clustering to, and split into new files
clusterorder = "%s/results/clustering-string_%s-%s.csv" %(dbfolder,species,clustermeth) #this file will be made: contains group data




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
# define clusters given the threshold
L = sch.fcluster(C,threshold,criterion=clustercrit)

# translate the cluster into a dictionary
cID = list(L)
cldict = {}
for n,cluster in enumerate(cID):
	cldict[protnumbers[n]] = cluster

# save clusterdata separately
orderfile = open(clusterorder, "w")
orderfile.write("Gene_stable_ID,Gene_name,Protein_stable_ID,Cluster\n")
for n,gene in enumerate(gID):
	orderfile.write("%s,%s,%s,%s\n" %(gID[0],gID[1],gID[2],cID[n]))
orderfile.close()

'''
Applying the clusters to data: open file that needs to be ordered according to the
clusters, and split it into however many clusters exist.
'''


for n in range(1,max(cID)+1):
	for p,gene in enumerate(gID):
		if cldict[gene[2]] == n:
			print n, gene[1], strings[p]

#print cldict

'''
# make a matrix of the result file
# then transpose the matrix and print to results
# this way, any number of headers and clusterdata can be added without specifying a number.

#L = sch.fcluster(C,cutoff,criterion='maxclust')
L8 = sch.fcluster(C,800,criterion='maxclust')
L9 = sch.fcluster(C,300,criterion='maxclust')
L0 = sch.fcluster(C,200,criterion='maxclust')
L1 = sch.fcluster(C,150,criterion='maxclust')
L2 = sch.fcluster(C,100,criterion='maxclust')
L3 = sch.fcluster(C,50,criterion='maxclust')
L4 = sch.fcluster(C,20,criterion='maxclust')
S = set(L0) #turns the clustering into a set so as to remove duplicates

#Llist = list(L) #turns the clustering into a list, so it may be indexed
Llist8 = list(L8)
Llist9 = list(L9)
Llist0 = list(L0)
Llist1 = list(L1)
Llist2 = list(L2)
Llist3 = list(L3)
Llist4 = list(L4)
'''

'''
Reporting on the results:

print "Number of categories: ", len(S)
print "Instances per category: ", Counter(Llist0) #counts instances per category
#print "Number of genenumbers/genenames: ", len(genenumbers), len(genenames)
print "Number of strings for clustering: ", len(strings)
print "Number of items assigned categories: ", len(L0), len(Llist0)

#cldict = {}
for i,ID in enumerate(gID):
	orderfile.write("%s,%s,%s,%s,%s,%s,%s,%s\n" %(ID,strings[i],Llist9[i],Llist0[i],Llist1[i],Llist2[i],Llist3[i],Llist4[i]))
	#cldict[genenumbers[i]] = Llist0[i]
orderfile.close()

sch.dendrogram(C,labels=genenames)
show()
'''
'''
#for each category, go into GO and motif files, and save the genes seperately
for c in range(1,cutoff+1):
outfile = infile[:-4] + "-cluster" + str(c) + '.csv'
selection = open(outfile, "w") #resultsfile (file that is cluster-specific)
to_select = open(infile) #file to read
for line in to_select:
l = line.strip().split(',')
if l[0] == "Gene_stable_ID":
selection.write(line)
elif l[0] in cldict:
if cldict[l[0]] == c:
selection.write(line)
selection.close()
to_select.close()
'''
