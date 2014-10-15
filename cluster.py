import scipy, pylab
import scipy.cluster.hierarchy as sch
from collections import Counter
from jellyfish import levenshtein_distance as jeljar
from numpy import triu_indices as nti
from numpy import apply_along_axis as naaa
from matplotlib.pyplot import show

species = "dmel"
dbfolder = "/home/barbara/Dropbox/shared_work/zinc_finger_data/data"
motiffile = "%s/results/motifseq_%s.fasta" %(dbfolder,species) #the file used for the clustering
infile = "%s/databases/GO_%s_old.csv" %(dbfolder,species) #the file to apply the clustering to, and split into new files
inputfile = open(motiffile)
clustermeth = "weighted"
clusterorder = "%s/results/clustering-string_%s-%s.csv" %(dbfolder,species,clustermeth)
cutoff = 150

orderfile = open(clusterorder, "w")

'''
Read the input file and extract:
- string with motif data (sequence of motifs in the protein)
- labels of genes (both numbers and names)
'''
observations = [] #will collect the actual matrix of observations
genenumbers = [] #will collect gene numbers (column 1)
genenames = [] #will collect gene names

for line in inputfile:
	if line[0] == ">":
		header = line.strip()[1:].split('|')
		genenumbers.append(header[0])
		genenames.append(header[1])
		flag = 1
	elif flag == 1:
		observations.append(line.strip())
		flag = 0

# set up the measurement of jaro distance of a combination of words
def wordcomp(coord):
	i,j = coord
	return jeljar(observations[i],observations[j])

# make the array of coordinates for reciprocal comparisons
ar = nti(len(observations),1)

# actually run the comparisons over the entire length of the array
matrix = naaa(wordcomp,0,ar)

'''
Calculate the clusters. Uses a cutoff value to define the number
of clustered categories that will be made.
'''
C = sch.linkage(matrix, method=clustermeth)
L = sch.fcluster(C,cutoff,criterion='maxclust')
sch.dendrogram(C,labels=observations)
show()

S = set(L) #turns the clustering into a set so as to remove duplicates
Llist = list(L) #turns the clustering into a list, so it may be indexed

'''
Reporting on the results:
'''
print "Number of categories: ", len(S)
print "Instances per category: ", Counter(Llist) #counts instances per category
print "Number of genenumbers/genenames: ", len(genenumbers), len(genenames)
print "Number of strings for clustering: ", len(observations)
print "Number of items assigned categories: ", len(L), len(Llist)

cldict = {}
for i in range(len(genenumbers)):
	orderfile.write("%s,%s,%s,%s\n" %(genenumbers[i],genenames[i],observations[i],Llist[i]))
	cldict[genenumbers[i]] = Llist[i]
orderfile.close()

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
