import scipy, pylab, csv
import scipy.cluster.hierarchy as sch
from collections import Counter

dbfolder = "/home/barbara/Dropbox/zinc_finger_data/data"
inputfile = csv.reader(open("%s/results/motifhits_dmel-count.csv" %(dbfolder))) 

'''
Read the input file and extract:
- matrix of observations (GOmatrix)
- labels of GOterms or motifs
- labels of genes (both numbers and names)
provides options to transpose the matrix
'''
observations = [] #will collect the actual matrix of observations
genenumbers = [] #will collect gene numbers (column 1)
genenames = [] #will collect gene names
for line in inputfile:
        if line[0] == "Gene_stable_ID": # gather GOterm/motif labels
            termlabels = line[3:] # skipping the first 3 columns
            termlabels = termlabels[:-1] # ignoring the last column (which is '' due to the trailing comma)
            continue
	genenumbers.append(line[0])
	genenames.append(line[1])
        GOline = [] #collects the observations per line
        for i in range(3,len(line)-1): #line ends in comma, so this removes the last empty item, and removes the first identifiers
            GOline.append(float(line[i]))
        observations.append(GOline)
GOmatrix = scipy.array(observations) # matrixify the observations
#GOmatrix = scipy.transpose(GOmatrix) # flips the matrix --- comment out if necessary!

'''
Calculate the clusters. Uses a cutoff value to define the number
of clustered categories that will be made.
'''
cutoff = 25 # determine cutoff: number of categories to be formed(1.1547)
L = sch.fclusterdata(GOmatrix, cutoff, criterion='maxclust', method='weighted') 
S = set(L) #turns the clustering into a set so as to remove duplicates
Llist = list(L) #turns the clustering into a list, so it may be indexed


print len(S)
print Counter(Llist)
'''
for i in range(len(S)):
	if Llist.count(i+1) > 10:
		print i+1, Llist.count(i+1)

'''
